# Run main locally
import sys
from pathlib import Path  # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

from algos.base import CompressionAlgorithm
from algos.predictors.predictor_base import PredictorFunction
from algos.predictors.simple_predictor import SimplePredictor
from algos.predictors.prev_delta_predictor import PrevDeltaPredictor
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import math
import numpy as np
from shapely import GeometryType as GT
import bisect
from bitarray import bitarray, util, bits2bytes

#####                 #####
# --- FP-DELTA BASELINE ---#
#####                 #####

# char: 8 bits
# float: 32 bits
# double: 64 bits
# int: 32 bits
# long: 64 bits

D_CNT_SIZE = 16
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16

MAX_NUM_DELTAS = 10  # Max number of deltas in a chunk before split


class FpdPredictor(CompressionAlgorithm):
    offset = 0  # Used when parsing
    predictor = PrevDeltaPredictor()

# ---- HELPER METHODS

    def double_as_long(self, num):
        return struct.unpack('!q', struct.pack('!d', num))[0]

    def long_as_double(self, num):
        return struct.unpack('!d', struct.pack('!q', num))[0]

    def double_to_bytes(self, x):
        return struct.pack("!d", x)

    # Inline refactorized from https://github.com/ilanschnell/bitarray/blob/master/bitarray/util.py

    def uint_to_ba(self, x, length):
        if x == 0:
            return util.zeros(length or 1, "big")

        a = bitarray(0, "big")

        a.frombytes(x.to_bytes(bits2bytes(x.bit_length()), byteorder="big"))
        la = len(a)
        if la == length:
            return a

        return a[-length:] if la > length else util.zeros(length - la, "big") + a

    def uchar_to_bytes(self, x):
        return x.to_bytes(1, 'big')

    def zz_encode(self, num):
        return -2 * num - 1 if num < 0 else 2 * num

    def zz_decode(self, num):
        return - num // 2 if num & 1 == 1 else num // 2

    def get_zz_encoded_delta(self, prev_coords, curr_coord):
        #Uses normal fp-delta encoding if no predictor has been set
        pred_coord = prev_coords[-1] if self.predictor == None else self.predictor.predict(prev_coords)
        prev_coords.append(curr_coord)
        return self.zz_encode(self.double_as_long(curr_coord) - self.double_as_long(pred_coord))

    def deltas_fit_in_bits(self, d_x, d_y, max_bits):
        return (d_x == 0 or math.log2(d_x) < max_bits) and (d_y == 0 or math.log2(d_y) < max_bits)

    # Returns the number of coords in each ring for a polygon. If multipolygon, index on polygon first
    def point_count(self, geometry):
        coords = shapely.get_coordinates(geometry)

        if shapely.get_type_id(geometry) == GT.LINESTRING:
            return [0, len(coords)]

        ring_count = deque([])
        if shapely.get_type_id(geometry) == GT.POLYGON:
            ring_count.append(len(geometry.exterior.coords))
            for i in range(len(geometry.interiors)):
                ring_count.append(len(geometry.interiors[i].coords))

        elif shapely.get_type_id(geometry) == GT.MULTIPOLYGON:
            for polygon in list(geometry.geoms):
                poly_ring_count = deque()
                poly_ring_count.append(len(list(polygon.exterior.coords)))
                for i in range(len(polygon.interiors)):
                    poly_ring_count.append(len(polygon.interiors[i].coords))
                ring_count.append(poly_ring_count)
        return ring_count

    def append_header(self, bits, geometry, d_size):
        # Meta data
        bits.frombytes(self.uchar_to_bytes(d_size))
        bits.frombytes(self.uchar_to_bytes(int(shapely.get_type_id(geometry))))  # 1 byte is enough for storing type

        # Bounding Box
        # print(shapely.bounds(geometry))

    def decode_header(self, bin):
        delta_size, type = struct.unpack_from('!BB', bin)
        type = GT(type)
        self.offset += 2 * 8  # Offset is 2 bytes for BB
        return delta_size, type

    def append_delta_pair(self, bits, d_x_zig, d_y_zig, d_size):
        x_bytes = self.uint_to_ba(d_x_zig, d_size)
        y_bytes = self.uint_to_ba(d_y_zig, d_size)
        bits.extend(x_bytes)
        bits.extend(y_bytes)

    def fp_delta_encoding(self, geometry, d_size):
        # List of resulting bytes.
        bits = bitarray(endian='big')
        # Init with 'd_size', 'geom_type'
        self.append_header(bits, geometry, d_size)

        # Type specific variables
        geo_type = shapely.get_type_id(geometry)
        is_linestring = geo_type == GT.LINESTRING
        is_multipolygon = geo_type == GT.MULTIPOLYGON
        is_polygon = geo_type == GT.POLYGON

        # Fetches number of points in each ring, nestled for multipoly
        poly_buffer = self.point_count(geometry)
        ring_buffer = poly_buffer if is_polygon else []  # Not nestled for poly, else overwritten below

        prev_xs, prev_ys = [0], [0]  # Absolute value of previous coord
        chk_deltas = 0  # Cnt of 'deltas in chunk'
        chk_deltas_idx = 0  # Pointer to 'deltas of chunk'
        num_chks_ring = 0  # Cnt of 'number of chunks for current ring'
        num_chks_ring_idx = 0  # Pointer to latest 'number of chunks for ring'
        rem_points_ring = 0  # Cnt of 'points left to process in current ring'

        # Loop all coordinates
        for x, y in shapely.get_coordinates(geometry):
            if not is_linestring and rem_points_ring == 1:  # Is the whole ring processed? We skip last coordinate
                # Store number of chunks used for the ring
                rem_points_ring = 0
                bits[num_chks_ring_idx:num_chks_ring_idx + RING_CHK_CNT_SIZE] = self.uint_to_ba(num_chks_ring, RING_CHK_CNT_SIZE)
                num_chks_ring = 0
                continue  # Skip last coordinate
            d_x_zig = self.get_zz_encoded_delta(prev_xs, x)  # Calculated delta based on previous iteration
            d_y_zig = self.get_zz_encoded_delta(prev_ys, y)
        
            # ---- CREATE NEW CHUNK? If 'first chunk', 'delta doesn't fit', 'new ring', or 'reached max deltas'
            if chk_deltas_idx == 0 or not self.deltas_fit_in_bits(d_x_zig, d_y_zig, d_size) or rem_points_ring == 0 or chk_deltas == MAX_NUM_DELTAS:
                # If not 'first chunk' -> save previous chunk's size
                if chk_deltas_idx != 0:
                    bits[chk_deltas_idx:chk_deltas_idx + D_CNT_SIZE] = self.uint_to_ba(chk_deltas, D_CNT_SIZE)

                ###### ---- INITIALIZE NEW CHUNK ----- ######
                chk_deltas = 0

                ### __ RING/MULTI-POLYGON META-DATA __ ###
                if not is_linestring:
                    # Ran out of points -> fetch number of points in next ring
                    if rem_points_ring == 0:
                        # Check if we ran out of rings -> fetch rings of NEXT POLYGON
                        if is_multipolygon and len(ring_buffer) == 0:
                            ring_buffer = poly_buffer.popleft()
                            bits.extend(self.uint_to_ba(len(ring_buffer), POLY_RING_CNT_SIZE))  # Append 'nbr of rings in poly'
                        # Set 'remaining points' to cnt in new ring
                        rem_points_ring = ring_buffer.popleft()
                        num_chks_ring_idx = len(bits)
                        bits.extend(self.uint_to_ba(0, RING_CHK_CNT_SIZE))  # Reserve space for number of chunks for current ring
                        num_chks_ring = 1
                    else:
                        num_chks_ring += 1
                ### __ ------------ END ------------- __ ###

                # Preparing chunk size (number of deltas)
                chk_deltas_idx = len(bits)
                bits.extend(self.uint_to_ba(0, D_CNT_SIZE))  # Reserve space for size

                # Add full coordinates
                bits.frombytes(self.double_to_bytes(x))
                bits.frombytes(self.double_to_bytes(y))
                prev_xs, prev_ys = [x], [y] # Reset points used to predict future points. Makes predictions chunk isolated

            else:
                # Delta fits, append it
                self.append_delta_pair(bits, d_x_zig, d_y_zig, d_size)

                chk_deltas += 1

            # Coord has been processed, remove it
            rem_points_ring -= 1

        # All points processed. Update size of final chunk
        bits[chk_deltas_idx:chk_deltas_idx + D_CNT_SIZE] = self.uint_to_ba(chk_deltas, D_CNT_SIZE)

        # util.pprint(bits)
        # print([int.from_bytes(i, 'big') for i in bytes], '\n')
        return bits.tobytes()

    def bytes_to_decoded_coord(self, bin, prev_coords, input_size=64):
        bin = bin[self.offset: self.offset + input_size]
        val = util.ba2int(bin, signed=False)
        pred_coord = prev_coords[-1] if self.predictor == None else self.predictor.predict(prev_coords)  #Uses normal fp-delta encoding if no predictor has been set
        val = self.zz_decode(val) + self.double_as_long(pred_coord)
        val = self.long_as_double(val)
        prev_coords.append(val)
        self.offset += input_size
        return val

    def bytes_to_double(self, bin):
        bin = bin[self.offset:self.offset + 64]
        val = struct.unpack('!d', bin)[0]
        self.offset += 8 * 8
        return val

    def bytes_to_uint(self, bin, len):
        val = util.ba2int(bin[self.offset:self.offset + len], signed=False)
        self.offset += len
        return val

    def sequence_decoder(self, bin, seq_list, delta_size):
        chk_size = self.bytes_to_uint(bin, D_CNT_SIZE)
        # Extract reset point
        x = self.bytes_to_double(bin)
        y = self.bytes_to_double(bin)
        prev_xs = [x]
        prev_ys = [y]
        seq_list.append((x, y))
        # Loop through deltas in chunk
        for _ in range(chk_size):
            x = self.bytes_to_decoded_coord(bin, prev_xs, delta_size)
            y = self.bytes_to_decoded_coord(bin, prev_ys, delta_size)
            seq_list.append((x, y))
          
    def ring_decoder(self, bin, polygon_list, delta_size):
        # Extract number of chunks for a ring
        chks_in_ring = self.bytes_to_uint(bin, RING_CHK_CNT_SIZE)
        ring_coords = []
        # Loop through chunks in ring
        for i in range(chks_in_ring):
            self.sequence_decoder(bin, ring_coords, delta_size)
        polygon_list.append(ring_coords)

    def polygon_decoder(self, bin, multipolygon_coords, delta_size):
        # Extract number of rings for a polygon
        rings_in_poly = self.bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        polygon_coords = []
        # Loop through rings in polygon
        for _ in range(rings_in_poly):
            self.ring_decoder(bin, polygon_coords, delta_size)
        multipolygon_coords.append(shapely.Polygon(shell=polygon_coords[0], holes=polygon_coords[1:]))

    def fp_delta_decoding(self, bin_in):
        self.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.decode_header(bin)

        binary_length = len(bin)
        coords = []
        if type == GT.LINESTRING:
            while (self.offset + 8 < binary_length):  # While != EOF
                self.sequence_decoder(bin, coords, delta_size)
            geometry = shapely.LineString(coords)

        elif type == GT.POLYGON:
            while (self.offset + 8 < binary_length):  # While != EOF, i.e. at least one byte left
                self.ring_decoder(bin, coords, delta_size)
            geometry = shapely.Polygon(shell=coords[0], holes=coords[1:])

        elif type == GT.MULTIPOLYGON:
            while (self.offset + 8 < binary_length):  # While != EOF
                self.polygon_decoder(bin, coords, delta_size)
            geometry = shapely.MultiPolygon(coords)
        return geometry

    def calculate_delta_size(self, geometry, return_deltas = False):
        deltas = []
        RESET_POINT_SIZE = 64 * 2 + D_CNT_SIZE
        coords = shapely.get_coordinates(geometry)
        prev = [[0], [0]] #prev[0] is x-coordinates and prev[1] is y-coordinates. Format as matrix to satisfy self.get_zz_encoded_delta
        bit_cnts = {}
        for coord in coords:
            bit_cnt = 0
            for i in range(2):
                if return_deltas:
                    deltas.append(prev[i][-1] if self.predictor == None else self.predictor.predict(prev[i]))

                d = self.get_zz_encoded_delta(prev[i], coord[i])
                d_bit_cnt = 1 if d == 0 else math.ceil(math.log2(d))
                bit_cnt = max(bit_cnt, d_bit_cnt)
                
            if bit_cnt not in bit_cnts:
                bit_cnts[bit_cnt] = 1
            else:
                bit_cnts[bit_cnt] += 1
        bit_cnts = dict(sorted(bit_cnts.items(), reverse=True))

        tot_size = {}
        upper_cnt = 0
        lower_cnt = len(coords)
        for n in bit_cnts.keys():
            tot_size[n] = n * lower_cnt * 2 + RESET_POINT_SIZE * upper_cnt
            lower_cnt -= bit_cnts[n]
            upper_cnt += bit_cnts[n]

        # print({i: tot_size[i] // 8 for i in tot_size.keys()})
        return min(tot_size, key=tot_size.get), bit_cnts, deltas

    def compress(self, geometry):
        s = time.perf_counter()
        optimal_size, _ = self.calculate_delta_size(geometry)
        bin = self.fp_delta_encoding(geometry, optimal_size)
        t = time.perf_counter()
        return t - s, bin

    def decompress(self, bin):
        s = time.perf_counter()
        geometry = self.fp_delta_decoding(bin)
        t = time.perf_counter()
        return t - s, geometry


# ---- UNARY ---- #

    def vertices(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords

    def type(self, bin):
        s = time.perf_counter()
        type = struct.unpack_from('!B', bin, offset=1)[0]  # 1 Byte offset
        if type == GT.LINESTRING:
            type = 'LineString'
        elif type == GT.POLYGON:
            type = 'Polygon'
        elif type == GT.MULTIPOLYGON:
            type = 'MultiPolygon'
        t = time.perf_counter()
        return t - s, type

    def bounding_box(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds

    def add_vertex(self, args):
        bin, insert_idx, pos = args
        s = time.perf_counter()

        _, geometry = self.decompress(bin)
        ragged = shapely.to_ragged_array([geometry])
        points = ragged[1]
        is_linestring = shapely.get_type_id(geometry) == shapely.GeometryType.LINESTRING

        # Use binary search O(log n) to find the index of the first element greater than insert_idx
        end_idx = min(len(ragged[2][0]) - 1, bisect.bisect_right(ragged[2][0], insert_idx))

        # Is the first coordinate in the ring?
        if ragged[2][0][end_idx - 1] == insert_idx and not is_linestring:
                points = np.delete(points, ragged[2][0][end_idx] - 1, axis=0)
        else:
            for i in range(end_idx, len(ragged[2][0])):
                    ragged[2][0][i] += 1
                    
        points = np.insert(points, insert_idx, pos, axis=0)        
 
        geometry = shapely.from_ragged_array(geometry_type=shapely.get_type_id(geometry), coords=points, offsets=ragged[2])[0]
        _, bin = self.compress(geometry)
        
        t = time.perf_counter()
        return t - s, bin

    def is_intersecting(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        _, r_geo = self.decompress(r_bin)
        res = shapely.intersects(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, res

    def intersection(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        _, r_geo = self.decompress(r_bin)
        res = shapely.intersection(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, res


def main():
    import random
    import json
    import pandas as pd
    import tqdm
    from shapely.geometry import shape
    from alg_fpd import Fpd
    x = FpdPredictor()

    geom1 = shapely.wkt.loads("MULTIPOLYGON (((13.193709 55.7021381, 13.1937743 55.7021279, 13.1938355 55.7021184, 13.1938461 55.702109, 13.1938566 55.7020984, 13.1938611 55.7020902, 13.1938655 55.7020774, 13.1938655 55.7020633, 13.1938583 55.7020408, 13.1938402 55.7020014, 13.1937184 55.7017259, 13.1937008 55.7016876, 13.1936836 55.7016654, 13.1936537 55.7016428, 13.1936223 55.7016242, 13.1935741 55.7016036, 13.1935354 55.7015911, 13.1935006 55.701584, 13.1934829 55.701598, 13.1934673 55.7016115, 13.1934736 55.7016164, 13.1934776 55.7016216, 13.1934875 55.7016633, 13.1934985 55.7016898, 13.1935196 55.7017337, 13.1935659 55.7018353, 13.1936162 55.7018282, 13.1936551 55.7019155, 13.1936651 55.7019377, 13.1936955 55.7020047, 13.1936497 55.7020119, 13.193709 55.7021381)), ((13.1938175 55.7017126, 13.1938602 55.7017068, 13.1939048 55.7017007, 13.1938998 55.7016861, 13.193892 55.7016685, 13.1938831 55.7016589, 13.193871 55.701651, 13.1938602 55.701646, 13.1938405 55.7016438, 13.193822 55.7016456, 13.1938062 55.7016517, 13.1937985 55.7016571, 13.1937953 55.7016646, 13.1937979 55.7016746, 13.1938017 55.7016836, 13.1938052 55.7016908, 13.1938175 55.7017126)), ((13.1940245 55.7019788, 13.19398 55.7019848, 13.1939372 55.7019907, 13.1939585 55.7020383, 13.1939692 55.7020479, 13.1939841 55.7020512, 13.1939975 55.7020519, 13.1940079 55.702051, 13.1940198 55.7020497, 13.1940317 55.7020463, 13.1940395 55.7020422, 13.1940435 55.7020369, 13.1940452 55.7020314, 13.1940457 55.7020218, 13.1940245 55.7019788)), ((13.1939779 55.7015541, 13.1939529 55.701555, 13.1939622 55.7015658, 13.1939755 55.7015942, 13.194075 55.7018201, 13.1941382 55.7019637, 13.1941483 55.7019866, 13.194164 55.7020087, 13.1941899 55.7020304, 13.1942142 55.7020424, 13.1942291 55.7020486, 13.1942638 55.702042, 13.195019 55.7018988, 13.1948681 55.7018923, 13.1944181 55.7018687, 13.1944172 55.7018717, 13.194395 55.7018706, 13.1942164 55.7018622, 13.194172 55.7017564, 13.1941218 55.701761, 13.1941279 55.7017262, 13.1941357 55.7016818, 13.1940872 55.7015737, 13.1940769 55.7015503, 13.1939779 55.7015541), (13.1942341 55.7020059, 13.1942075 55.7020095, 13.1941895 55.7019673, 13.1941696 55.701921, 13.1941936 55.7019177, 13.1941884 55.7019055, 13.19426 55.7018958, 13.1942645 55.7019063, 13.1943172 55.7018991, 13.1943567 55.7019912, 13.1943041 55.7019984, 13.1943086 55.7020089, 13.1942394 55.7020183, 13.1942341 55.7020059)))")
    geom2 = shapely.wkt.loads(
        "LINESTRING (13.199378 55.7034667, 13.1999441 55.7033986, 13.200125 55.7033882, 13.2002723 55.7033936, 13.2004383 55.7034097, 13.2005935 55.7034211, 13.2007699 55.703423, 13.2011275 55.7034136, 13.2012413 55.7034103, 13.2012947 55.7034088)")
    geom3 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848603 55.7056619, 13.1846238 55.7056422, 13.1846085 55.7057159, 13.1846356 55.7057179, 13.1848537 55.7057363), (13.1846694 55.705714, 13.1846543 55.7057128, 13.1846563 55.705705, 13.1846714 55.7057062, 13.1846694 55.705714), (13.1847425 55.7057123, 13.1847405 55.7057201, 13.1847254 55.7057188, 13.1847274 55.705711, 13.1847425 55.7057123), (13.1848001 55.7057179, 13.1848152 55.7057192, 13.1848131 55.705727, 13.1847981 55.7057258, 13.1848001 55.7057179), (13.1848068 55.7056929, 13.1848088 55.7056851, 13.1848239 55.7056863, 13.1848218 55.7056941, 13.1848068 55.7056929), (13.1847507 55.7056878, 13.1847356 55.7056865, 13.1847377 55.7056787, 13.1847528 55.70568, 13.1847507 55.7056878), (13.1846811 55.7056732, 13.184679 55.705681, 13.184664 55.7056798, 13.184666 55.705672, 13.1846811 55.7056732))')
    geom4 = shapely.wkt.loads(
        'POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.184812 55.705646, 13.1848537 55.7057363))')
    geom5 = shapely.wkt.loads('MULTIPOLYGON (((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848861 55.705646, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1847300 55.705712, 13.1847425 55.7057123)), ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848841 55.705626, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1848861 55.705646, 13.1847425 55.7057123)))')
    geom_list = [geom4, geom1, geom2, geom3, geom5]
    for geom in geom_list:
       t, bin3 = x.compress(geom)
       #t, bin3 = x.add_vertex((bin3, 2, [0.2, 0.2]))
       t, decomp = x.decompress(bin3)
       print(decomp == geom)

    # faulty_add = shapely.wkt.loads(
    #     'POLYGON ((13.2288713 55.713951, 13.2287496 55.7139225, 13.2286938 55.7139984, 13.2287396 55.7140091, 13.228704 55.7140575, 13.2287798 55.7140753, 13.2288713 55.713951))')
    # faulty_add = geom2
    # t, bin3 = x.compress(faulty_add)
    # # binary = bitarray()
    # # binary.frombytes(bin3)
    # # util.pprint(binary)
    # t, bin3 = x.add_vertex((bin3, 5, (13.2, 55.7)))
    # t, decomp = x.decompress(bin3)
    # print(decomp == faulty_add)
    # # print(x.bytes_to_float32(b'0xa70x300x530x41'))
    # # print(bin3.hex(sep=' '))
    # # print("-DELTASIZE-_TY_POINTSINCHK_XFIRST-----------------_YFIRST-----------------")
    # # print("-DELTASIZE-_TY_RINGS------_POINTSINCHK_XFIRST-----------------_YFIRST-----------------_--XD1---------_--YD1---------_--XD2---------_--YD2---------_--XD3---------_--YD3---------_--XD4---------_--YD4---------_--XD5---------_--YD5---------_--XD6---------_--YD6---------_--XD7---------_--YD7---------_--XD8---------_--YD8---------_--XD9---------_--YD9---------")
    # # t = shapely.to_ragged_array([geom1])
    # # print(t[2])

    # # print(x.bytes_to_float(bin3, size=16))
    # # print(x.calculate_delta_size(geom1))
    # # t, bin = x.decompress(bin)

    # # print(bin(0b01100111 | 0b0000000000000000))
    # # t = '{0:0{1}b}'.format(0b101, 10)
    # # print(t)
    # # print(int(t, 2))
    # # print(x.get_poly_ring_count(geom5))

    # # print(len(bytes(shapely.to_wkt(geom4), 'utf-8')), shapely.to_wkt(geom4))
    # print(shapely.to_wkt(faulty_add, rounding_precision=-1))
    # print(shapely.to_wkt(decomp, rounding_precision=-1))
    # print("LEN COMP", len(bin3))

    # return

    # DATASET_PATH = "data/lund_building_highway.json"
    # # DATASET_PATH = "data/world.json"
    # NBR_ITER = 16000

    # SEED = 123  # If we want to enforce the same ordering and indexes for multiple runs, else None
    # random.seed(SEED)  # Init random

    # # Extract the nested feature attribute of the geo_json file containing the geometries
    # with open(DATASET_PATH, 'r') as f:
    #     data = json.loads(f.read())
    # file_df: pd.DataFrame = pd.json_normalize(data, record_path=['features'])
    # # Create a dataframe suitable for the WKT format for easy convertion to shapely objects
    # df = pd.DataFrame(
    #     {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

    # max_idx = len(df) - 1
    # unary_idxs = [random.randint(0, max_idx) for i in range(NBR_ITER)]  # Generate list of indexes to query on
    # random.seed(SEED)  # Reset random

    # # Compress files, benchmark unaries
    # for idx in tqdm.tqdm(unary_idxs):  # List of single idxs
    #     t, comp = x.compress(shape(df.iloc[idx]))
    #     t, decomp = x.decompress(comp)
    #     if (shape(df.iloc[idx]) != decomp):
    #         print("FAIL")


if __name__ == "__main__":
    main()
