# Run main locally
import sys
from pathlib import Path  # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

from algos.base import CompressionAlgorithm
from algos.fpd_extended_lib.intersection_chunk_bbox_wrapper import *
from algos.fpd_extended_lib.add_vertex import AddVertex
from algos.fpd_extended_lib.functions import *

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
from var_float import VarFloat


#####                 #####
# --- FP-DELTA BASELINE ---#
#####                 #####

# char: 8 bits
# float: 32 bits
# double: 64 bits
# int: 32 bits
# long: 64 bits


class FpdExtended(CompressionAlgorithm):
    USE_DEFAULT_DOUBLE = False
    FLOAT_SIZE = 40
    EXPONENT = 6
    D_CNT_SIZE = 16
    POLY_RING_CNT_SIZE = 16
    RING_CHK_CNT_SIZE = 16
    MAX_NUM_DELTAS = 15  # Max number of deltas in a chunk before split
    EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing

    offset = 0  # Used when parsing
# ---- HELPER METHODS

    var_float = VarFloat(EXPONENT, FLOAT_SIZE)

    def double_as_long(self, num):
        if self.USE_DEFAULT_DOUBLE:
            return struct.unpack('!q', struct.pack('!d', num))[0]
        else:
            return self.var_float.bits_to_long(self.var_float.float_to_bin(num))

    def long_as_double(self, num):
        if self.USE_DEFAULT_DOUBLE:
            return struct.unpack('!d', struct.pack('!q', num))[0]
        else:
            return self.var_float.bin_to_float(self.var_float.long_to_bits(num))
        
    def double_to_bytes(self, x):
        if self.USE_DEFAULT_DOUBLE:
            return struct.pack("!d", x)
        else:
            return self.var_float.float_to_bin(x)
    
    def bin_to_double(self, bin):
        if self.USE_DEFAULT_DOUBLE:
            return struct.unpack('!d', bin)[0]
        else:
            return self.var_float.bin_to_float(bin)

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
        return - (num + 1) // 2 if num & 1 == 1 else num // 2

    def get_zz_encoded_delta(self, prev_coord, curr_coord):
        return self.zz_encode(self.double_as_long(curr_coord) - self.double_as_long(prev_coord))

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
        bounds = shapely.bounds(geometry)
        bits.frombytes(self.double_to_bytes(bounds[0]))
        bits.frombytes(self.double_to_bytes(bounds[1]))
        bits.frombytes(self.double_to_bytes(bounds[2]))
        bits.frombytes(self.double_to_bytes(bounds[3]))

        append_intersection_header(self, bits, geometry)
   
    def decode_header(self, bin, get_idxs=False):
        delta_size, type = struct.unpack_from('!BB', bin)
        type = GT(type)
        self.offset += 2 * 8 + 4 * self.FLOAT_SIZE  # Offset is 2 bytes for BB + 64 * 4 for bounding box
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

        prev_x, prev_y = 0, 0  # Absolute value of previous coord
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
                bits[num_chks_ring_idx:num_chks_ring_idx + self.RING_CHK_CNT_SIZE] = self.uint_to_ba(num_chks_ring, self.RING_CHK_CNT_SIZE)
                num_chks_ring = 0
                continue  # Skip last coordinate
            d_x_zig = self.get_zz_encoded_delta(prev_x, x)  # Calculated delta based on previous iteration
            d_y_zig = self.get_zz_encoded_delta(prev_y, y)
            prev_x, prev_y = (x, y)

            # ---- CREATE NEW CHUNK? If 'first chunk', 'delta doesn't fit', 'new ring', or 'reached max deltas'
            if chk_deltas_idx == 0 or not self.deltas_fit_in_bits(d_x_zig, d_y_zig, d_size) or rem_points_ring == 0 or chk_deltas == self.MAX_NUM_DELTAS:
                # If not 'first chunk' -> save previous chunk's size
                if chk_deltas_idx != 0:
                    bits[chk_deltas_idx:chk_deltas_idx + self.D_CNT_SIZE] = self.uint_to_ba(chk_deltas, self.D_CNT_SIZE)

                ###### ---- INITIALIZE NEW CHUNK ----- ######
                chk_deltas = 0

                ### __ RING/MULTI-POLYGON META-DATA __ ###
                if not is_linestring:
                    # Ran out of points -> fetch number of points in next ring
                    if rem_points_ring == 0:
                        # Check if we ran out of rings -> fetch rings of NEXT POLYGON
                        if is_multipolygon and len(ring_buffer) == 0:
                            ring_buffer = poly_buffer.popleft()
                            bits.extend(self.uint_to_ba(len(ring_buffer), self.POLY_RING_CNT_SIZE))  # Append 'nbr of rings in poly'
                        # Set 'remaining points' to cnt in new ring
                        rem_points_ring = ring_buffer.popleft()
                        num_chks_ring_idx = len(bits)
                        bits.extend(self.uint_to_ba(0, self.RING_CHK_CNT_SIZE))  # Reserve space for number of chunks for current ring
                        num_chks_ring = 1
                    else:
                        num_chks_ring += 1
                ### __ ------------ END ------------- __ ###

                # Preparing chunk size (number of deltas)
                chk_deltas_idx = len(bits)
                bits.extend(self.uint_to_ba(0, self.D_CNT_SIZE))  # Reserve space for size

                # Add full coordinates
                bits.frombytes(self.double_to_bytes(x))
                bits.frombytes(self.double_to_bytes(y))
            else:
                # Delta fits, append it
                self.append_delta_pair(bits, d_x_zig, d_y_zig, d_size)
                chk_deltas += 1

            # Coord has been processed, remove it
            rem_points_ring -= 1

        # All points processed. Update size of final chunk
        bits[chk_deltas_idx:chk_deltas_idx + self.D_CNT_SIZE] = self.uint_to_ba(chk_deltas, self.D_CNT_SIZE)

        # util.pprint(bits)
        # print([int.from_bytes(i, 'big') for i in bytes], '\n')
        return bits.tobytes()

    def bytes_to_decoded_coord(self, bin, prev_coord, inputsize=32):
        bin = bin[self.offset: self.offset + inputsize]
        val = util.ba2int(bin, signed=False)
        val = self.zz_decode(val) + self.double_as_long(prev_coord)
        val = self.long_as_double(val)
        self.offset += inputsize
        return val

    def bytes_to_double(self, bin, offset = None):
        bin = bin[self.offset:self.offset + self.FLOAT_SIZE]
        val = self.bin_to_double(bin)
        self.offset += self.FLOAT_SIZE
        return val

    def bytes_to_uint(self, bin, len):
        val = util.ba2int(bin[self.offset:self.offset + len], signed=False)
        self.offset += len
        return val

    def sequence_decoder(self, bin, seq_list, delta_size):
        chk_size = self.bytes_to_uint(bin, self.D_CNT_SIZE)
        # Extract reset point
        x = self.bytes_to_double(bin)
        y = self.bytes_to_double(bin)
        seq_list.append((x, y))
        # Loop through deltas in chunk
        for _ in range(chk_size):
            x = self.bytes_to_decoded_coord(bin, x, delta_size)
            y = self.bytes_to_decoded_coord(bin, y, delta_size)

            seq_list.append((float("{:.7f}".format(x)), float("{:.7f}".format(y))))

    def ring_decoder(self, bin, polygon_list, delta_size):
        # Extract number of chunks for a ring
        chks_in_ring = self.bytes_to_uint(bin, self.RING_CHK_CNT_SIZE)
        ring_coords = []
        # Loop through chunks in ring
        for i in range(chks_in_ring):
            self.sequence_decoder(bin, ring_coords, delta_size)
        polygon_list.append(ring_coords)

    def polygon_decoder(self, bin, multipolygon_coords, delta_size):
        # Extract number of rings for a polygon
        rings_in_poly = self.bytes_to_uint(bin, self.POLY_RING_CNT_SIZE)
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
            while (self.offset + self.EOF_THRESHOLD <= binary_length):  # While != EOF
                self.sequence_decoder(bin, coords, delta_size)
            geometry = shapely.LineString(coords)

        elif type == GT.POLYGON:
            while (self.offset + self.EOF_THRESHOLD <= binary_length):  # While != EOF, i.e. at least one byte left
                self.ring_decoder(bin, coords, delta_size)
                
            geometry = shapely.Polygon(shell=coords[0], holes=coords[1:])

        elif type == GT.MULTIPOLYGON:
            while (self.offset + self.EOF_THRESHOLD <= binary_length):  # While != EOF
                self.polygon_decoder(bin, coords, delta_size)
            geometry = shapely.MultiPolygon(coords)
        return geometry

    def calculate_delta_size(self, geometry, return_deltas=False):
        deltas = [[], []]
        RESET_POINT_SIZE = self.FLOAT_SIZE * 2 + self.D_CNT_SIZE
        coords = shapely.get_coordinates(geometry)
        prev = [0, 0]
        bit_cnts = {}
        for coord in coords:
            bit_cnt = 0
            for i in range(2):
                d = self.get_zz_encoded_delta(prev[i], coord[i])
                d_bit_cnt = 1 if d == 0 else math.ceil(math.log2(d))
                bit_cnt = max(bit_cnt, d_bit_cnt)
                if return_deltas:
                    deltas[0].append(coord[i] - prev[i])
                    deltas[1].append(d)

            if bit_cnt not in bit_cnts:
                bit_cnts[bit_cnt] = 1
            else:
                bit_cnts[bit_cnt] += 1
            prev = coord
        bit_cnts = dict(sorted(bit_cnts.items(), reverse=True))

        tot_size = {}
        upper_cnt = 0
        lower_cnt = len(coords)
        for n in bit_cnts.keys():
            tot_size[n] = n * lower_cnt * 2 + RESET_POINT_SIZE * upper_cnt
            lower_cnt -= bit_cnts[n]
            upper_cnt += bit_cnts[n]

        return min(tot_size, key=tot_size.get), bit_cnts, deltas

    def compress(self, geometry):
        s = time.perf_counter()
        optimal_size, _, _ = self.calculate_delta_size(geometry)
        bin = self.fp_delta_encoding(geometry, optimal_size)
        t = time.perf_counter()
        return t - s, bin

    def decompress(self, bin):
        s = time.perf_counter()
        geometry = self.fp_delta_decoding(bin)
        t = time.perf_counter()
        return t - s, geometry

    # Export helper functions
    get_chunks = get_chunks
# ---- UNARY ---- #

    def vertices(self, bin):
        return Funcs(self).vertices(bin)

    def type(self, bin):
        return Funcs(self).type(bin)

    def bounding_box(self, bin):
        return Funcs(self).bounding_box(bin)

    # TODO: Handle bounding box
    def add_vertex(self, args):
        return AddVertex(self).add_vertex(args)
 
# ---- BINARY ---- #

    def is_intersecting(self, args):
        return Intersection(self).is_intersecting(args)

    def intersection(self, args):
        return Intersection(self).intersection(args)


def main():
    import random
    import json
    import pandas as pd
    import tqdm
    from shapely.geometry import shape
    from alg_fpd import Fpd
    x = FpdExtended()
    geom1 = shapely.wkt.loads('POLYGON ((-24.3 10.48, -19.32 12.44, -15.3 14.2, -15.3 13.78, -15.3 13.9, -15.06 10.4, -17.44 11.38, -19.18 11.46, -14.82 9.08, -12.9 10.14, -12.08 7.86, -14.36 5.94, -15.92 8.34, -16.86 3.48, -19.38 4.4, -18.2 6.52, -20.08 7.4, -24.34 6.68, -24.24 8.66, -27.52 11.1,  -27.0 11.1, -24.3 10.48))')


    t, bin = x.compress(geom1)
    t, geomx = x.decompress(bin)
    print(x.bin2float(x.float2bin(-180.12223)))
    


if __name__ == "__main__":
    main()

