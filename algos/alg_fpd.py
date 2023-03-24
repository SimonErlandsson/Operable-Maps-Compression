# Run main locally
import sys
from pathlib import Path # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

from algos.base import CompressionAlgorithm
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
#--- FP-DELTA BASELINE ---#
#####                 #####

## char: 8 bits
## float: 32 bits
## double: 64 bits
## int: 32 bits
## long: 64 bits

D_CNT_SIZE = 6
POLY_RING_CNT_SIZE = 9
RING_CHK_CNT_SIZE = 11

class Fpd(CompressionAlgorithm):
    MAX_NUM_DELTAS = 31
    offset = 0 # Used when parsing

# ---- HELPER METHODS
    
    def double_as_long(self, num):
        return struct.unpack('!q', struct.pack('!d', num))[0]
      
    def long_as_double(self, num):
        return struct.unpack('!d', struct.pack('!q', num))[0]
    
    def double_to_bytes(self, x):
        return struct.pack("!d", x)
    
     #Inline refactorized from https://github.com/ilanschnell/bitarray/blob/master/bitarray/util.py
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
        bits.frombytes(self.uchar_to_bytes(int(shapely.get_type_id(geometry)))) # 1 byte is enough for storing type
    
    def decode_header(self, bin):
        delta_size, type = struct.unpack_from('!BB', bin)
        type = GT(type)
        self.offset += 2 * 8 # Offset is 2 bytes for BB  
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
        geotype = shapely.get_type_id(geometry)
        is_linestring = geotype == GT.LINESTRING  
        is_multipolygon = geotype == GT.MULTIPOLYGON
        is_polygon = geotype == GT.POLYGON
        
        # Fetches number of points in each ring, nestled for multipoly
        poly_buffer = self.point_count(geometry)
        ring_buffer = poly_buffer if is_polygon else [] # Not nestled for poly, else overwritten below

        prev_x, prev_y = 0, 0 # Absolute value of previous coord
        chk_deltas = 0 # Cnt of 'deltas in chunk'
        chk_deltas_idx = 0  # Pointer to 'deltas of chunk'
        num_chks_ring = 0 # Cnt of 'number of chunks for current ring'
        num_chks_ring_idx = 0 # Pointer to latest 'number of chunks for ring'
        rem_points_ring = 0 # Cnt of 'points left to process in current ring'
        
        # Loop all coordinates
        for x, y in shapely.get_coordinates(geometry):
            d_x_zig = self.get_zz_encoded_delta(prev_x, x) # Calculated delta based on previous iteration
            d_y_zig = self.get_zz_encoded_delta(prev_y, y)
            prev_x, prev_y = (x, y)
            
            # ---- CREATE NEW CHUNK? If 'first chunk', 'delta doesn't fit', 'new ring', or 'reached max deltas'
            if chk_deltas_idx == 0 or not self.deltas_fit_in_bits(d_x_zig, d_y_zig, d_size) or rem_points_ring == 0 or chk_deltas == self.MAX_NUM_DELTAS:
                # If not 'first chunk' -> save previous chunk's size
                if chk_deltas_idx != 0:
                    bits[chk_deltas_idx:chk_deltas_idx+D_CNT_SIZE] = self.uint_to_ba(chk_deltas, D_CNT_SIZE)

                ###### ---- INITIALIZE NEW CHUNK ----- ######
                chk_deltas = 0

                ### __ RING/MULTI-POLYGON META-DATA __ ###
                if not is_linestring:
                    # Ran out of points -> fetch number of points in next ring
                    if rem_points_ring == 0:
                        # Check if we ran out of rings -> fetch rings of NEXT POLYGON
                        if is_multipolygon and len(ring_buffer) == 0:
                            ring_buffer = poly_buffer.popleft()
                            bits.extend(self.uint_to_ba(len(ring_buffer), POLY_RING_CNT_SIZE)) # Append 'nbr of rings in poly'
                        # Set 'remaining points' to cnt in new ring
                        rem_points_ring = ring_buffer.popleft()
                        num_chks_ring_idx = len(bits)
                        bits.extend(self.uint_to_ba(0, RING_CHK_CNT_SIZE)) # Reserve space for number of chunks for current ring
                        num_chks_ring = 1
                    else:
                        num_chks_ring += 1
                ### __ ------------ END ------------- __ ###
                    
                # Preparing chunk size (number of deltas)
                chk_deltas_idx = len(bits)
                bits.extend(self.uint_to_ba(0, D_CNT_SIZE)) # Reserve space for size

                # Add full coordinates
                bits.frombytes(self.double_to_bytes(x))
                bits.frombytes(self.double_to_bytes(y))
            else:
                # Delta fits, append it
                self.append_delta_pair(bits, d_x_zig, d_y_zig, d_size)
                chk_deltas += 1
                    
            # Coord has been processed, remove it
            rem_points_ring -= 1
            if not is_linestring and rem_points_ring == 0: # Is the whole ring processed?
                # Store number of chunks used for the ring
                bits[num_chks_ring_idx:num_chks_ring_idx+RING_CHK_CNT_SIZE] = self.uint_to_ba(num_chks_ring, RING_CHK_CNT_SIZE)
                num_chks_ring = 0

        # All points processed. Update size of final chunk
        bits[chk_deltas_idx:chk_deltas_idx+D_CNT_SIZE] = self.uint_to_ba(chk_deltas, D_CNT_SIZE)

        #util.pprint(bits)
        #print([int.from_bytes(i, 'big') for i in bytes], '\n')
        return bits.tobytes()
    
    def bytes_to_decoded_coord(self, bin, prev_coord, input_size=64):
        bin = bin[self.offset : self.offset + input_size]
        val = util.ba2int(bin, signed=False)
        val = self.zz_decode(val) + self.double_as_long(prev_coord)
        val = self.long_as_double(val)
        self.offset += input_size
        return val
    
    def bytes_to_double(self, bin):
        bin = bin[self.offset:self.offset+64]
        val = struct.unpack('!d', bin)[0]
        self.offset += 8 * 8
        return val
    
    def bytes_to_uint(self, bin, len):
        val = util.ba2int(bin[self.offset:self.offset+len], signed=False)
        self.offset += len
        return val
    
    def sequence_decoder(self, bin, seq_list, delta_size):
        chk_size = self.bytes_to_uint(bin, D_CNT_SIZE)
        # Extract reset point
        x = self.bytes_to_double(bin)
        y = self.bytes_to_double(bin)
        seq_list.append((x, y))
        # Loop through deltas in chunk
        for _ in range(chk_size):
            x = self.bytes_to_decoded_coord(bin, x, delta_size)
            y = self.bytes_to_decoded_coord(bin, y, delta_size)
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
            self.ring_decoder(bin,polygon_coords,delta_size)
        multipolygon_coords.append(shapely.Polygon(shell=polygon_coords[0], holes=polygon_coords[1:]))
        
    def fp_delta_decoding(self, bin_in):
        self.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.decode_header(bin)

        binary_length = len(bin)
        coords = []
        if type == GT.LINESTRING:
            while (self.offset + 8 < binary_length): # While != EOF  
                self.sequence_decoder(bin, coords, delta_size)
            geometry = shapely.LineString(coords)

        elif type == GT.POLYGON:
            while (self.offset + 8 < binary_length): # While != EOF, i.e. at least one byte left
                self.ring_decoder(bin, coords, delta_size)
            geometry = shapely.Polygon(shell=coords[0], holes=coords[1:])
            
        elif type == GT.MULTIPOLYGON:
            while (self.offset + 8 < binary_length): # While != EOF
                self.polygon_decoder(bin, coords, delta_size)
            geometry = shapely.MultiPolygon(coords)
        return geometry
    
    def calculate_delta_size(self, geometry, return_deltas = False):
        deltas = [[], []]
        RESET_POINT_SIZE = 64 * 2 + D_CNT_SIZE
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

        #print({i: tot_size[i] // 8 for i in tot_size.keys()})
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


# ---- UNARY ---- #

    def vertices(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords

    def type(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        type = geometry.geom_type
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

    geom4 = shapely.wkt.loads('MULTIPOLYGON (((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848861 55.705646, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1847300 55.705712, 13.1847425 55.7057123)), ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848841 55.705626, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1848861 55.705646, 13.1847425 55.7057123)))')
    geom4 = shapely.wkt.loads('POLYGON ((13.1969457 55.6906953, 13.197008 55.6906875, 13.1971208 55.6906735, 13.1971913 55.6906647, 13.1971772 55.6906286, 13.1971097 55.690637, 13.1970943 55.6905975, 13.1969816 55.6906115, 13.1969958 55.6906479, 13.1969304 55.6906561, 13.1969457 55.6906953))')
    
    fpd = Fpd()
    _, org = fpd.compress(geom4)
    _, org = fpd.add_vertex((org, 0, (13.19693570000000, 55.6907853)))
    _, de = fpd.decompress(org)   
    print(shapely.to_wkt(de, rounding_precision=-1))
    print('POLYGON ((13.196936 55.690785, 13.196946 55.690695, 13.197008 55.690688, 13.197121 55.690674, 13.197191 55.690665, 13.197177 55.690629, 13.19711 55.690637, 13.197094 55.690598, 13.196982 55.690612, 13.196996 55.690648, 13.19693 55.690656, 13.196946 55.690695, 13.196936 55.690785))')
    print(shapely.to_wkt(geom4, rounding_precision=-1))

    x = Fpd()
    return

    DATASET_PATH = "data/lund_building_highway.json"
    #DATASET_PATH = "data/world.json"
    NBR_ITER = 16000


    SEED = 123 # If we want to enforce the same ordering and indexes for multiple runs, else None
    random.seed(SEED) # Init random

    # Extract the nested feature attribute of the geo_json file containing the geometries
    with open(DATASET_PATH, 'r') as f:
        data = json.loads(f.read())
    file_df: pd.DataFrame = pd.json_normalize(data, record_path=['features'])
    # Create a dataframe suitable for the WKT format for easy convertion to shapely objects
    df = pd.DataFrame(
        {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

    max_idx = len(df) - 1
    unary_idxs = [random.randint(0, max_idx) for i in range(NBR_ITER)] # Generate list of indexes to query on
    random.seed(SEED) # Reset random

    # Compress files, benchmark unaries
    for idx in tqdm.tqdm(unary_idxs): # List of single idxs
        t, comp = x.compress(shape(df.iloc[idx]))
        t, decomp = x.decompress(comp)
        if(shape(df.iloc[idx]) != decomp):
            print("FAIL")
        


    


if __name__ == "__main__":
    main()
