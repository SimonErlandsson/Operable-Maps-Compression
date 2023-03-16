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

D_CNT_SIZE = 16
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16

MAX_NUM_DELTAS = 2 # Max number of deltas in a chunk before split

class FpdExtended(CompressionAlgorithm):
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
        return - num // 2 if num & 1 == 1 else num // 2

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

        # Bounding Box
        #(shapely.bounds(geometry))
    
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
        # TODO: REMOVE!!!
        d_size = 50
        # List of resulting bytes.
        bits = bitarray(endian='big')
        # Init with 'd_size', 'geom_type'
        self.append_header(bits, geometry, d_size)

        # Type specific variables
        is_linestring = shapely.get_type_id(geometry) == GT.LINESTRING  
        is_multipolygon = shapely.get_type_id(geometry) == GT.MULTIPOLYGON
        is_polygon = shapely.get_type_id(geometry) == GT.POLYGON
        
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
            if chk_deltas_idx == 0 or not self.deltas_fit_in_bits(d_x_zig, d_y_zig, d_size) or rem_points_ring == 0 or chk_deltas == MAX_NUM_DELTAS:
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

    def calculate_delta_size(self, geometry):
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
        return min(tot_size, key=tot_size.get)

    def compress(self, geometry):
        s = time.perf_counter()
        optimal_size = self.calculate_delta_size(geometry)
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
        type = struct.unpack_from('!B', bin, offset=1)[0] # 1 Byte offset
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
        bin_in, insert_idx, pos = args
        s = time.perf_counter()

        self.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.decode_header(bin)
        # Type specific variables
        is_linestring = type == GT.LINESTRING  
        is_multipolygon = type == GT.MULTIPOLYGON
        is_polygon = type == GT.POLYGON

        p_idx = 0
        chunks_in_ring_left = 0 # Used for iteration
        chunks_in_ring = 0 # Cache for store later
        chunks_in_ring_idx = 0
        rings_left = 0
        while (p_idx <= insert_idx):
            if is_multipolygon and rings_left == 0:
                rings_left = self.bytes_to_uint(bin, POLY_RING_CNT_SIZE)
            if not is_linestring and chunks_in_ring_left == 0:
                chunks_in_ring_offset = self.offset
                chunks_in_ring_left = self.bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            deltas_in_chunk_offset = self.offset
            deltas_in_chunk = self.bytes_to_uint(bin, D_CNT_SIZE)
            #print(p_idx, deltas_in_chunk, insert_idx)
            # Found chunk to append/prepend?
            if p_idx <= insert_idx and insert_idx <= p_idx + deltas_in_chunk:
                # Is first (reset)?
                if p_idx == insert_idx:
                    # Does the delta to the old reset point fit? -> The added point becomes the reset
                    reset_point = (self.bytes_to_double(bin), self.bytes_to_double(bin))
                    d_x = self.get_zz_encoded_delta(pos[0], reset_point[0])
                    d_y = self.get_zz_encoded_delta(pos[1], reset_point[1])
                    if self.deltas_fit_in_bits(d_x, d_y, delta_size):
                        left = bin[0:deltas_in_chunk_offset]
                        middle = bitarray(endian='big')
                        
                        middle.extend(self.uint_to_ba(deltas_in_chunk + 1, D_CNT_SIZE))
                        # Add full coordinates
                        middle.frombytes(self.double_to_bytes(pos[0]))
                        middle.frombytes(self.double_to_bytes(pos[1]))
                        self.append_delta_pair(middle, d_x, d_y, delta_size)

                        right = bin[deltas_in_chunk_offset + D_CNT_SIZE + 64 * 2:]
                        bin = left
                        bin.extend(middle)
                        bin.extend(right)
                    # Does not fit, create new chunk using the old chunk's counter
                    else:
                        left = bin[0:deltas_in_chunk_offset]
                        middle = bitarray(endian='big')
                        
                        middle.extend(self.uint_to_ba(0, D_CNT_SIZE))
                        # Add full coordinates
                        middle.frombytes(self.double_to_bytes(pos[0]))
                        middle.frombytes(self.double_to_bytes(pos[1]))
                        right = bin[deltas_in_chunk_offset + D_CNT_SIZE:]
                        bin = left
                        bin.extend(middle)
                        bin.extend(right)
                        bin[chunks_in_ring_offset:chunks_in_ring_offset+RING_CHK_CNT_SIZE] = self.uint_to_ba(chunks_in_ring + 1, RING_CHK_CNT_SIZE)
                else:
                    print("") 
                break
            else:
                # Jump to next chunk
                p_idx += 1 + deltas_in_chunk
                self.offset += 64 * 2 + delta_size * 2 * deltas_in_chunk
                chunks_in_ring_left -= 1

            # A B C D
              
              # 1 Första, delta får plats: blir nya start av gamla chunken, bildar kedja efter
              # 2 Första, delta får inte plats. Ny chunk prependas
              # 3 Mitten: delta får plats, lägg til direkt, men ändra den efter
              # 4 Mitten: delta får inte plats, dela gamla chunken i två nya
              # 5 Sista: delta får plats, lägg till direkt
              # 6 Sista: delta får inte plats, skapa ny chunk och appenda efter

              # 2, 6 Samma skiljer bara offset där chunk ska läggas till
              # 3, 5 Samma förutom att 3 behöver ändra den efter
              # 
              # Bör fixa så inte start-koordinat finns i både början och slut
    
        bin = bin.tobytes()
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
    
    #x_bytes = util.int2ba(d_x_zig, d_size, 'big', signed=False)


def main():
    import random
    import json
    import pandas as pd
    import tqdm
    from shapely.geometry import shape

    x = FpdExtended()
    

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
