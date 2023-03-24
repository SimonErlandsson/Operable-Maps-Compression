# Run main locally
import sys
from pathlib import Path  # if you haven't already done so
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

EOF_THRESHOLD = D_CNT_SIZE + 64 * 2 # Number of bits required to continue parsing


class FpdExtended(CompressionAlgorithm):
    offset = 0  # Used when parsing

# ---- HELPER METHODS

    def get_bounds_intersect(self,  bound_min1, bound_min2, bound_max1, bound_max2):
        if bound_min1 <= bound_min2 and bound_max1 <= bound_max2 and bound_max1 >= bound_min2:
            return (bound_min2, bound_max1)
        elif bound_min1 >= bound_min2 and bound_max1 <= bound_max2:
            return (bound_min1, bound_max1)
        if bound_min2 <= bound_min1 and bound_max2 <= bound_max1 and bound_max2 >= bound_min1:
            return (bound_min1, bound_max2)
        elif bound_min2 >= bound_min1 and bound_max2 <= bound_max1:
            return (bound_min2, bound_max2)
        else:
            return (None, None)

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

        coords = shapely.get_coordinates(geometry)
        bits.extend(self.uint_to_ba(int(len(coords)), 4 * 8)) #size of integer
        sorted_idxs = [np.argsort([coord[0] for coord in coords]), np.argsort([coord[1] for coord in coords])]
        idx_bits = math.ceil(math.log2(len(coords)))
        for i in range(2):
            for idx in range(len(coords)):
                bits.extend(self.uint_to_ba(int(sorted_idxs[i][idx]), idx_bits))

        
    def decode_header(self, bin, get_idxs = False):
        delta_size, type = struct.unpack_from('!BB', bin)
        type = GT(type)
        self.offset += 2 * 8 + 4 * 64  # Offset is 2 bytes for BB + 64 * 4 for bounding box
        
        #Code segment needed for extracting sorted indexes
        coord_count = self.bytes_to_uint(bin, 4 * 8)
        idx_sizes = math.ceil(math.log2(coord_count))
        sorted_idxs = [[],[]]
        if get_idxs:
            for i in range(2):
                for _ in range(coord_count):
                    sorted_idxs[i].append(self.bytes_to_uint(bin, idx_sizes))
            return delta_size, type, sorted_idxs
        else:
            self.offset += idx_sizes * 2 * coord_count
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
                bits[num_chks_ring_idx:num_chks_ring_idx + RING_CHK_CNT_SIZE] = self.uint_to_ba(num_chks_ring, RING_CHK_CNT_SIZE)
                num_chks_ring = 0
                continue  # Skip last coordinate
            d_x_zig = self.get_zz_encoded_delta(prev_x, x)  # Calculated delta based on previous iteration
            d_y_zig = self.get_zz_encoded_delta(prev_y, y)
            prev_x, prev_y = (x, y)

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

    def bytes_to_decoded_coord(self, bin, prev_coord, input_size=64):
        bin = bin[self.offset: self.offset + input_size]
        val = util.ba2int(bin, signed=False)
        val = self.zz_decode(val) + self.double_as_long(prev_coord)
        val = self.long_as_double(val)
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
            while (self.offset + EOF_THRESHOLD <= binary_length):  # While != EOF
                self.sequence_decoder(bin, coords, delta_size)
            geometry = shapely.LineString(coords)

        elif type == GT.POLYGON:
            while (self.offset + EOF_THRESHOLD <= binary_length):  # While != EOF, i.e. at least one byte left
                self.ring_decoder(bin, coords, delta_size)
            geometry = shapely.Polygon(shell=coords[0], holes=coords[1:])

        elif type == GT.MULTIPOLYGON:
            while (self.offset + EOF_THRESHOLD <= binary_length):  # While != EOF
                self.polygon_decoder(bin, coords, delta_size)
            geometry = shapely.MultiPolygon(coords)
        return geometry

    def calculate_delta_size(self, geometry, return_deltas = False):
        deltas = []
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
                    deltas.append(coord[i] - prev[i])

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

        # print({i: tot_size[i] // 8 for i in tot_size.keys()})
        return min(tot_size, key=tot_size.get), bit_cnt, deltas

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

    def vertices(self, bin_in):
        s = time.perf_counter()
        
        coords = []
        self.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.decode_header(bin)
        # Type specific variables
        is_linestring = type == GT.LINESTRING
        is_multipolygon = type == GT.MULTIPOLYGON

        chunks_in_ring_left = 0  # Used for iteration
        chunks_in_ring = 0
        rings_left = 0
        bin_len = len(bin)
        while (self.offset + EOF_THRESHOLD <= bin_len):
            if is_multipolygon and rings_left == 0:
                rings_left = self.bytes_to_uint(bin, POLY_RING_CNT_SIZE)
            if not is_linestring and chunks_in_ring_left == 0:
                chunks_in_ring_left = self.bytes_to_uint(bin, RING_CHK_CNT_SIZE)
                chunks_in_ring = chunks_in_ring_left

            # Go through chunk (inlined sequence decode)
            deltas_in_chunk = self.bytes_to_uint(bin, D_CNT_SIZE)
            # Extract reset point
            x = self.bytes_to_double(bin)
            y = self.bytes_to_double(bin)
            if chunks_in_ring_left == chunks_in_ring:
                x_ring, y_ring = (x, y)
            coords.append([x, y])
            # Loop through deltas in chunk
            for _ in range(deltas_in_chunk):
                x = self.bytes_to_decoded_coord(bin, x, delta_size)
                y = self.bytes_to_decoded_coord(bin, y, delta_size)
                coords.append([x, y])
            chunks_in_ring_left -= 1
            if chunks_in_ring_left == 0:
                coords.append([x_ring, y_ring])
                rings_left -= 1

        coords = np.array(coords)

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
        
        bounds = list(struct.unpack_from('!dddd', bin, offset=2)) # Skip first part of header

        t = time.perf_counter()
        return t - s, bounds
    
    # Supply cache if used repeatedly for same shape
    def access_vertex(self, bin_in, access_idx, cache=None):
        old_offset = self.offset
        self.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.decode_header(bin)
        # Type specific variables
        is_linestring = type == GT.LINESTRING
        is_multipolygon = type == GT.MULTIPOLYGON

        p_idx = 0
        chunks_in_ring_left = 0  # Used for iteration
        rings_left = 0
        while (p_idx <= access_idx):
            if is_multipolygon and rings_left == 0:
                rings_left = self.bytes_to_uint(bin, POLY_RING_CNT_SIZE)
            if not is_linestring and chunks_in_ring_left == 0:
                chunks_in_ring_left = self.bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            deltas_in_chunk_offset = self.offset
            deltas_in_chunk = self.bytes_to_uint(bin, D_CNT_SIZE)
            
            # Found chunk containing vertex?
            if p_idx <= access_idx and access_idx <= p_idx + deltas_in_chunk:
                x, y, cache = self.access_vertex_chk(bin, deltas_in_chunk_offset, access_idx - p_idx, delta_size, cache)
                break
            else:
                # Jump to next chunk
                p_idx += 1 + deltas_in_chunk
                self.offset += 64 * 2 + delta_size * 2 * deltas_in_chunk
                chunks_in_ring_left -= 1
                if (chunks_in_ring_left == 0):
                    rings_left -= 1
        self.offset = old_offset
        return (x, y), cache
    
    # Supply the offset to D_CNT, and idx is the index within the chunk
    def access_vertex_chk(self, bin, chk_offset, idx, delta_size, cache=None):
        old_offset = self.offset
        self.offset = chk_offset + D_CNT_SIZE
        # Extract reset point
        x, y = (self.bytes_to_double(bin), self.bytes_to_double(bin))
        # Loop through deltas in chunk
        for _ in range(idx):
            x = self.bytes_to_decoded_coord(bin, x, delta_size)
            y = self.bytes_to_decoded_coord(bin, y, delta_size)
        self.offset = old_offset
        return (x, y), cache

    # TODO: Handle bounding box
    def add_vertex(self, args):
        def create_chunk(reset_point, delta_cnt=0):
                middle = bitarray(endian='big')
                middle.extend(self.uint_to_ba(delta_cnt, D_CNT_SIZE))
                # Add full coordinates
                middle.frombytes(self.double_to_bytes(reset_point[0]))
                middle.frombytes(self.double_to_bytes(reset_point[1]))
                return middle

        bin_in, insert_idx, pos = args
        s = time.perf_counter()

        self.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.decode_header(bin)
        # Type specific variables
        is_linestring = type == GT.LINESTRING
        is_multipolygon = type == GT.MULTIPOLYGON

        p_idx = 0
        chunks_in_ring_left = 0  # Used for iteration
        chunks_in_ring = 0  # Cache for store later
        rings_left = 0
        bin_len = len(bin)
        while (p_idx <= insert_idx):
            if is_multipolygon and rings_left == 0:
                rings_left = self.bytes_to_uint(bin, POLY_RING_CNT_SIZE)
            if not is_linestring and chunks_in_ring_left == 0:
                chunks_in_ring_offset = self.offset
                chunks_in_ring_left = self.bytes_to_uint(bin, RING_CHK_CNT_SIZE)
                chunks_in_ring = chunks_in_ring_left
            deltas_in_chunk_offset = self.offset
            deltas_in_chunk = self.bytes_to_uint(bin, D_CNT_SIZE)
            
            # print(p_idx, deltas_in_chunk, insert_idx)
            # Found chunk to append/prepend?
            if p_idx <= insert_idx and insert_idx <= p_idx + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0):
                deltas_left = max(insert_idx - p_idx - 1, 0)
                deltas_right = deltas_in_chunk - deltas_left
                #print(deltas_left, deltas_right)
                # Handle left
                if p_idx != insert_idx:  # Has a left coordinate?
                    split = deltas_in_chunk_offset + D_CNT_SIZE + 64 * 2 + delta_size * 2 * deltas_left
                    # Update delta cnt
                    bin[deltas_in_chunk_offset:deltas_in_chunk_offset + D_CNT_SIZE] = self.uint_to_ba(deltas_left, D_CNT_SIZE)
                else:
                    split = deltas_in_chunk_offset

                middle = create_chunk(pos)
                chunks_in_ring += 1

                # Handle chunk tail
                if deltas_right > 0 and p_idx != insert_idx:
                    # Get the absolute coordinate for first right coordinate
                    rst_p, _ = self.access_vertex_chk(bin, deltas_in_chunk_offset, insert_idx - p_idx, delta_size)
                    right = create_chunk(rst_p, deltas_right - 1)
                    # Append old tail, without the one extracted point
                    right.extend(bin[split + delta_size * 2:])
                    chunks_in_ring += 1
                else:
                    right = bin[split:]

                left = bin[0:split]
                if not is_linestring:
                    left[chunks_in_ring_offset:chunks_in_ring_offset + RING_CHK_CNT_SIZE] = self.uint_to_ba(chunks_in_ring, RING_CHK_CNT_SIZE)
                bin = left
                bin.extend(middle)
                bin.extend(right)

                break
            else:
                # Jump to next chunk
                p_idx += 1 + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0)
                self.offset += 64 * 2 + delta_size * 2 * deltas_in_chunk
                chunks_in_ring_left -= 1
                if (chunks_in_ring_left == 0):
                    rings_left -= 1

                if self.offset >= bin_len and is_linestring:
                    # Reached end without appending: is linestring!
                    new = create_chunk(pos)
                    bin.extend(new)
                    break

            # A B C D

            # 1 Första, delta får plats: blir nya start av gamla chunken, bildar kedja efter
            # - Skapa chunk: <left: counter r pos> counter s: 1 + r, s, delta s, r <right: counter r pos + D_SIZE_CHUNK + 2 * 64>
            # 2 Första, delta får inte plats. Ny chunk prependas
            # - Skapa chunk: <left: counter r pos> counter s (0), s <right: counter r pos>
            # - Uppdater count av antal chunks för ring
            # 3 Mitten: delta innan får plats, lägg til direkt, men ändra den efter
            # - Två deltas <left: counter r pos + D_SIZE_CHUNK + 64 * 2 + D_SIZE * n> s, d: s_(n+1) <right: left pos + D_SIZE>
            # - Uppdater CNT av antal deltas i chunk
            # 4 Mitten: delta får inte plats, dela gamla chunken i två nya
            # - Skapa chunk: <left: counter r pos+ D_SIZE_CHUNK + 64 * 2 + D_SIZE * n> counter s (0), s <right: counter r pos>
            # Måste även hantera om deltat efter inte får plats...
            # 5 Sista: delta får plats, lägg till direkt
            # 6 Sista: delta får inte plats, skapa ny chunk och appenda efter

            # 2, 6 Samma skiljer bara offset där chunk ska läggas till
            # 3, 5 Samma förutom att 3 behöver ändra den efter
            #
            # Bör fixa så inte start-koordinat finns i både början och slut

            # def add_chunk(bin, left_pos, right_pos, try_merge_left, try_merge_left)

            # MVP:
            # Lägg till en ny chunk som består av enbart den nya koordinaten. Delar den gamla i två nya chunks (fall 4)
            # Skapa bara ny del-chunk om inte slut på höger eller vänster
        bin = bin.tobytes()
        t = time.perf_counter()
        return t - s, bin

    def is_intersecting(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        l_bounds = list(struct.unpack_from('!dddd', l_bin, offset=2)) # Skip first part of header
        r_bounds = list(struct.unpack_from('!dddd', r_bin, offset=2)) # Skip first part of header
        if (self.get_bounds_intersect(l_bounds[0], r_bounds[0], l_bounds[2], r_bounds[2]) == (None, None) 
            or self.get_bounds_intersect(l_bounds[1], r_bounds[1], l_bounds[3], r_bounds[3]) == (None, None)):
            t = time.perf_counter()
            return t - s, False
        
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
    x = FpdExtended()

    geom1 = shapely.wkt.loads("MULTIPOLYGON (((13.193709 55.7021381, 13.1937743 55.7021279, 13.1938355 55.7021184, 13.1938461 55.702109, 13.1938566 55.7020984, 13.1938611 55.7020902, 13.1938655 55.7020774, 13.1938655 55.7020633, 13.1938583 55.7020408, 13.1938402 55.7020014, 13.1937184 55.7017259, 13.1937008 55.7016876, 13.1936836 55.7016654, 13.1936537 55.7016428, 13.1936223 55.7016242, 13.1935741 55.7016036, 13.1935354 55.7015911, 13.1935006 55.701584, 13.1934829 55.701598, 13.1934673 55.7016115, 13.1934736 55.7016164, 13.1934776 55.7016216, 13.1934875 55.7016633, 13.1934985 55.7016898, 13.1935196 55.7017337, 13.1935659 55.7018353, 13.1936162 55.7018282, 13.1936551 55.7019155, 13.1936651 55.7019377, 13.1936955 55.7020047, 13.1936497 55.7020119, 13.193709 55.7021381)), ((13.1938175 55.7017126, 13.1938602 55.7017068, 13.1939048 55.7017007, 13.1938998 55.7016861, 13.193892 55.7016685, 13.1938831 55.7016589, 13.193871 55.701651, 13.1938602 55.701646, 13.1938405 55.7016438, 13.193822 55.7016456, 13.1938062 55.7016517, 13.1937985 55.7016571, 13.1937953 55.7016646, 13.1937979 55.7016746, 13.1938017 55.7016836, 13.1938052 55.7016908, 13.1938175 55.7017126)), ((13.1940245 55.7019788, 13.19398 55.7019848, 13.1939372 55.7019907, 13.1939585 55.7020383, 13.1939692 55.7020479, 13.1939841 55.7020512, 13.1939975 55.7020519, 13.1940079 55.702051, 13.1940198 55.7020497, 13.1940317 55.7020463, 13.1940395 55.7020422, 13.1940435 55.7020369, 13.1940452 55.7020314, 13.1940457 55.7020218, 13.1940245 55.7019788)), ((13.1939779 55.7015541, 13.1939529 55.701555, 13.1939622 55.7015658, 13.1939755 55.7015942, 13.194075 55.7018201, 13.1941382 55.7019637, 13.1941483 55.7019866, 13.194164 55.7020087, 13.1941899 55.7020304, 13.1942142 55.7020424, 13.1942291 55.7020486, 13.1942638 55.702042, 13.195019 55.7018988, 13.1948681 55.7018923, 13.1944181 55.7018687, 13.1944172 55.7018717, 13.194395 55.7018706, 13.1942164 55.7018622, 13.194172 55.7017564, 13.1941218 55.701761, 13.1941279 55.7017262, 13.1941357 55.7016818, 13.1940872 55.7015737, 13.1940769 55.7015503, 13.1939779 55.7015541), (13.1942341 55.7020059, 13.1942075 55.7020095, 13.1941895 55.7019673, 13.1941696 55.701921, 13.1941936 55.7019177, 13.1941884 55.7019055, 13.19426 55.7018958, 13.1942645 55.7019063, 13.1943172 55.7018991, 13.1943567 55.7019912, 13.1943041 55.7019984, 13.1943086 55.7020089, 13.1942394 55.7020183, 13.1942341 55.7020059)))")
    geom2 = shapely.wkt.loads(
        "LINESTRING (13.199378 55.7034667, 13.1999441 55.7033986, 13.200125 55.7033882, 13.2002723 55.7033936, 13.2004383 55.7034097, 13.2005935 55.7034211, 13.2007699 55.703423, 13.2011275 55.7034136, 13.2012413 55.7034103, 13.2012947 55.7034088)")
    geom3 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848603 55.7056619, 13.1846238 55.7056422, 13.1846085 55.7057159, 13.1846356 55.7057179, 13.1848537 55.7057363), (13.1846694 55.705714, 13.1846543 55.7057128, 13.1846563 55.705705, 13.1846714 55.7057062, 13.1846694 55.705714), (13.1847425 55.7057123, 13.1847405 55.7057201, 13.1847254 55.7057188, 13.1847274 55.705711, 13.1847425 55.7057123), (13.1848001 55.7057179, 13.1848152 55.7057192, 13.1848131 55.705727, 13.1847981 55.7057258, 13.1848001 55.7057179), (13.1848068 55.7056929, 13.1848088 55.7056851, 13.1848239 55.7056863, 13.1848218 55.7056941, 13.1848068 55.7056929), (13.1847507 55.7056878, 13.1847356 55.7056865, 13.1847377 55.7056787, 13.1847528 55.70568, 13.1847507 55.7056878), (13.1846811 55.7056732, 13.184679 55.705681, 13.184664 55.7056798, 13.184666 55.705672, 13.1846811 55.7056732))')
    geom4 = shapely.wkt.loads(
        'POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.184812 55.705646, 13.1848537 55.7057363))')
    geom5 = shapely.wkt.loads('MULTIPOLYGON (((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848861 55.705646, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1847300 55.705712, 13.1847425 55.7057123)), ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848841 55.705626, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1848861 55.705646, 13.1847425 55.7057123)))')
    geom_list = [geom2, geom4]
    t, bin3 = x.compress(geom5)
    print(x.access_vertex(bin3, 3))
    t, bin4 = x.compress(geom4)
    #print(x.intersection((bin3, bin4)))
    #t, bin3 = x.add_vertex((bin3, 2, [0.2, 0.2]))
      
    t, decomp = x.decompress(bin3)

#     #    print(decomp == geom)

#     faulty_add = shapely.wkt.loads('POLYGON ((13.1948358 55.7202146, 13.1948222 55.7201619, 13.19485 55.7201605, 13.1950544 55.72015, 13.1950566 55.7201627, 13.1950679 55.7202269, 13.1950494 55.720227, 13.195037 55.7202706, 13.1950151 55.7202973, 13.1949752 55.7203157, 13.1949257 55.7203314, 13.1948712 55.7203365, 13.1948402 55.7203393, 13.1948 55.720343, 13.1947972 55.7203333, 13.1946787 55.7203443, 13.1945924 55.7203522, 13.1945976 55.7203702, 13.1944086 55.7203873, 13.1943905 55.7203184, 13.1944804 55.7203102, 13.1944776 55.7203004, 13.1945399 55.7202947, 13.19453 55.7202607, 13.1945931 55.7202549, 13.194588 55.7202375, 13.1948358 55.7202146))')
#     faulty_vertices = shapely.wkt.loads('POLYGON ((-7.989662638999988 10.161990662000065, -7.999842895999933 10.150699361000093, -8.007232625999961 10.13832285600003, -8.012348591999853 10.124861145000068, -8.015707560999942 10.110365906000126, -8.034879516999865 10.085793763000055, -8.057875529999905 10.073081360000103, -8.108776814999942 10.054193624000078, -8.127535359999854 10.037114563000031, -8.14774084499993 10.01003611300007, -8.15848954299986 9.984921367000055, -8.148877726999928 9.973888449000057, -8.136268676999919 9.969444275000043, -8.14009273299996 9.959703268000098, -8.152391723399319 9.949988098372046, -8.21339600502165 9.95651264114548, -8.261466525659614 9.97754349403698, -8.33958112057229 10.00458316088418, -8.390656047738446 10.05866249637711, -8.420700122912365 10.12776386900731, -8.429713345194727 10.178838796173466, -8.477783864933429 10.199869649065022, -8.57392490531015 10.217896093629804, -8.673070353114326 10.232918131666452, -8.724145280280482 10.283993058832607, -8.73015409513539 10.332063578571308, -8.748180540599492 10.392151728919146, -8.78423343062832 10.41919139666561, -8.85633921068603 10.434213433802938, -8.898400915569823 10.44923547094021, -8.946471435308524 10.503314806433139, -8.955484657590887 10.560398548454259, -8.96449788077257 10.665552811113344, -9.159784368054147 10.680574849149934, -9.213863702647757 10.704610108569625, -9.300991520742116 10.722636554033727, -9.373097300799827 10.749676220880872, -9.370092892472996 10.791737926663984, -9.364084077618088 10.824786408366037, -9.343053225625852 10.86684811414915, -9.313009150451933 10.899896595851203, -9.237898962966767 10.953975931344132, -9.150771145771785 10.984020006518051, -9.069652142532334 11.017068489119481, -8.976515510482443 11.080161045995453, -8.922436174989514 11.131235974060928, -8.88638328496063 11.230381420965784, -8.874365655250813 11.338540091052323, -8.853334803258576 11.428672316574136, -8.811350064154169 11.526719468228862, -8.78718014204992 11.557928777890766, -8.778635220139279 11.585516668759226, -8.767160610973406 11.61005280199521, -8.746815558591152 11.622870184411454, -8.726516886999889 11.649470927000024, -8.712874307999897 11.64006581600006, -8.690705118999858 11.568028870000092, -8.687346150999957 11.55001963300009, -8.679026245999893 11.529839986000113, -8.667967488999892 11.510874736000062, -8.656753702999936 11.496431173000133, -8.602183389999908 11.472091573000029, -8.57520829299989 11.470282898000036, -8.549576782999907 11.490126648000029, -8.543995727999913 11.477646790000037, -8.539424991999908 11.456312043000139, -8.535520792999904 11.43808848100008, -8.527562622999938 11.423386536000095, -8.51578039599994 11.418968201000084, -8.488133503999933 11.418193055000074, -8.427723754999903 11.400080465000059, -8.397648071999924 11.385895285000103, -8.375323851999951 11.367653504000145, -8.393358927999884 11.363002625000078, -8.414132852999927 11.353468323000087, -8.426070108999852 11.340471700000066, -8.41795690899994 11.325511373000097, -8.39956009999986 11.322152405000082, -8.38126664299989 11.325097961000083, -8.3710347089999 11.320550436000048, -8.376719115999947 11.294660543000106, -8.38839798999993 11.281405538000115, -8.4066914469999 11.275281881000055, -8.426896931999892 11.27419667600006, -8.44457027199988 11.276263733000107, -8.4832242439999 11.284919536000075, -8.492577677999918 11.27778818800003, -8.503119669999904 11.255412293000049, -8.50234452399988 11.251639913000119, -8.499243936999932 11.246704814000083, -8.497021850999943 11.240710348000064, -8.498778848999905 11.233501486000094, -8.502137817999937 11.231176046000073, -8.512059692999912 11.23029754700012, -8.515625365999881 11.229057312000094, -8.52410030099989 11.223760478000145, -8.544615844999953 11.214975485000068, -8.551798868999924 11.210738017000025, -8.563942830999906 11.193426412000107, -8.569678914999884 11.157795512000035, -8.578153848999875 11.142060039000057, -8.586990518999897 11.137460836000102, -8.607221842999877 11.13715077700003, -8.6154125569999 11.134050191000071, -8.624610961999906 11.122035421000106, -8.624921019999903 11.112191061000132, -8.622388875999889 11.102605083000128, -8.623112344999896 11.091262106000059, -8.63368017599987 11.071289165000024, -8.679801391999888 11.025193786000102, -8.688689737999908 11.01255889900007, -8.696906290999948 10.996048279000078, -8.700730346999876 10.978452454000092, -8.696234496999864 10.962742819000127, -8.680938272999896 10.952304179000095, -8.665228637999945 10.956955058000062, -8.64796870899994 10.96581756600007, -8.62822831199989 10.967832947000105, -8.611381794999915 10.965145773000089, -8.593346720999904 10.964706522000071, -8.5756217039999 10.966747742000024, -8.559550333999937 10.97160532600013, -8.541566935999924 10.983775127000044, -8.531645059999931 10.999045512000095, -8.505755167999922 11.052789001000093, -8.489063679999958 11.05273732500008, -8.469013224999912 11.046381124000064, -8.448704385999918 11.04710459500005, -8.438162393999932 11.049740092000064, -8.429067341999883 11.04829315200007, -8.410463825999926 11.04191111300014, -8.400645303999852 11.04333221500012, -8.381318318999888 11.053925883000105, -8.369329385999947 11.054235941000087, -8.347315225999864 11.043125509000077, -8.326592976999848 11.023850199000037, -8.311141723999924 11.000802511000074, -8.305095581999893 10.978245748000049, -8.30364864099991 10.874892883000058, -8.305870727999888 10.856289368000105, -8.311813516999905 10.847245993000058, -8.318893188999937 10.839520366000059, -8.339822142999935 10.783864848000036, -8.341940877999889 10.76324595200012, -8.315585896999949 10.75309153200007, -8.311141723999924 10.73528900200003, -8.295535440999885 10.53744578000007, -8.284425007999886 10.511245830000064, -8.25388423699988 10.469827169000013, -8.22887284399988 10.423421733000097, -8.21057938699991 10.41406829900005, -8.189030313999865 10.413629048000132, -8.165982624999884 10.418099060000145, -8.150583047999902 10.424532776000092, -8.142056436999894 10.429648743000087, -8.135441853999907 10.427323303000051, -8.125416625999975 10.411665345000102, -8.122109333999958 10.39949554400009, -8.120920776999952 10.373037212000014, -8.115494750999886 10.360247295000121, -8.096116088999906 10.346733907000058, -8.070071166999924 10.341876322000132, -8.015707560999942 10.339344178000147, -7.996277220999986 10.328259583000133, -7.981549438999878 10.307020569000073, -7.971214151999959 10.281802470000116, -7.964599568999887 10.258935649000094, -7.962584187999937 10.233252462000024, -7.968010213999889 10.209998068000104, -7.989662638999988 10.161990662000065))')
#     t, bin3 = x.compress(faulty_vertices)
#     # binary = bitarray()
#     # binary.frombytes(bin3)
#     # util.pprint(binary)
#     # print(bin3.hex(sep=' '))
#     add_ifx = 17
#     # for i in range(11):
#     #     t, bin3 = x.add_vertex((bin3, 280, (13.1945576, 55.7203302)))
#     #     print(bin3.hex(sep=' '))
#     #     print("---")
    
#     print(x.vertices(bin3)) 
#     t, decomp = x.decompress(bin3)
#     print(decomp)
#     print(bin3)

#     # print(x.bytes_to_float32(b'0xa70x300x530x41'))
#     # print(bin3.hex(sep=' '))
#     # print("-DELTASIZE-_TY_POINTSINCHK_XFIRST-----------------_YFIRST-----------------")
#     # print("-DELTASIZE-_TY_RINGS------_POINTSINCHK_XFIRST-----------------_YFIRST-----------------_--XD1---------_--YD1---------_--XD2---------_--YD2---------_--XD3---------_--YD3---------_--XD4---------_--YD4---------_--XD5---------_--YD5---------_--XD6---------_--YD6---------_--XD7---------_--YD7---------_--XD8---------_--YD8---------_--XD9---------_--YD9---------")
#     # t = shapely.to_ragged_array([geom1])
#     # print(t[2])

#     # print(x.bytes_to_float(bin3, size=16))
#     # print(x.calculate_delta_size(geom1))
#     # t, bin = x.decompress(bin)

#     # print(bin(0b01100111 | 0b0000000000000000))
#     # t = '{0:0{1}b}'.format(0b101, 10)
#     # print(t)
#     # print(int(t, 2))
#     # print(x.get_poly_ring_count(geom5))

#     # print(len(bytes(shapely.to_wkt(geom4), 'utf-8')), shapely.to_wkt(geom4))
#     # print(shapely.to_wkt(faulty_add, rounding_precision=-1))
#     # print(shapely.to_wkt(decomp, rounding_precision=-1))
#    #print("LEN COMP", len(bin3))

#     # fpd = Fpd()
#     # _, org = fpd.compress(faulty_add)
#     # _, org = fpd.add_vertex((org, add_idx, (13.196935700000001, 55.6907853)))

#     # _, de = fpd.decompress(org)   
#     # print(decomp == de)
#     # print(shapely.to_wkt(de, rounding_precision=-1))

#     # print(x.bounding_box(bin3))

#     return

#     DATASET_PATH = "data/lund_building_highway.json"
#     # DATASET_PATH = "data/world.json"
#     NBR_ITER = 16000

#     SEED = 123  # If we want to enforce the same ordering and indexes for multiple runs, else None
#     random.seed(SEED)  # Init random

#     # Extract the nested feature attribute of the geo_json file containing the geometries
#     with open(DATASET_PATH, 'r') as f:
#         data = json.loads(f.read())
#     file_df: pd.DataFrame = pd.json_normalize(data, record_path=['features'])
#     # Create a dataframe suitable for the WKT format for easy convertion to shapely objects
#     df = pd.DataFrame(
#         {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

#     max_idx = len(df) - 1
#     unary_idxs = [random.randint(0, max_idx) for i in range(NBR_ITER)]  # Generate list of indexes to query on
#     random.seed(SEED)  # Reset random

#     # Compress files, benchmark unaries
#     for idx in tqdm.tqdm(unary_idxs):  # List of single idxs
#         t, comp = x.compress(shape(df.iloc[idx]))
#         t, decomp = x.decompress(comp)
#         if (shape(df.iloc[idx]) != decomp):
#             print("FAIL")


if __name__ == "__main__":
    main()
