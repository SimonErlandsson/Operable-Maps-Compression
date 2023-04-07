from collections import deque
import shapely
import math
from bitarray import bitarray, util, bits2bytes
import time
from shapely import GeometryType as GT
from algos.fpd_extended_lib.intersection_chunk_bbox_wrapper import *
from algos.fpd_extended_lib.low_level import *

def get_zz_encoded_delta(prev_coord, curr_coord):
    return zz_encode(double_as_long(curr_coord) - double_as_long(prev_coord))

def deltas_fit_in_bits(d_x, d_y, max_bits):
    return (d_x == 0 or math.log2(d_x) < max_bits) and (d_y == 0 or math.log2(d_y) < max_bits)

# Returns the number of coords in each ring for a polygon. If multipolygon, index on polygon first
def point_count(geometry):
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

def append_header(bits, geometry, d_size):
    # Meta data
    bits.frombytes(uchar_to_bytes(d_size))
    bits.frombytes(uchar_to_bytes(int(shapely.get_type_id(geometry))))  # 1 byte is enough for storing type
    # Bounding Box
    bounds = shapely.bounds(geometry)
    bits.frombytes(double_to_bytes(bounds[0]))
    bits.frombytes(double_to_bytes(bounds[1]))
    bits.frombytes(double_to_bytes(bounds[2]))
    bits.frombytes(double_to_bytes(bounds[3]))
    intersection_reserve_header(bits)


def append_delta_pair(bits, d_x_zig, d_y_zig, d_size):
        x_bytes = uint_to_ba(d_x_zig, d_size)
        y_bytes = uint_to_ba(d_y_zig, d_size)
        bits.extend(x_bytes)
        bits.extend(y_bytes)

def fp_delta_encoding(geometry, d_size):
    # List of resulting bytes.
    bits = bitarray(endian='big')
    # Init with 'd_size', 'geom_type'
    append_header(bits, geometry, d_size)

    # Type specific variables
    geo_type = shapely.get_type_id(geometry)
    is_linestring = geo_type == GT.LINESTRING
    is_multipolygon = geo_type == GT.MULTIPOLYGON
    is_polygon = geo_type == GT.POLYGON

    # Fetches number of points in each ring, nestled for multipoly
    poly_buffer = point_count(geometry)
    ring_buffer = poly_buffer if is_polygon else []  # Not nestled for poly, else overwritten below

    prev_x, prev_y = 0, 0  # Absolute value of previous coord
    chk_deltas = 0  # Cnt of 'deltas in chunk'
    chk_deltas_idx = 0  # Pointer to 'deltas of chunk'
    num_chks_ring = 0  # Cnt of 'number of chunks for current ring'
    num_chks_ring_idx = 0  # Pointer to latest 'number of chunks for ring'
    rem_points_ring = 0  # Cnt of 'points left to process in current ring'

    intersection_chunk_bboxes = []
    def intersection_add_point(x, y, previous_chunk=False):
        i = -2 if previous_chunk else -1
        x_l, y_b, x_r, y_t = intersection_chunk_bboxes[i]
        intersection_chunk_bboxes[i] = [min(x, x_l), min(y, y_b), max(x, x_r), max(y, y_t)]

    # Loop all coordinates
    for x, y in shapely.get_coordinates(geometry):
        if not is_linestring and rem_points_ring == 1:  # Is the whole ring processed? We skip last coordinate
            # Store number of chunks used for the ring
            rem_points_ring = 0
            bits[num_chks_ring_idx:num_chks_ring_idx + RING_CHK_CNT_SIZE] = uint_to_ba(num_chks_ring, RING_CHK_CNT_SIZE)
            num_chks_ring = 0
            intersection_add_point(*intersection_first_coord_ring)
            continue  # Skip last coordinate
        d_x_zig = get_zz_encoded_delta(prev_x, x)  # Calculated delta based on previous iteration
        d_y_zig = get_zz_encoded_delta(prev_y, y)
        prev_x, prev_y = (x, y)

        # ---- CREATE NEW CHUNK? If 'first chunk', 'delta doesn't fit', 'new ring', or 'reached max deltas'
        if chk_deltas_idx == 0 or not deltas_fit_in_bits(d_x_zig, d_y_zig, d_size) or rem_points_ring == 0 or chk_deltas == MAX_NUM_DELTAS:
            # If not 'first chunk' -> save previous chunk's size
            if chk_deltas_idx != 0:
                bits[chk_deltas_idx:chk_deltas_idx + D_CNT_SIZE] = uint_to_ba(chk_deltas, D_CNT_SIZE)

            ###### ---- INITIALIZE NEW CHUNK ----- ######
            chk_deltas = 0
            intersection_chunk_bboxes.append([x, y, x, y])

            ### __ RING/MULTI-POLYGON META-DATA __ ###
            if not is_linestring:
                # Ran out of points -> fetch number of points in next ring
                if rem_points_ring == 0:
                    # Check if we ran out of rings -> fetch rings of NEXT POLYGON
                    if is_multipolygon and len(ring_buffer) == 0:
                        ring_buffer = poly_buffer.popleft()
                        bits.extend(uint_to_ba(len(ring_buffer), POLY_RING_CNT_SIZE))  # Append 'nbr of rings in poly'
                    # Set 'remaining points' to cnt in new ring
                    rem_points_ring = ring_buffer.popleft()
                    num_chks_ring_idx = len(bits)
                    bits.extend(uint_to_ba(0, RING_CHK_CNT_SIZE))  # Reserve space for number of chunks for current ring
                    num_chks_ring = 1
                    intersection_first_coord_ring = (x, y)
                else:
                    num_chks_ring += 1
                    intersection_add_point(x, y, previous_chunk=True)
            ### __ ------------ END ------------- __ ###

            # Preparing chunk size (number of deltas)
            chk_deltas_idx = len(bits)
            bits.extend(uint_to_ba(0, D_CNT_SIZE))  # Reserve space for size

            # Add full coordinates
            bits.frombytes(double_to_bytes(x))
            bits.frombytes(double_to_bytes(y))
        else:
            # Delta fits, append it
            append_delta_pair(bits, d_x_zig, d_y_zig, d_size)
            chk_deltas += 1
            intersection_add_point(x, y)

        # Coord has been processed, remove it
        rem_points_ring -= 1

    # All points processed. Update size of final chunk
    bits[chk_deltas_idx:chk_deltas_idx + D_CNT_SIZE] = uint_to_ba(chk_deltas, D_CNT_SIZE)
    bits = intersection_append_header(bits, intersection_chunk_bboxes)

    # util.pprint(bits)
    # print([int.from_bytes(i, 'big') for i in bytes], '\n')
    return bits.tobytes()

def calculate_delta_size(geometry, return_deltas=False):
    deltas = [[], []]
    RESET_POINT_SIZE = FLOAT_SIZE * 2 + D_CNT_SIZE
    coords = shapely.get_coordinates(geometry)
    prev = [0, 0]
    bit_cnts = {}
    for coord in coords:
        bit_cnt = 0
        for i in range(2):
            d = get_zz_encoded_delta(prev[i], coord[i])
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
    optimal_size, _, _ = calculate_delta_size(geometry)
    bin = fp_delta_encoding(geometry, optimal_size)
    t = time.perf_counter()
    return t - s, bin