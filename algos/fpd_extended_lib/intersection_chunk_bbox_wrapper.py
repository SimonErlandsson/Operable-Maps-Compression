
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import math
import numpy as np
from shapely import GeometryType as GT
from bitarray import bitarray, util, bits2bytes
import algos.fpd_extended_lib.cfg as cfg
from algos.fpd_extended_lib.low_level import *

chunk_bboxes = []
chunk_bounds_offset = -1

INTERSECTION_CHK_CNT_SIZE = 32

def get_chunk_bboxes_len():
    return len(chunk_bboxes)

def intersection_reserve_header(bits):
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        return
    global chunk_bounds_offset
    global chunk_bboxes
    chunk_bounds_offset = len(bits)
    chunk_bboxes = []

def intersection_new_chunk():
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        return
    global chunk_bboxes
    chunk_bboxes.append([999.0, 999.0, -999.0, -999.0])

def intersection_add_point(x, y, previous_chunk=False):
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        return
    i = -2 if previous_chunk else -1
    x_l, y_b, x_r, y_t = chunk_bboxes[i]
    chunk_bboxes[i] = [min(x, x_l), min(y, y_b), max(x, x_r), max(y, y_t)]

def intersection_append_header(bits):
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        return bits
        
    left = bits[0:chunk_bounds_offset]
    if cfg.DELTA_ENCODE_CHUNK_BBOXES:
        from algos.fpd_extended_lib.compress import append_delta_pair, deltas_fit_in_bits, get_zz_encoded_delta, calculate_delta_size
        coords = [[bbox[2 * i], bbox[2 * i + 1]] for i in range(2) for bbox in chunk_bboxes]
        d_size = calculate_delta_size(coords=coords)[0]
        chk_dt_cnt = 0
        prev_x, prev_y = 0, 0
        chk_hdr_offset = 0
        left.extend(uint_to_ba(0, INTERSECTION_CHK_CNT_SIZE)) # Reserve space for 'total len'
        total_len_offset = len(left)
        left.frombytes(uchar_to_bytes(d_size)) # Append delta size
        for x, y in coords:
            d_x_zig = get_zz_encoded_delta(prev_x, x)  # Calculated delta based on previous iteration
            d_y_zig = get_zz_encoded_delta(prev_y, y)
            prev_x, prev_y = x, y

            # ---- CREATE NEW CHUNK? If 'first chunk', 'delta doesn't fit', 'new ring', or 'reached max deltas'
            if chk_hdr_offset == 0 or not deltas_fit_in_bits(d_x_zig, d_y_zig, d_size):
                # If not 'first chunk' -> save previous chunk's size
                if chk_hdr_offset != 0:
                    left[chk_hdr_offset:chk_hdr_offset + cfg.D_CNT_SIZE] = uint_to_ba(chk_dt_cnt, cfg.D_CNT_SIZE)
                chk_dt_cnt = 0

                # Preparing chunk size (number of deltas)
                chk_hdr_offset = len(left)
                left.extend(uint_to_ba(0, cfg.D_CNT_SIZE)) # Reserve space for 'chk_dt_cnt'
                # Add full coordinates
                left.frombytes(double_to_bytes(x))
                left.frombytes(double_to_bytes(y))
            else:
                # Delta fits, append it
                append_delta_pair(left, d_x_zig, d_y_zig, d_size)
                chk_dt_cnt += 1

        # All points processed. Update size of final chunk
        left[chk_hdr_offset:chk_hdr_offset + cfg.D_CNT_SIZE] = uint_to_ba(chk_dt_cnt, cfg.D_CNT_SIZE)
        # Store final length
        left[total_len_offset-INTERSECTION_CHK_CNT_SIZE:total_len_offset] = uint_to_ba(len(left) - total_len_offset, INTERSECTION_CHK_CNT_SIZE)

        right = bits[chunk_bounds_offset:]
        left.extend(right)
    else:
        left.extend(uint_to_ba(len(chunk_bboxes), INTERSECTION_CHK_CNT_SIZE)) # Store nbr of chunks
        for bbox in chunk_bboxes:
            for i in range(4):
                left.frombytes(double_to_bytes(bbox[i]))
        right = bits[chunk_bounds_offset:]
        left.extend(right)
    return left

def intersection_skip_header(bin):
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        return
    if cfg.DELTA_ENCODE_CHUNK_BBOXES:
        cfg.offset += struct.unpack_from('!I', bin, offset=cfg.offset//8)[0]
    else:
        chk_cnt = struct.unpack_from('!I', bin, offset=cfg.offset//8)[0]
        cfg.offset += INTERSECTION_CHK_CNT_SIZE + 4 * cfg.FLOAT_SIZE * chk_cnt

def get_chunk_bounds(bin_in):
    pre_header_len = 3 * 8 if cfg.USE_ENTROPY else 2 * 8
    pre_header_len += (4 * cfg.FLOAT_SIZE if not cfg.DISABLE_OPTIMIZED_BOUNDING_BOX else 0) # Skip normal header
    cfg.offset = pre_header_len 
    total_len = struct.unpack_from('!I', bin_in, offset=cfg.offset//8)[0]
    cfg.offset += INTERSECTION_CHK_CNT_SIZE
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)
    start_offset = cfg.offset
    bounds = []

    if cfg.DELTA_ENCODE_CHUNK_BBOXES:
        def add_coord(x, y, prev):
            if prev == None:
                return (x,  y)
            else:
                bounds.append([prev[0], prev[1], x, y])
                return None
            
        delta_size = struct.unpack_from('!B', bin_in, offset=cfg.offset//8)[0]
        cfg.offset += 8

        prev = None
        while cfg.offset < start_offset + total_len:
            chk_size = bytes_to_uint(bin, cfg.D_CNT_SIZE)
            # Extract reset point
            x = bytes_to_double(bin)
            y = bytes_to_double(bin)
            prev = add_coord(x, y, prev)
            print(x,  y)
            

            # Loop through deltas in chunk
            for _ in range(chk_size):
                x = bytes_to_decoded_coord(bin, x, delta_size)
                y = bytes_to_decoded_coord(bin, y, delta_size)
                print(x,y)
                prev = add_coord(x, y, prev)
    else:
        for _ in range(total_len):
            bounds.append([bytes_to_double(bin), bytes_to_double(bin), bytes_to_double(bin), bytes_to_double(bin)])
    return bounds

def is_intersecting(self, args):
    from intersection.chunk_bbox_intersection import is_intersecting as intersects # Prevent circular import
    s = time.perf_counter()
    res = intersects(args)
    t = time.perf_counter()
    return t - s, res


def intersection(self, args):
    from intersection.chunk_bbox_intersection import intersection as intersect # Prevent circular import
    s = time.perf_counter()
    res = intersect(args)
    t = time.perf_counter()

    return t - s, res
