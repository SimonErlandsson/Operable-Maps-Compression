
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import math
import numpy as np
from shapely import GeometryType as GT
from bitarray import bitarray, util, bits2bytes
from algos.fpd_extended_lib.cfg import *
from algos.fpd_extended_lib.low_level import *

chunk_bboxes = []
chunk_bounds_offset = -1

INTERSECTION_CHK_CNT_SIZE = 32

def get_chunk_bboxes_len():
    return len(chunk_bboxes)

def intersection_reserve_header(bits):
    if DISABLE_OPTIMIZED_INTERSECTION:
        return
    global chunk_bounds_offset
    global chunk_bboxes
    chunk_bounds_offset = len(bits)
    chunk_bboxes = []

def intersection_new_chunk():
    if DISABLE_OPTIMIZED_INTERSECTION:
        return
    global chunk_bboxes
    chunk_bboxes.append([999.0, 999.0, -999.0, -999.0])

def intersection_add_point(x, y, previous_chunk=False):
    if DISABLE_OPTIMIZED_INTERSECTION:
        return
    i = -2 if previous_chunk else -1
    x_l, y_b, x_r, y_t = chunk_bboxes[i]
    chunk_bboxes[i] = [min(x, x_l), min(y, y_b), max(x, x_r), max(y, y_t)]

def intersection_append_header(bits):
    if DISABLE_OPTIMIZED_INTERSECTION:
        return bits
        
    left = bits[0:chunk_bounds_offset]
    left.extend(uint_to_ba(len(chunk_bboxes), INTERSECTION_CHK_CNT_SIZE)) # Store nbr of chunks
    for bbox in chunk_bboxes:
        for i in range(4):
            left.frombytes(double_to_bytes(bbox[i]))
    right = bits[chunk_bounds_offset:]
    left.extend(right)
    return left

def intersection_skip_header(bin):
    if DISABLE_OPTIMIZED_INTERSECTION:
        return
    chk_cnt = struct.unpack_from('!I', bin, offset=cfg.offset//8)[0]
    cfg.offset += INTERSECTION_CHK_CNT_SIZE + 4 * FLOAT_SIZE * chk_cnt

def get_chunk_bounds(bin_in):
    init_offset = 3 * 8 if cfg.USE_ENTROPY else 2 * 8
    cfg.offset = init_offset + (4 * FLOAT_SIZE if not cfg.DISABLE_OPTIMIZED_BOUNDING_BOX else 0) # Skip normal header
    chk_cnt = struct.unpack_from('!I', bin_in, offset=cfg.offset//8)[0]
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)
    cfg.offset += INTERSECTION_CHK_CNT_SIZE
    bounds = []
    for _ in range(chk_cnt):
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
