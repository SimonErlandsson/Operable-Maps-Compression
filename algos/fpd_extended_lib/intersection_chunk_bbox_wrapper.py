
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

def intersection_reserve_header(bits):
    global chunk_bounds_offset, chunk_bboxes
    chunk_bounds_offset = len(bits)
    chunk_bboxes = []

def intersection_new_chunk():
    global chunk_bboxes
    chunk_bboxes.append([999.0, 999.0, -999.0, -999.0])

def intersection_add_point(x, y, previous_chunk=False):
    i = -2 if previous_chunk else -1
    x_l, y_b, x_r, y_t = chunk_bboxes[i]
    chunk_bboxes[i] = [min(x, x_l), min(y, y_b), max(x, x_r), max(y, y_t)]

def intersection_append_header(bits):
    left = bits[0:chunk_bounds_offset]
    left.extend(uint_to_ba(len(chunk_bboxes), 32))
    for bbox in chunk_bboxes:
        for i in range(4):
            left.frombytes(double_to_bytes(bbox[i]))
    right = bits[chunk_bounds_offset:]
    left.extend(right)
    return left

def intersection_skip_header(bin):
    chk_cnt = struct.unpack_from('!I', bin, offset=cfg.offset//8)[0]
    cfg.offset += 32 + 4 * FLOAT_SIZE * chk_cnt

def get_chunk_bounds(bin_in):
    cfg.offset = 3 * 8 + 4 * FLOAT_SIZE # Skip normal header
    chk_cnt = struct.unpack_from('!I', bin_in, offset=cfg.offset//8)[0]
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)
    cfg.offset += 32
    bounds = []
    for _ in range(chk_cnt):
        bounds.append([bytes_to_double(bin), bytes_to_double(bin), bytes_to_double(bin), bytes_to_double(bin)])
    return bounds

def is_intersecting(self, args):
    from intersection.chunk_bbox_intersection import is_intersecting as intersects # Prevent circular import
    l_bin, r_bin = args
    s = time.perf_counter()

    res = intersects(args)
    # _, l_geo = self.decompress(l_bin)
    # _, r_geo = self.decompress(r_bin)
    # res = shapely.intersects(l_geo, r_geo)

    t = time.perf_counter()
    return t - s, res


def intersection(self, args):
    from intersection.chunk_bbox_intersection import intersection as intersect # Prevent circular import
    l_bin, r_bin = args
    s = time.perf_counter()

    #res = intersect(args)
    _, l_geo = self.decompress(l_bin)
    _, r_geo = self.decompress(r_bin)
    res = shapely.intersection(l_geo, r_geo)

    t = time.perf_counter()

    return t - s, res
