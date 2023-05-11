from collections import deque
import time
import struct
import math
import numpy as np
from shapely import GeometryType as GT
import bisect
from bitarray import bitarray, util
from algos.fpd_extended_lib.cfg import *
from algos.fpd_extended_lib.decompress import *
from algos.fpd_extended_lib.helpers import decompress_chunk


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

def vertices(self, bin_in):
    s = time.perf_counter()

    coords = []
    cfg.offset = 0
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)

    delta_size, type = decode_header(bin)
    # Type specific variables
    is_linestring = type == GT.LINESTRING
    is_multipolygon = type == GT.MULTIPOLYGON

    chunks_in_ring_left = 0  # Used for iteration
    chunks_in_ring = 0
    rings_left = 0
    cfg.binary_length = len(bin)
    while (cfg.offset + EOF_THRESHOLD <= cfg.binary_length):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            chunks_in_ring = chunks_in_ring_left

        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        if cfg.COMPRESS_CHUNK or cfg.USE_ENTROPY:
            delta_bytes_size = deltas_in_chunk * delta_size * 2 - bytes_to_int(bin, D_BITSIZE_SIZE)

        # Extract reset point
        x = bytes_to_double(bin)
        y = bytes_to_double(bin)

        if COMPRESS_CHUNK:
            chk_deltas_offset = cfg.offset
            bin, coord_bit_len = decompress_chunk(bin, chk_deltas_offset, delta_bytes_size) 
   
        if chunks_in_ring_left == chunks_in_ring:
            x_ring, y_ring = (x, y)
        coords.append([x, y])
        # Loop through deltas in chunk
        for _ in range(deltas_in_chunk):
            x = bytes_to_decoded_coord(bin, x, delta_size)
            y = bytes_to_decoded_coord(bin, y, delta_size)
            coords.append([x, y])
        if COMPRESS_CHUNK:
            cfg.offset = chk_deltas_offset + coord_bit_len

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
    if not DISABLE_OPTIMIZED_BOUNDING_BOX:
        res = bitarray()
        res.frombytes(bin)
        bounds = [bin_to_double(res[3 * 8 + FLOAT_SIZE * i: 3 * 8 + FLOAT_SIZE * (i + 1)]) for i in range(4)]
    else:
        _, geometry = self.decompress(bin)
        bounds = shapely.bounds(geometry)
    t = time.perf_counter()
    return t - s, bounds
