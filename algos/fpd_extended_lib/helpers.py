from bitarray import bitarray, util
from shapely import GeometryType as GT
from algos.fpd_extended_lib.cfg import *
from algos.fpd_extended_lib.low_level import *
from algos.fpd_extended_lib.decompress import *

# Only used for non-timed operations, i.e. debugging/testing implementations
# NOTE: Also returns the first ring coordinate when reaching end of ring!
def get_chunks(bin_in):
    chunks = []
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

    is_last_ring_chunk = []
    bin_len = len(bin)
    while (cfg.offset + EOF_THRESHOLD <= bin_len):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            chunks_in_ring = chunks_in_ring_left

        # Go through chunk (inlined sequence decode)
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)
        # Extract reset point
        x = bytes_to_double(bin)
        y = bytes_to_double(bin)
        if chunks_in_ring_left == chunks_in_ring:
            x_ring, y_ring = (x, y)

        chunk = [[x, y]]
        # Loop through deltas in chunk
        for _ in range(deltas_in_chunk):
            x = bytes_to_decoded_coord(bin, x, delta_size)
            y = bytes_to_decoded_coord(bin, y, delta_size)
            chunk.append([x, y])
        chunks.append(chunk)
        is_last_ring_chunk.append(False)
        chunks_in_ring_left -= 1
        if chunks_in_ring_left == 0:
            chunks[-1].append([x_ring, y_ring])
            rings_left -= 1
            is_last_ring_chunk[-1] = True
    return chunks, is_last_ring_chunk


def access_vertex(bin_in, access_idx, cache=[], getBoundsData=False):
    if not getBoundsData and access_idx in cache:
        return cache[access_idx], cache, None
    old_offset = cfg.offset
    cfg.offset = 0
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)

    if 'header' in cache:
        delta_size, type = cache['header']
        cfg.offset = cache['offset']
    else:
        delta_size, type = decode_header(bin)
    # Type specific variables
    is_linestring = type == GT.LINESTRING
    is_multipolygon = type == GT.MULTIPOLYGON
    is_polygon = type == GT.POLYGON

    idx_found = False
    ring_beg_idx, ring_end_idx, poly_beg_idx, poly_end_idx = None, None, None, None

    p_idx = 0
    chunks_in_ring_left = 0  # Used for iteration
    rings_left = 0
    while (p_idx <= access_idx or getBoundsData):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
            if getBoundsData:
                if not idx_found:
                    poly_beg_idx = p_idx

        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            if not idx_found:
                ring_beg_idx = p_idx

        deltas_in_chunk_offset = cfg.offset
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        # Found chunk containing vertex?
        if not idx_found and (p_idx <= access_idx and access_idx <= p_idx + deltas_in_chunk):
            idx_found = True
            (x, y), cache = access_vertex_chk(bin, deltas_in_chunk_offset, access_idx - p_idx, delta_size, cache, offset_idx=p_idx)
            if not getBoundsData:
                break

        # Jump to next chunk
        p_idx += 1 + deltas_in_chunk
        cfg.offset += FLOAT_SIZE * 2 + delta_size * 2 * deltas_in_chunk
        chunks_in_ring_left -= 1
        if (chunks_in_ring_left == 0):
            if idx_found:
                ring_end_idx = p_idx - 1
            rings_left -= 1
        if (rings_left == 0):
            if idx_found and is_multipolygon:
                poly_end_idx = p_idx - 1

        if getBoundsData and ((is_multipolygon and poly_end_idx != None) or (is_polygon and ring_end_idx != None) or is_linestring):
            break
    cfg.offset = old_offset
    return (x, y), cache, None if not getBoundsData else (ring_beg_idx, ring_end_idx, poly_beg_idx, poly_end_idx)

    # Supply the offset to D_CNT, and idx is the index within the chunk
    # @profile


def access_vertex_chk(bin, chk_offset, idx, delta_size, cache=None, offset_idx=0):
    if cache != None and idx + offset_idx in cache:
        return cache[idx + offset_idx], cache
    old_offset = cfg.offset
    cfg.offset = chk_offset + D_CNT_SIZE
    # Extract reset point
    x, y = (bytes_to_double(bin), bytes_to_double(bin))
    # Loop through deltas in chunk
    for idx in range(idx):
        if cache != None and idx + offset_idx in cache:
            cfg.offset += delta_size * 2
            (x, y) = cache[idx + offset_idx]
        else:
            x = bytes_to_decoded_coord(bin, x, delta_size)
            y = bytes_to_decoded_coord(bin, y, delta_size)
            if cache != None:
                cache[idx + offset_idx] = (x, y)
    cfg.offset = old_offset
    return (x, y), cache

def access_vertex(bin_in, access_idx, cache=[], getBoundsData=False):
    if not getBoundsData and access_idx in cache:
        return cache[access_idx], cache, None
    old_offset = cfg.offset
    cfg.offset = 0
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)

    if 'header' in cache:
        delta_size, type = cache['header']
        cfg.offset = cache['offset']
    else:
        delta_size, type = decode_header(bin)
    # Type specific variables
    is_linestring = type == GT.LINESTRING
    is_multipolygon = type == GT.MULTIPOLYGON
    is_polygon = type == GT.POLYGON

    idx_found = False
    ring_beg_idx, ring_end_idx, poly_beg_idx, poly_end_idx = None, None, None, None

    p_idx = 0
    chunks_in_ring_left = 0  # Used for iteration
    rings_left = 0
    while (p_idx <= access_idx or getBoundsData):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
            if getBoundsData:
                if not idx_found:
                    poly_beg_idx = p_idx

        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            if not idx_found:
                ring_beg_idx = p_idx

        deltas_in_chunk_offset = cfg.offset
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        # Found chunk containing vertex?
        if not idx_found and (p_idx <= access_idx and access_idx <= p_idx + deltas_in_chunk):
            idx_found = True
            (x, y), cache = access_vertex_chk(bin, deltas_in_chunk_offset, access_idx - p_idx, delta_size, cache, offset_idx=p_idx)
            if not getBoundsData:
                break

        # Jump to next chunk
        p_idx += 1 + deltas_in_chunk
        cfg.offset += FLOAT_SIZE * 2 + delta_size * 2 * deltas_in_chunk
        chunks_in_ring_left -= 1
        if (chunks_in_ring_left == 0):
            if idx_found:
                ring_end_idx = p_idx - 1
            rings_left -= 1
        if (rings_left == 0):
            if idx_found and is_multipolygon:
                poly_end_idx = p_idx - 1

        if getBoundsData and ((is_multipolygon and poly_end_idx != None) or (is_polygon and ring_end_idx != None) or is_linestring):
            break
    cfg.offset = old_offset
    return (x, y), cache, None if not getBoundsData else (ring_beg_idx, ring_end_idx, poly_beg_idx, poly_end_idx)

    # Supply the offset to D_CNT, and idx is the index within the chunk
def access_vertex_chk(bin, chk_offset, idx, delta_size, cache=None, offset_idx=0):
    if cache != None and idx + offset_idx in cache:
        return cache[idx + offset_idx], cache
    old_offset = cfg.offset
    cfg.offset = chk_offset + D_CNT_SIZE
    # Extract reset point
    x, y = (bytes_to_double(bin), bytes_to_double(bin))
    # Loop through deltas in chunk
    for idx in range(idx):
        if cache != None and idx + offset_idx in cache:
            cfg.offset += delta_size * 2
            (x, y) = cache[idx + offset_idx]
        else:
            x = bytes_to_decoded_coord(bin, x, delta_size)
            y = bytes_to_decoded_coord(bin, y, delta_size)
            if cache != None:
                cache[idx + offset_idx] = (x, y)
    cfg.offset = old_offset
    return (x, y), cache