import math
import zlib
import gzip
from bitarray import bitarray, util
from shapely import GeometryType as GT
from algos.fpd_extended_lib.cfg import *
from algos.fpd_extended_lib.low_level import *
from algos.fpd_extended_lib.decompress import *


def get_chunks(bin_in, include_next_start=True):
    """
    # DEBUG ONLY!
    Only used for non-timed operations, i.e. debugging/testing implementations
    # NOTE: Can also return the first ring coordinate when reaching end of ring!
    # DO NOT USE ENTROPY OR CHUNK COMPRESSION WHEN CALLING METHOD.
    """
    assert(not(cfg.USE_ENTROPY or cfg.COMPRESS_CHUNK))
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

        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        # Extract reset point
        x = bytes_to_double(bin)
        y = bytes_to_double(bin)
        if chunks_in_ring_left == chunks_in_ring:
            x_ring, y_ring = (x, y)
        elif include_next_start:
            # Append first coord next chunk
            chunks[-1].append([x, y])

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
            if include_next_start:
                chunks[-1].append([x_ring, y_ring])
            rings_left -= 1
            is_last_ring_chunk[-1] = True
    return chunks, is_last_ring_chunk


def access_vertex(bin_in, idx, cache=[]):
    return random_access(bin_in, idx, cache, get_chunk=False)

def access_chunk(bin_in, idx, offset_cache=None, cache=[]):
    """
    Get a chunk based on chunk index. Note that
    chunks are connected. I.e. each chunk contains the
    next point in sequence.
    """
    return random_access(bin_in, idx, cache, offset_cache=offset_cache, get_chunk=True)


# Supply cache if used repeatedly for same shape
#@profile
def random_access(bin_in, idx, cache, offset_cache=None, get_chunk=False):
    old_offset = cfg.offset
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)
    cfg.binary_length = len(bin)

    #Use offset_cache if possible
    if offset_cache != None and  len(offset_cache)  != 0:
        #Go to last index before the requested one, else go directly to it 
        delta_size, type, cur_idx, chunks_in_ring_left, rings_left, cfg.offset, ring_start_offset = offset_cache[idx if len(offset_cache) > idx else len(offset_cache)  - 1] 

    else:
            cfg.offset = 0
            ring_start_offset = None
            delta_size, type = decode_header(bin)
            cur_idx = 0
            chunks_in_ring_left = 0  # Used for iteration
            rings_left = 0

    # Type specific variables
    is_linestring = type == GT.LINESTRING
    is_multipolygon = type == GT.MULTIPOLYGON
    
   
    while (cur_idx <= idx):
        if offset_cache != None and cur_idx not in offset_cache:
            offset_cache[cur_idx] = delta_size, type, cur_idx, chunks_in_ring_left, rings_left, cfg.offset, ring_start_offset

        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            ring_start_offset = cfg.offset
        deltas_in_chunk_offset = cfg.offset
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)


        delta_bytes_size = deltas_in_chunk * delta_size * 2
        if cfg.COMPRESS_CHUNK or cfg.USE_ENTROPY:
            delta_bytes_size -= bytes_to_int(bin, D_BITSIZE_SIZE)

        if get_chunk and cur_idx == idx:
            # Looking for whole chunk, found it
            chk, cache = access_vertex_chk(bin, deltas_in_chunk_offset, delta_size, cache=cache, list_vertices=True)

            # Chunks should contain next point for overlapping bboxes
            if not is_linestring and chunks_in_ring_left == 1:
                # Was last chunk in ring, append start
                next_chk_offset = ring_start_offset
            else:
                # Append next chunk start
                next_chk_offset = cfg.offset + FLOAT_SIZE * 2 + delta_bytes_size
            # Avoid reading if EOF
            if next_chk_offset + cfg.EOF_THRESHOLD <= cfg.binary_length:
                next_vert = get_upcoming_vertex(bin, next_chk_offset)
                chk.append(next_vert)
            cfg.offset = old_offset
            return chk, cache
        elif not get_chunk and cur_idx <= idx and idx <= cur_idx + deltas_in_chunk:
            # Found chunk containing vertex
            (x, y), cache = access_vertex_chk(bin, deltas_in_chunk_offset, delta_size, idx - cur_idx, cache)
            cfg.offset = old_offset
            return (x, y), cache

        # Jump to next chunk
        cur_idx += 1 + (deltas_in_chunk if not get_chunk else 0)
        cfg.offset += FLOAT_SIZE * 2 + delta_bytes_size
        chunks_in_ring_left -= 1
        if (chunks_in_ring_left == 0):
            rings_left -= 1
    raise Exception("Out of bounds!")

def get_upcoming_vertex(bin, chk_offset):
    cfg.offset = chk_offset + D_CNT_SIZE
    if cfg.COMPRESS_CHUNK or cfg.USE_ENTROPY:
        cfg.offset += D_BITSIZE_SIZE
   
    return (bytes_to_double(bin), bytes_to_double(bin))

#@profile
def access_vertex_chk(bin, chk_offset, delta_size, idx=None, cache=None, list_vertices=False):
    """
    Can be used if the chunk location in bin is already known, and a vertex within the chunk is needed.
    Supply the offset to D_CNT, and idx is the index within the chunk
    List_Vertices can be used to toggle between only returning the vertex at idx, or a list of all vertices
    within the chunk, up until and including idx.
    DO NOT USE ENTROPY OR CHUNK COMPRESSION WHEN CALLING METHOD.
    """
    old_offset = cfg.offset
    cfg.offset = chk_offset
    deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)
    if cfg.COMPRESS_CHUNK or cfg.USE_ENTROPY:
        delta_bytes_size = deltas_in_chunk * delta_size * 2 - bytes_to_int(bin, D_BITSIZE_SIZE)
    
    if idx == None:
        idx = deltas_in_chunk

    # Extract reset point
    x, y = (bytes_to_double(bin), bytes_to_double(bin))
    if list_vertices:
        vertices = [(x, y)]

    if cfg.COMPRESS_CHUNK:
        bin, _ = decompress_chunk(bin, cfg.offset, delta_bytes_size) 

    #Optimized for common case
    if not cfg.USE_ENTROPY and list_vertices:
        for idx in range(idx):
            x, y = bytes_to_decoded_coord_pair(bin, (x, y), delta_size)
            vertices.append((x, y))
    else:
        for idx in range(idx):
            x = bytes_to_decoded_coord(bin, x, delta_size)
            y = bytes_to_decoded_coord(bin, y, delta_size)
            if list_vertices:
                vertices.append((x, y))

    cfg.offset = old_offset
    return ((x, y) if not list_vertices else vertices), cache