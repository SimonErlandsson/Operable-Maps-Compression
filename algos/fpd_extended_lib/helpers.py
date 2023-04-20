import math
import zlib
import gzip
from bitarray import bitarray, util
from shapely import GeometryType as GT
from algos.fpd_extended_lib.cfg import *
from algos.fpd_extended_lib.low_level import *
from algos.fpd_extended_lib.decompress import *


def get_chunks(bin_in, include_ring_start=True):
    """
    Only used for non-timed operations, i.e. debugging/testing implementations
    # NOTE: Can also return the first ring coordinate when reaching end of ring!
    """
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
        deltas_bytes_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)
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
            if include_ring_start:
                chunks[-1].append([x_ring, y_ring])
            rings_left -= 1
            is_last_ring_chunk[-1] = True
    return chunks, is_last_ring_chunk


def access_vertex(bin_in, idx, cache=[]):
    return random_access(bin_in, idx, cache, get_chunk=False)


def access_chunk(bin_in, idx, cache=[]):
    """
    Get a chunk based on chunk index. Note that
    chunks are connected. I.e. each chunk contains the
    next point in sequence.
    """
    return random_access(bin_in, idx, cache, get_chunk=True)


# Supply cache if used repeatedly for same shape
def random_access(bin_in, idx, cache, get_chunk=False):
    old_offset = cfg.offset
    cfg.offset = 0
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)

    delta_size, type = decode_header(bin)

    # Type specific variables
    is_linestring = type == GT.LINESTRING
    is_multipolygon = type == GT.MULTIPOLYGON

    cur_idx = 0
    chunks_in_ring_left = 0  # Used for iteration
    rings_left = 0
    while (cur_idx <= idx):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            ring_start_offset = cfg.offset
        deltas_in_chunk_offset = cfg.offset
        deltas_bytes_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        if get_chunk and cur_idx == idx:
            # Looking for whole chunk, found it
            chk, cache = access_vertex_chk(bin, deltas_in_chunk_offset, delta_size, cache=cache, list_vertices=True)

            # Chunks should contain next point for overlapping bboxes
            if not is_linestring:
                if chunks_in_ring_left == 1:
                    # Was last chunk in ring, append start
                    next_chk_offset = ring_start_offset
                else:
                    # Append next chunk start
                    next_chk_offset = cfg.offset + FLOAT_SIZE * 2 + delta_size * 2 * deltas_in_chunk
                next_vert, cache = access_vertex_chk(bin, next_chk_offset, delta_size, 0, cache)
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
        cfg.offset += FLOAT_SIZE * 2 + deltas_bytes_in_chunk
        chunks_in_ring_left -= 1
        if (chunks_in_ring_left == 0):
            rings_left -= 1
    raise Exception("Out of bounds!")


def access_vertex_chk(bin, chk_offset, delta_size, idx=None, cache=None, list_vertices=False, get_delta_offsets=False):
    """
    Can be used if the chunk location in bin is already known, and a vertex within the chunk is needed.
    Supply the offset to D_CNT, and idx is the index within the chunk
    List_Vertices can be used to toggle between only returning the vertex at idx, or a list of all vertices
    within the chunk, up until and including idx.
    """
    old_offset = cfg.offset
    cfg.offset = chk_offset + D_CNT_SIZE #skips delta bytes
    deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)
    if idx == None:
        idx = deltas_in_chunk

    # Extract reset point
    x, y = (bytes_to_double(bin), bytes_to_double(bin))
    if list_vertices:
        vertices = [(x, y)]
    delta_offsets = []
    # Loop through deltas in chunk
    for idx in range(idx):
        delta_offsets.append(cfg.offset)
        x = bytes_to_decoded_coord(bin, x, delta_size)
        y = bytes_to_decoded_coord(bin, y, delta_size)
        if list_vertices:
            vertices.append((x, y))
    delta_offsets.append(cfg.offset)

    cfg.offset = old_offset
    if get_delta_offsets:
        return ((x, y) if not list_vertices else vertices), cache, delta_offsets
    
    return ((x, y) if not list_vertices else vertices), cache

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

def get_zz_encoded_delta(prev_coord, curr_coord):
    return zz_encode(double_as_long(curr_coord) - double_as_long(prev_coord))


def compress_chunk(bits, chk_hdr_offset, delta_bytes_size):
    chk_coord_offset = chk_hdr_offset + D_CNT_SIZE * 2
    before_bits = bits[:chk_coord_offset]
    coords_bits = bits[chk_coord_offset: chk_coord_offset + delta_bytes_size + 2 * FLOAT_SIZE]
    if cfg.COMPRESSION_METHOD == "zlib":
        comp_coords_bytes = zlib.compress(coords_bits.tobytes())
    elif cfg.COMPRESSION_METHOD == "gzip":
        comp_coords_bytes = gzip.compress(coords_bits.tobytes())
    else:
        comp_coords_bytes = coords_bits.tobytes()
    comp_coords_bits = bitarray()
    comp_coords_bits.frombytes(comp_coords_bytes)
    before_bits += comp_coords_bits
    return before_bits, len(comp_coords_bits)

def decompress_chunk(bits, chk_coord_offset, coords_bytes_size):
    before_bits = bits[:chk_coord_offset]
    comp_coords_bits = bits[chk_coord_offset: chk_coord_offset + coords_bytes_size]
    after_bits = bits[chk_coord_offset + coords_bytes_size:]

    if cfg.COMPRESSION_METHOD == "zlib":
        coords_bits = zlib.decompress(comp_coords_bits.tobytes())
    elif cfg.COMPRESSION_METHOD == "gzip":
        coords_bits = gzip.decompress(comp_coords_bits.tobytes())
    else:
        coords_bits = comp_coords_bits.tobytes()
    decompressed_bits = bitarray()
    decompressed_bits.frombytes(coords_bits)

    before_bits[chk_coord_offset - D_CNT_SIZE * 2:chk_coord_offset - D_CNT_SIZE] =  uint_to_ba(len(decompressed_bits) - 2 * FLOAT_SIZE, D_CNT_SIZE) 

    before_bits += decompressed_bits
    before_bits += after_bits
    cfg.binary_length = len(before_bits)
    return before_bits, len(decompressed_bits)
