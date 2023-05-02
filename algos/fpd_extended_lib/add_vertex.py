import time
from shapely import GeometryType as GT
from bitarray import bitarray
from algos.fpd_extended_lib.helpers import *
from algos.fpd_extended_lib.cfg import *
import algos.fpd_extended_lib.cfg
from algos.fpd_extended_lib.low_level import *
from algos.fpd_extended_lib.decompress import *

def add_vertex(self, args):
    def create_chunk(reset_point, delta_cnt=0, delta_bytes=0):
        middle = bitarray(endian='big')
        middle.extend(uint_to_ba(delta_bytes, D_CNT_SIZE))
        middle.extend(uint_to_ba(delta_cnt, D_CNT_SIZE))
        # Add full coordinates
        middle.frombytes(double_to_bytes(reset_point[0]))
        middle.frombytes(double_to_bytes(reset_point[1]))
        return middle

    bin_in, insert_idx, pos = args
    s = time.perf_counter()

    cfg.offset = 0
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)

    delta_size, type = decode_header(bin)
    # Type specific variables
    is_linestring = type == GT.LINESTRING
    is_multipolygon = type == GT.MULTIPOLYGON

    p_idx = 0
    chunks_in_ring_left = 0  # Used for iteration
    chunks_in_ring = 0  # Cache for store later
    rings_left = 0
    cfg.binary_length = len(bin)
    while (p_idx <= insert_idx):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_offset = cfg.offset
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            chunks_in_ring = chunks_in_ring_left
        deltas_in_chunk_offset = cfg.offset
        deltas_bytes_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        # print(p_idx, deltas_in_chunk, insert_idx)
        # Found chunk to append/prepend?
        if p_idx <= insert_idx and insert_idx <= p_idx + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0):
            rst_p, _, delta_offsets = access_vertex_chk(bin, deltas_in_chunk_offset, delta_size, insert_idx - p_idx, get_delta_offsets=True)
            deltas_left = max(insert_idx - p_idx - 1, 0)
            deltas_right = deltas_in_chunk - deltas_left
            # print(deltas_left, deltas_right)
            # Handle left
            if p_idx != insert_idx:  # Has a left coordinate?
                split = deltas_in_chunk_offset + D_CNT_SIZE * 2 + FLOAT_SIZE * 2 + delta_offsets[deltas_left] - delta_offsets[0]
                # Update delta cnt
                bin[deltas_in_chunk_offset :deltas_in_chunk_offset + D_CNT_SIZE] = uint_to_ba(delta_offsets[deltas_left] - delta_offsets[0], D_CNT_SIZE)
                bin[deltas_in_chunk_offset + D_CNT_SIZE :deltas_in_chunk_offset + D_CNT_SIZE * 2] = uint_to_ba(deltas_left, D_CNT_SIZE)
                
            else:
                split = deltas_in_chunk_offset

            middle = create_chunk(pos)
            chunks_in_ring += 1

            # Handle chunk tail
            if deltas_right > 0 and p_idx != insert_idx:
                # Get the absolute coordinate for first right coordinate
                right = create_chunk(rst_p, deltas_right - 1, delta_bytes=delta_offsets[-1] - delta_offsets[deltas_left + 1])
                # Append old tail, without the one extracted point
                right.extend(bin[split + delta_offsets[deltas_left + 1] - delta_offsets[deltas_left]:])
                chunks_in_ring += 1
            else:
                right = bin[split:]

            left = bin[0:split]
            if not is_linestring:
                left[chunks_in_ring_offset:chunks_in_ring_offset + RING_CHK_CNT_SIZE] = uint_to_ba(chunks_in_ring, RING_CHK_CNT_SIZE)
            bin = left
            bin.extend(middle)
            bin.extend(right)

            break
        else:
            # Jump to next chunk
            p_idx += 1 + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0)
            cfg.offset += FLOAT_SIZE * 2 + deltas_bytes_in_chunk
            chunks_in_ring_left -= 1
            if (chunks_in_ring_left == 0):
                rings_left -= 1

            if cfg.offset >= cfg.binary_length and is_linestring:
                # Reached end without appending: is linestring!
                new = create_chunk(pos)
                bin.extend(new)
                break

    bin = bin.tobytes()
    t = time.perf_counter()
    return t - s, bin
