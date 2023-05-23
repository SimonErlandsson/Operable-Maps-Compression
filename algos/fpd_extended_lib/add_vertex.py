import time
import bisect
from shapely import GeometryType as GT
from bitarray import bitarray
from algos.fpd_extended_lib.helpers import *
from algos.fpd_extended_lib.cfg import *
import algos.fpd_extended_lib.cfg
from algos.fpd_extended_lib.low_level import *

# Not implemented to work with COMPRESS_CHUNK==TRUE or USE_ENTROPY==True
def add_vertex(self, args):
    bin_in, insert_idx, pos = args
    if DISABLE_OPTIMIZED_ADD_VERTEX:
        from algos.alg_fpd_extended import decompress, compress
        s = time.perf_counter()

        _, geometry = decompress(None, bin_in)
        ragged = shapely.to_ragged_array([geometry])
        points = ragged[1]
        is_linestring = shapely.get_type_id(geometry) == shapely.GeometryType.LINESTRING

        # Use binary search O(log n) to find the index of the first element greater than insert_idx
        end_idx = min(len(ragged[2][0]) - 1, bisect.bisect_right(ragged[2][0], insert_idx))

        # Is the first coordinate in the ring?
        if ragged[2][0][end_idx - 1] == insert_idx and not is_linestring:
                points = np.delete(points, ragged[2][0][end_idx] - 1, axis=0)
        else:
            for i in range(end_idx, len(ragged[2][0])):
                    ragged[2][0][i] += 1
                    
        points = np.insert(points, insert_idx, pos, axis=0)        
 
        geometry = shapely.from_ragged_array(geometry_type=shapely.get_type_id(geometry), coords=points, offsets=ragged[2])[0]
        _, bin_in = compress(None, geometry)
        
        t = time.perf_counter()
        return t - s, bin_in

    """
    DO NOT USE ENTROPY OR CHUNK COMPRESSION WHEN CALLING METHOD.
    """
    assert(not(USE_ENTROPY or COMPRESS_CHUNK))
    def create_chunk(reset_point, delta_cnt=0):
        middle = bitarray(endian='big')
        middle.extend(uint_to_ba(delta_cnt, D_CNT_SIZE))
        # Add full coordinates
        middle.frombytes(double_to_bytes(reset_point[0]))
        middle.frombytes(double_to_bytes(reset_point[1]))
        return middle

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
    bin_len = len(bin)
    while (p_idx <= insert_idx):
        if is_multipolygon and rings_left == 0:
            rings_left = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
        if not is_linestring and chunks_in_ring_left == 0:
            chunks_in_ring_offset = cfg.offset
            chunks_in_ring_left = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
            chunks_in_ring = chunks_in_ring_left
        deltas_in_chunk_offset = cfg.offset
        deltas_in_chunk = bytes_to_uint(bin, D_CNT_SIZE)

        # print(p_idx, deltas_in_chunk, insert_idx)
        # Found chunk to append/prepend?
        if p_idx <= insert_idx and insert_idx <= p_idx + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0):
            deltas_left = max(insert_idx - p_idx - 1, 0)
            deltas_right = deltas_in_chunk - deltas_left
            # print(deltas_left, deltas_right)
            # Handle left
            if p_idx != insert_idx:  # Has a left coordinate?
                split = deltas_in_chunk_offset + D_CNT_SIZE + FLOAT_SIZE * 2 + delta_size * 2 * deltas_left
                # Update delta cnt
                bin[deltas_in_chunk_offset:deltas_in_chunk_offset + D_CNT_SIZE] = uint_to_ba(deltas_left, D_CNT_SIZE)
            else:
                split = deltas_in_chunk_offset

            middle = create_chunk(pos)
            chunks_in_ring += 1

            # Handle chunk tail
            if deltas_right > 0 and p_idx != insert_idx:
                # Get the absolute coordinate for first right coordinate
                rst_p, _ = access_vertex_chk(bin, deltas_in_chunk_offset, delta_size, insert_idx - p_idx)
                right = create_chunk(rst_p, deltas_right - 1)
                # Append old tail, without the one extracted point
                right.extend(bin[split + delta_size * 2:])
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
            cfg.offset += FLOAT_SIZE * 2 + delta_size * 2 * deltas_in_chunk
            chunks_in_ring_left -= 1
            if (chunks_in_ring_left == 0):
                rings_left -= 1

            if cfg.offset >= bin_len and is_linestring:
                # Reached end without appending: is linestring!
                new = create_chunk(pos)
                bin.extend(new)
                break

    bin = bin.tobytes()
    t = time.perf_counter()
    return t - s, bin
