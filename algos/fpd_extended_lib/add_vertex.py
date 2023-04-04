import time
from shapely import GeometryType as GT
import bisect
from bitarray import bitarray


class AddVertex:
    ALG = None
    funcs = None

    def __init__(self, ALG) -> None:
        self.ALG = ALG

    def add_vertex(self, args):
        def create_chunk(reset_point, delta_cnt=0):
            middle = bitarray(endian='big')
            middle.extend(self.ALG.uint_to_ba(delta_cnt, self.ALG.D_CNT_SIZE))
            # Add full coordinates
            middle.frombytes(self.ALG.double_to_bytes(reset_point[0]))
            middle.frombytes(self.ALG.double_to_bytes(reset_point[1]))
            return middle

        bin_in, insert_idx, pos = args
        s = time.perf_counter()

        self.ALG.offset = 0
        bin = bitarray(endian='big')
        bin.frombytes(bin_in)

        delta_size, type = self.ALG.decode_header(bin)
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
                rings_left = self.ALG.bytes_to_uint(bin, self.ALG.POLY_RING_CNT_SIZE)
            if not is_linestring and chunks_in_ring_left == 0:
                chunks_in_ring_offset = self.ALG.offset
                chunks_in_ring_left = self.ALG.bytes_to_uint(bin, self.ALG.RING_CHK_CNT_SIZE)
                chunks_in_ring = chunks_in_ring_left
            deltas_in_chunk_offset = self.ALG.offset
            deltas_in_chunk = self.ALG.bytes_to_uint(bin, self.ALG.D_CNT_SIZE)

            # print(p_idx, deltas_in_chunk, insert_idx)
            # Found chunk to append/prepend?
            if p_idx <= insert_idx and insert_idx <= p_idx + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0):
                deltas_left = max(insert_idx - p_idx - 1, 0)
                deltas_right = deltas_in_chunk - deltas_left
                # print(deltas_left, deltas_right)
                # Handle left
                if p_idx != insert_idx:  # Has a left coordinate?
                    split = deltas_in_chunk_offset + self.ALG.D_CNT_SIZE + self.ALG.FLOAT_SIZE * 2 + delta_size * 2 * deltas_left
                    # Update delta cnt
                    bin[deltas_in_chunk_offset:deltas_in_chunk_offset + self.ALG.D_CNT_SIZE] = self.ALG.uint_to_ba(deltas_left, self.ALG.D_CNT_SIZE)
                else:
                    split = deltas_in_chunk_offset

                middle = create_chunk(pos)
                chunks_in_ring += 1

                # Handle chunk tail
                if deltas_right > 0 and p_idx != insert_idx:
                    # Get the absolute coordinate for first right coordinate
                    rst_p, _ = self.funcs.access_vertex_chk(bin, deltas_in_chunk_offset, insert_idx - p_idx, delta_size)
                    right = create_chunk(rst_p, deltas_right - 1)
                    # Append old tail, without the one extracted point
                    right.extend(bin[split + delta_size * 2:])
                    chunks_in_ring += 1
                else:
                    right = bin[split:]

                left = bin[0:split]
                if not is_linestring:
                    left[chunks_in_ring_offset:chunks_in_ring_offset + self.ALG.RING_CHK_CNT_SIZE] = self.ALG.uint_to_ba(chunks_in_ring, self.ALG.RING_CHK_CNT_SIZE)
                bin = left
                bin.extend(middle)
                bin.extend(right)

                break
            else:
                # Jump to next chunk
                p_idx += 1 + deltas_in_chunk + (1 if chunks_in_ring_left == 1 else 0)
                self.ALG.offset += self.ALG.FLOAT_SIZE * 2 + delta_size * 2 * deltas_in_chunk
                chunks_in_ring_left -= 1
                if (chunks_in_ring_left == 0):
                    rings_left -= 1

                if self.ALG.offset >= bin_len and is_linestring:
                    # Reached end without appending: is linestring!
                    new = create_chunk(pos)
                    bin.extend(new)
                    break

        bin = bin.tobytes()
        t = time.perf_counter()
        return t - s, bin
