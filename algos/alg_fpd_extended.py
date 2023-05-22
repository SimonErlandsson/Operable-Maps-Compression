# Run main locally
import sys
from pathlib import Path  # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

from algos.base import CompressionAlgorithm
from algos.fpd_extended_lib.intersection_chunk_bbox_wrapper import *
from algos.fpd_extended_lib.add_vertex import *
from algos.fpd_extended_lib.operations import *
from algos.fpd_extended_lib.helpers import *
from algos.fpd_extended_lib.compress import *
from algos.fpd_extended_lib.decompress import *

from collections import deque
import shapely
import shapely.wkt
from bitarray import bitarray, util, bits2bytes

#####                 #####
# --- FP-DELTA BASELINE ---#
#####                 #####

# char: 8 bits
# float: 32 bits
# double: 64 bits
# int: 32 bits
# long: 64 bits


class FpdExtended(CompressionAlgorithm):

    # ---- HELPER METHODS

    # Export some helper functions
    get_chunks = lambda self, bin, include_ring_start=True: get_chunks(bin, include_ring_start)
    access_vertex_chk = lambda self, bin, chk_offset, delta_size, idx, cache=None: access_vertex_chk(bin, chk_offset, delta_size, idx, cache=cache)
    access_vertex = lambda self, bin, access_idx, cache=[]: access_vertex(bin, access_idx, cache)
    get_chunk = lambda self, bin, access_idx, offset_cache=None, cache=[]: access_chunk(bin, access_idx, offset_cache=offset_cache, cache=cache)

    # Intersection
    get_chunk_bounds = lambda self, bin: get_chunk_bounds(bin)

# ---- UNARY ---- #
    compress = compress
    decompress = decompress

# ---- UNARY ---- #
    vertices = vertices
    type = type
    bounding_box = bounding_box
    add_vertex = add_vertex

# ---- BINARY ---- #
    is_intersecting = is_intersecting
    intersection = intersection


def main():
    import numpy as np
    import shapely
    import bisect
    import shapely.wkt
    import matplotlib.pyplot as plt
    import math
    import geopandas as gpd
    import json
    from bench_utils import parse_intersection_data 
    import intersection.first_bin_search
    import intersection.chunk_bbox_intersection
    # binary_intersection = intersection.first_bin_search.binary_intersection
    # chunk_bbox_is_intersecting = intersection.chunk_bbox_intersection.is_intersecting
    chunk_bbox_intersection = intersection.chunk_bbox_intersection.intersection
    lund_data, lund_data_stats = parse_intersection_data("lund.json", 10000)
    world_data, world_data_stats = parse_intersection_data("world.json", 10)
    special_cases, _ = parse_intersection_data("latest_export.json", strip_precision=True)


    from algos.alg_fpd_extended import FpdExtended
    fpd = FpdExtended()
    total = 0
    passed = 0
    for g1, g2 in world_data:
        #is_intersecting, intersect_points = binary_intersection(g1, g2)

        _, b1 = fpd.compress(g1)
        _, b2 = fpd.compress(g2)
        exp_shape = shapely.intersection(g1, g2)
        stats, real_shape = chunk_bbox_intersection((b1, b2),get_stats=True)
        total += 1
        if exp_shape.equals(real_shape):
            passed += 1

        # if total % 100 == 0:
        #     print(f"Passed {passed} of {total}. Total in set: {len(data)}")
    print("World: ", total, passed)

    # total = 0
    # passed = 0
    # for g1, g2 in lund_data:
    #     #is_intersecting, intersect_points = binary_intersection(g1, g2)

    #     _, b1 = fpd.compress(g1)
    #     _, b2 = fpd.compress(g2)
    #     exp_shape = shapely.intersection(g1, g2)
    #     real_shape = fpd.intersection((b1, b2))

    #     total += 1
    #     if exp_shape.equals(real_shape[1]):
    #         passed += 1

    #     # if total % 100 == 0:
    #     #     print(f"Passed {passed} of {total}. Total in set: {len(data)}")
    # print("Lund: ", total, passed)

    # total = 0
    # passed = 0
    # for g1, g2 in special_cases:
    #     #is_intersecting, intersect_points = binary_intersection(g1, g2)

    #     _, b1 = fpd.compress(g1)
    #     _, b2 = fpd.compress(g2)
    #     exp_shape = shapely.intersection(g1, g2)
    #     real_shape = fpd.intersection((b1, b2))

    #     total += 1
    #     if exp_shape.equals(real_shape[1]):
    #         passed += 1

    #     # if total % 100 == 0:
    #     #     print(f"Passed {passed} of {total}. Total in set: {len(data)}")
    # print("Special_cases: ", total, passed)


if __name__ == "__main__":
    main()
