# Run main locally
import sys
from pathlib import Path  # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

from algos.base import CompressionAlgorithm
from algos.fpd_extended_lib.intersection_chunk_bbox_wrapper import *
from algos.fpd_extended_lib.add_vertex import *
from algos.fpd_extended_lib.functions import *
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
    access_vertex_chk = lambda self, bin, chk_offset, delta_size, idx, cache=None: access_vertex_chk(bin, chk_offset, delta_size, idx, cache)
    access_vertex = lambda self, bin, access_idx, cache=[]: access_vertex(bin, access_idx, cache)
    get_chunk = lambda self, bin, access_idx, cache=[]: access_chunk(bin, access_idx, cache)

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
    import random
    import json
    import pandas as pd
    import tqdm
    from shapely.geometry import shape
    from alg_fpd import Fpd
    x = FpdExtended()

    algos.fpd_extended_lib.entropy_coder.set_entropy_code("data/lund_building_highway.json")
    # geom1 = shapely.wkt.loads('POLYGON ((-24.3 10.48, -19.32 12.44, -15.3 14.2, -15.3 13.78, -15.3 13.9, -15.06 10.4, -17.44 11.38, -19.18 11.46, -14.82 9.08, -12.9 10.14, -12.08 7.86, -14.36 5.94, -15.92 8.34, -16.86 3.48, -19.38 4.4, -18.2 6.52, -20.08 7.4, -24.34 6.68, -24.24 8.66, -27.52 11.1,  -27.0 11.1, -24.3 10.48))')
    geom_fail = shapely.wkt.loads('POLYGON ((0 0, 10 10, 20 0, 0 0))')
    geom1 = shapely.wkt.loads('LINESTRING (-24.3 10.48, -19.32 12.44, -15.3 14.2)')
    geom2 = shapely.wkt.loads('POLYGON ((-9.9 16.85, -5.95 17.67, -6.19 13.49, -9.81 12.74, -7.35 9.2, -6.82 6.19, -10 6, -12.36 5.75, -14.59 8.1, -12 10, -13.93 12.31, -17.35 12.45, -16.83 15.6, -20.45 14.6, -22.36 12, -22 9.37, -27.1 6.48, -30 11.7, -27.9 15.5, -21.46 17.26, -19.6 16.1, -14.77 17.6, -11.43 13.32, -9.9 16.85))')

    geom3 = shapely.wkt.loads("MULTIPOLYGON (((13.193709 55.7021381, 13.1937743 55.7021279, 13.1938355 55.7021184, 13.1938461 55.702109, 13.1938566 55.7020984, 13.1938611 55.7020902, 13.1938655 55.7020774, 13.1938655 55.7020633, 13.1938583 55.7020408, 13.1938402 55.7020014, 13.1937184 55.7017259, 13.1937008 55.7016876, 13.1936836 55.7016654, 13.1936537 55.7016428, 13.1936223 55.7016242, 13.1935741 55.7016036, 13.1935354 55.7015911, 13.1935006 55.701584, 13.1934829 55.701598, 13.1934673 55.7016115, 13.1934736 55.7016164, 13.1934776 55.7016216, 13.1934875 55.7016633, 13.1934985 55.7016898, 13.1935196 55.7017337, 13.1935659 55.7018353, 13.1936162 55.7018282, 13.1936551 55.7019155, 13.1936651 55.7019377, 13.1936955 55.7020047, 13.1936497 55.7020119, 13.193709 55.7021381)), ((13.1938175 55.7017126, 13.1938602 55.7017068, 13.1939048 55.7017007, 13.1938998 55.7016861, 13.193892 55.7016685, 13.1938831 55.7016589, 13.193871 55.701651, 13.1938602 55.701646, 13.1938405 55.7016438, 13.193822 55.7016456, 13.1938062 55.7016517, 13.1937985 55.7016571, 13.1937953 55.7016646, 13.1937979 55.7016746, 13.1938017 55.7016836, 13.1938052 55.7016908, 13.1938175 55.7017126)), ((13.1940245 55.7019788, 13.19398 55.7019848, 13.1939372 55.7019907, 13.1939585 55.7020383, 13.1939692 55.7020479, 13.1939841 55.7020512, 13.1939975 55.7020519, 13.1940079 55.702051, 13.1940198 55.7020497, 13.1940317 55.7020463, 13.1940395 55.7020422, 13.1940435 55.7020369, 13.1940452 55.7020314, 13.1940457 55.7020218, 13.1940245 55.7019788)), ((13.1939779 55.7015541, 13.1939529 55.701555, 13.1939622 55.7015658, 13.1939755 55.7015942, 13.194075 55.7018201, 13.1941382 55.7019637, 13.1941483 55.7019866, 13.194164 55.7020087, 13.1941899 55.7020304, 13.1942142 55.7020424, 13.1942291 55.7020486, 13.1942638 55.702042, 13.195019 55.7018988, 13.1948681 55.7018923, 13.1944181 55.7018687, 13.1944172 55.7018717, 13.194395 55.7018706, 13.1942164 55.7018622, 13.194172 55.7017564, 13.1941218 55.701761, 13.1941279 55.7017262, 13.1941357 55.7016818, 13.1940872 55.7015737, 13.1940769 55.7015503, 13.1939779 55.7015541), (13.1942341 55.7020059, 13.1942075 55.7020095, 13.1941895 55.7019673, 13.1941696 55.701921, 13.1941936 55.7019177, 13.1941884 55.7019055, 13.19426 55.7018958, 13.1942645 55.7019063, 13.1943172 55.7018991, 13.1943567 55.7019912, 13.1943041 55.7019984, 13.1943086 55.7020089, 13.1942394 55.7020183, 13.1942341 55.7020059)))")

    geom4 = shapely.wkt.loads(
        "LINESTRING (13.199378 55.7034667, 13.1999441 55.7033986, 13.200125 55.7033882, 13.2002723 55.7033936, 13.2004383 55.7034097, 13.2005935 55.7034211, 13.2007699 55.703423, 13.2011275 55.7034136, 13.2012413 55.7034103, 13.2012947 55.7034088)")

    geom5 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848603 55.7056619, 13.1846238 55.7056422, 13.1848537 55.7057363))')
    geom8 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848861 55.705646, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1847300 55.705712, 13.1847425 55.7057123))')
    geom6 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848861 55.705646, 13.1848537 55.7057363))')
    
    geom7 = shapely.wkt.loads('MULTIPOLYGON (((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848861 55.705646, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1847300 55.705712, 13.1847425 55.7057123)), ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848841 55.705626, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847274 55.705711, 13.1848861 55.705646, 13.1847425 55.7057123)))')

    t, bin3 = x.compress(geom8)
    t, geomx = x.decompress(bin3)
    print(bin3)
    _, bin3 = x.add_vertex((bin3, 2, (13.1848861, 55.705646)))
    t, geomx = x.decompress(bin3)
    print(geomx)


if __name__ == "__main__":
    main()
