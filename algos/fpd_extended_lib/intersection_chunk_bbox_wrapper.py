
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import math
import numpy as np
from shapely import GeometryType as GT
from bitarray import bitarray, util, bits2bytes


def append_intersection_header(self, bits, geometry):
    pass


def is_intersecting(self, args):
    l_bin, r_bin = args
    s = time.perf_counter()

    _, l_geo = self.ALG.decompress(l_bin)
    _, r_geo = self.ALG.decompress(r_bin)
    res = shapely.intersects(l_geo, r_geo)

    t = time.perf_counter()
    return t - s, res


def intersection(self, args):
    l_bin, r_bin = args
    s = time.perf_counter()
    _, l_geo = self.ALG.decompress(l_bin)
    _, r_geo = self.ALG.decompress(r_bin)
    res = shapely.intersection(l_geo, r_geo)
    t = time.perf_counter()

    return t - s, res
