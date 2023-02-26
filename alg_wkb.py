from base import CompressionAlgorithm
from shapely.geometry import shape
import shapely.wkt
import shapely
import time


class Wkb(CompressionAlgorithm):

    def compress(self, geometry):
        s = time.perf_counter()
        content = shapely.to_wkb(geometry)
        t = time.perf_counter()
        return t - s, content
    

    def decompress(self, bin):
        s = time.perf_counter()
        geometry = shapely.from_wkb(bin)
        t = time.perf_counter()
        return t - s, geometry     


# ---- UNARY ---- #

    def vertices(self, bin):
        s = time.perf_counter()
        geometry = shapely.wkb.loads(bin)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords

    def type(self, bin): 
        s = time.perf_counter()
        geometry = shapely.wkb.loads(bin)
        type = geometry.geom_type
        t = time.perf_counter()
        return t - s, type
    
    def bounding_box(self, bin):
        s = time.perf_counter()
        geometry = shapely.wkb.loads(bin)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds
    
    def add_vertex(self, args):
        bin, idx, pos = args
        s = time.perf_counter()
        geometry = shapely.wkb.loads(bin)
        _, bin = self.compress(geometry)
        t = time.perf_counter()
        return t - s, bin
    
    def is_intersecting(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        _, r_geo = self.decompress(r_bin)
        res = shapely.intersects(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, res

    def intersection(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        _, r_geo = self.decompress(r_bin)
        res = shapely.intersection(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, res

