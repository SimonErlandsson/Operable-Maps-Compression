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
        _, geometry = self.decompress(bin)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords

    def type(self, bin): 
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        type = geometry.geom_type
        t = time.perf_counter()
        return t - s, type
    
    def bounding_box(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds
    
    def add_vertex(self, args):
        bin, insert_idx, pos = args
        s = time.perf_counter()

        _, geometry = self.decompress(bin)
        wkt = shapely.to_wkt(geometry, rounding_precision=-1)
        point_idx = 0
        prev = ''
        for c_idx, char in enumerate(wkt):
            if char == ',' and prev != ')':
                if insert_idx == point_idx:
                    insert_string = f', {pos[0]} {pos[1]}'
                    wkt = wkt[:c_idx] + insert_string + wkt[c_idx:]
                    break
                point_idx += 1
            prev = char
        geometry = shapely.wkt.loads(wkt)
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


def main():
    x = Wkb()
    #geom = shapely.wkt.loads("POLYGON ((13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362))")
    geom = shapely.wkt.loads("LINESTRING (13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362)")
    t, bin = x.compress(geom)
    _, v = x.vertices(bin)
    print(v)
    _, add_bin = x.add_vertex((bin, 0, (0.2, 0.3)))
    _, v = x.vertices(add_bin)
    print(v)

if __name__ == "__main__":
    main()