import json
import os
import shutil
from base import CompressionAlgorithm
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape
import shapely.wkt
import shapely
import linecache
import gzip
import time
import numpy as np


class WktComp(CompressionAlgorithm):

    def compress(self, geometry):
        #Convert the geometry category names to numbers for smaller space (COMMENTED OUT FOR NOW) <- REF 1
        
        # type_comp = pd.factorize(file_df['geometry.type'])
        # self.num_to_type = dict({(idx, type_name) for idx, type_name in enumerate(type_comp[1])})
        # self.num_to_type.update({(type_name, idx) for idx, type_name in enumerate(type_comp[1])})

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        
        #Write the compressed geometry
        content = gzip.compress(bytes(shapely.to_wkt(geometry, rounding_precision=-1), 'utf-8'))
        t = time.perf_counter()
        return t - s, content
    

    def decompress(self, bin):
        s = time.perf_counter()
        decomp_geom = gzip.decompress(bin).decode('utf-8') # Decompressing data
        geometry = shapely.wkt.loads(decomp_geom)
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

        wkt = gzip.decompress(bin).decode('utf-8') # Decompressing data
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
        bin = gzip.compress(bytes(wkt, 'utf-8')) 

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
    import random
    import matplotlib.pyplot as plt

    x = WktComp()
    #geom = shapely.wkt.loads("POLYGON ((13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362))")
    geom = shapely.wkt.loads("LINESTRING (13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362)")
    t, bin = x.compress(geom)
    _, v = x.vertices(bin)
    print(v)
    #add_idx = random.randint(0, len(v) - 1)
    #add_point = (v[add_idx][0] + random.randint(-25, 25) * 0.00001, v[add_idx][1] + random.randint(-25, 25) * 0.00001)
    #_, add_bin = x.add_vertex((bin, add_idx, (add_point)))
    _, add_bin = x.add_vertex((bin, 2, (0.2, 0.3)))
    # print(add_point, add_idx)
    # #t, bin = x.add_vertex((bin, 3, (24.5, 12.3)))
    # t, add_geom = x.decompress(add_bin)
    #wkt = shapely.to_wkt(geom)
    # print(wkt)
    _, v = x.vertices(add_bin)
    print(v)
    # _, v = x.vertices(add_bin)

    # p = gpd.GeoSeries(geom)
    # p.plot()
    # p = gpd.GeoSeries(add_geom)
    # p.plot()
    # plt.show()
      

if __name__ == "__main__":
    main()
