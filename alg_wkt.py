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


class Wkt(CompressionAlgorithm):

    def compress(self, geometry):
        s = time.perf_counter()
        #Write the compressed geometry
        content = bytes(shapely.to_wkt(geometry, rounding_precision=-1), 'utf-8')
        t = time.perf_counter()
        return t - s, content
    

    def decompress(self, bin):
        s = time.perf_counter()
        decomp_geom = bin.decode('utf-8') # Is already in WKT
        geometry = shapely.from_wkt(decomp_geom)
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
        
        wkt = bin.decode('utf-8') # Decompressing data
        point_idx = 0
        for c_idx, char in enumerate(wkt):
            if char == ',':
                if insert_idx == point_idx:
                    insert_string = f', {pos[0]} {pos[1]}'
                    wkt = wkt[:c_idx] + insert_string + wkt[c_idx:]
                    break
                point_idx += 1
        bin = bytes(wkt, 'utf-8')

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
    x = Wkt()
    geom = shapely.wkt.loads("POLYGON ((13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362))")
    t, bin = x.compress(geom)
    x.add_vertex((bin, [0,1], (24.5, 12.3)))
    print(geom)


if __name__ == "__main__":
    main()
