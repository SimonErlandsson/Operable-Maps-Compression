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


class MetaWkt(CompressionAlgorithm):

    def compress(self, geometry):
        #Convert the geometry category names to numbers for smaller space (COMMENTED OUT FOR NOW) <- REF 1
        
        # type_comp = pd.factorize(file_df['geometry.type'])
        # self.num_to_type = dict({(idx, type_name) for idx, type_name in enumerate(type_comp[1])})
        # self.num_to_type.update({(type_name, idx) for idx, type_name in enumerate(type_comp[1])})

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        geo_type = geometry.geom_type #self.num_to_type[geometry.geom_type] <- FOR REF 1
        geo_vert_count = shapely.count_coordinates(geometry)
        geo_area = shapely.area(geometry)
        geo_length = shapely.length(geometry)

        #Write the metadata data for the operations as a special file per geometry
        content = (bytes(str(geo_type) + "\t" + str(geo_vert_count) + "\t" + str(geo_area) + "\t" + str(geo_length) + "\n", 'utf-8'))

        #Write the compressed geometry
        content += gzip.compress(bytes(shapely.to_wkt(geometry, rounding_precision=-1), 'utf-8'))
        t = time.perf_counter()
        return t - s, content
    

    def decompress(self, bin):
        s = time.perf_counter()
        data = bin.split(b'\n', 1)[1] # Split at newline byte
        decomp_geom = gzip.decompress(data).decode('utf-8') # Decompressing data
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
        data = str(bin.split(b'\t', 1)[0], 'utf-8') # Split at \t
        t = time.perf_counter()
        return t - s, data 
    
    def bounding_box(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds
    
    def add_vertex(self, args):
        bin, idx, pos = args
        s = time.perf_counter()
        data = bin.split(b'\n', 1)[1] # Split at newline byte
        decomp_geom = gzip.decompress(data).decode('utf-8')
        geometry = shapely.wkt.loads(decomp_geom)
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
    x = SetAsideCompression()
    geom = shapely.wkt.loads("POLYGON ((13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362))")
    t, bin = x.compress(geom)
    x.add_vertex((bin, [0,1], (24.5, 12.3)))
    print(geom)


if __name__ == "__main__":
    main()
