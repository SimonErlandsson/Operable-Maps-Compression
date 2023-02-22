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


class SetAsideCompression(CompressionAlgorithm):
    # Fields for saving most recent compressed/decompressed file
    num_to_type = {}
    max_geom_idx = 0

    #Expects a geoJSON file
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
        content += gzip.compress(bytes(str(geometry), 'utf-8'))
        t = time.perf_counter()
        return t - s, content
    

    def decompress(self, bin):
        s = time.perf_counter()
        data = bin.split(b'\n', 1)[1] # Split at newline byte
        decomp_geom = gzip.decompress(data).decode('utf-8') # Decompressing data
        t = time.perf_counter()
        return t - s, decomp_geom     


# ---- UNARY ---- #

    def vertices(self, bin):
        s = time.perf_counter()
        data = bin.split(b'\n', 1)[1] # Split at newline byte
        decomp_geom = gzip.decompress(data).decode('utf-8')
        geometry = shapely.wkt.loads(decomp_geom)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords

    def type(self, bin): 
        s = time.perf_counter()
        data = bin.split(b'\t', 1)[0] # Split at \t
        t = time.perf_counter()
        return t - s, data 
    
    def bounding_box(self, bin):
        s = time.perf_counter()
        data = bin.split(b'\n', 1)[1] # Split at newline byte
        decomp_geom = gzip.decompress(data).decode('utf-8')
        geometry = shapely.wkt.loads(decomp_geom)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds
    
    def add_vertex(self, args):
        bin, idx, pos = args
        s = time.perf_counter()
        data = bin.split(b'\n', 1)[1] # Split at newline byte
        decomp_geom = gzip.decompress(data).decode('utf-8')
        geometry = shapely.wkt.loads(decomp_geom)
        coords = shapely.get_coordinates(geometry)
        print(len(idx))
        print(coords)
        # Fix add
        _, bin = self.compress(geometry)
        t = time.perf_counter()
        return t - s, bin
    
    def is_intersecting(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        l_geo = shapely.wkt.loads(l_geo)

        _, r_geo = self.decompress(r_bin)
        r_geo = shapely.wkt.loads(r_geo)
        shapely.intersects(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, True

    def intersection(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        l_geo = shapely.wkt.loads(l_geo)

        _, r_geo = self.decompress(r_bin)
        r_geo = shapely.wkt.loads(r_geo)
        shapely.intersection(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, []


def main():
    x = SetAsideCompression()
    x.compress('data/lund_building_highway.json', 'data/meta_aside_output')
    x.decompress('data/meta_aside_output', 'data/testbench_decompressed')
    print(x.vertices(2404))


if __name__ == "__main__":
    main()
