import json
import os
import shutil
from base import CompressionAlgorithm
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape
import shapely.wkt
import linecache
import gzip
import time


class SetAsideCompression(CompressionAlgorithm):
    # Fields for saving most recent compressed/decompressed file
    num_to_type = {}
    max_geom_idx = 0


    #Expects a geoJSON file
    def compress(self, file_comp, geometry):
        #Convert the geometry category names to numbers for smaller space (COMMENTED OUT FOR NOW) <- REF 1
        
        # type_comp = pd.factorize(file_df['geometry.type'])
        # self.num_to_type = dict({(idx, type_name) for idx, type_name in enumerate(type_comp[1])})
        # self.num_to_type.update({(type_name, idx) for idx, type_name in enumerate(type_comp[1])})

        #Create pre computed values to store as metadata
        geo_type = geometry.geom_type #self.num_to_type[geometry.geom_type] <- FOR REF 1
        geo_vert_count = shapely.count_coordinates(geometry)
        geo_area = shapely.area(geometry)
        geo_length = shapely.length(geometry)

        #Write the metadata data for the operations as a special file per geometry
        f = open(file_comp, "w")
        f.write(str(geo_type) + "\t" + str(geo_vert_count) + "\t" + str(geo_area) + "\t" + str(geo_length) + "\n")

        #Write the compressed geometry
        compressed_geometry = gzip.compress(bytes(str(geometry), 'utf-8'))
        f = open(file_comp, "ab")
        f.write(compressed_geometry)
        f.close() 
             
        

    def decompress(self, file_comp, file_decomp):
        #Extract the number of files contained in the compression map
        nbr_of_geoms = len([entry for entry in os.listdir(file_comp) if os.path.isfile(os.path.join(file_comp, entry))])
        feauture_list = []
        
        #Iterate through each file, decompress geometry and append to decompression structure
        for idx in range(nbr_of_geoms):
            f = open('%s/%s.txt'%(file_comp, idx), "rb")
            next(f)  #Skip first line containing operation information
            decompressed_data = gzip.decompress(f.read()).decode('utf-8') #Decompressing data
            f.close()

            #Add decompressed data as JSON structure to the feature list containing all decompressed geometries
            feauture_list.append({'type': 'Feature', 'properties': {
            }, 'geometry': shapely.geometry.mapping(shapely.wkt.loads(decompressed_data))})
        
        # Create the FeatureCollection wrapper standard to GeoJson
        geojson_dict = {"type": "FeatureCollection", "features": feauture_list}
       
        # Write to file
        with open(file_decomp, "w") as file:
            json.dump(geojson_dict, file)

        self.file_decomp = file_decomp




# ---- UNARY ---- #

    def vertices(self, bin_file):
        s = time.perf_counter()
        data = data.split(b'\n', 1)[1] # Split at newline byte
        decompressed_data = gzip.decompress(data).decode('utf-8')
        geometry = shapely.wkt.loads(decompressed_data)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords


    def type(self, idx):
        idx_file = open('%s/%s.txt'%(self.file_comp_path, idx), "rb")
        data = idx_file.readline().decode('utf-8').split('\t')
        return data[0]
    
    def vertex_count(self, idx):
        idx_file = open('%s/%s.txt'%(self.file_comp_path, idx), "rb")
        data = idx_file.readline().decode('utf-8').split('\t')
        return data[1]

    # For Polygon
    def area(self, idx):
        idx_file = open('%s/%s.txt'%(self.file_comp_path, idx), "rb")
        data = idx_file.readline().decode('utf-8').split('\t')
        return data[2]

    def length(self, idx):
        idx_file = open('%s/%s.txt'%(self.file_comp_path, idx), "rb")
        data = idx_file.readline().decode('utf-8').split('\t')
        return data[3]


def main():
    x = SetAsideCompression()
    x.compress('data/lund_building_highway.json', 'data/meta_aside_output')
    x.decompress('data/meta_aside_output', 'data/testbench_decompressed')
    print(x.vertices(2404))


if __name__ == "__main__":
    main()
