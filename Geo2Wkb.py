import json
import shutil
from base import CompressionAlgorithm
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape
import shapely.wkt
import linecache



class Geo2Wkb(CompressionAlgorithm):
    #Fields for saving most recent compressed/decompressed file
    file_comp = None
    file_decomp = None 

    def compress(self, file_uncomp, file_comp):
        #Extract the nested feature attribute of the geo_json file containing the geometries
        with open(file_uncomp,'r') as f:
            data = json.loads(f.read())        
        file_df:pd.DataFrame = pd.json_normalize(data, record_path=['features'])
        df = pd.DataFrame({'type':file_df['geometry.type'],'coordinates':file_df['geometry.coordinates']})
        
        #Fill an array of all the geometries. Done for only writing once to the list
        rows = []
        for _, row in df.iterrows():
            rows.append(str(shape(row).wkt) + '\n')
            
        #Write the result to file
        f = open(file_comp, 'w')
        f.writelines(rows)
        f.close()

        self.file_comp = file_comp
                

    def decompress(self, file_comp, file_decomp):
        #Extract the file in WKT format where each line is a geometry
        f = open(file_comp, 'r')
        lines = f.readlines()
        #For each geometry convert to corresponding json format
        feauture_list = []
        for line in lines:
            feauture_list.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(shapely.wkt.loads(line))})

        #Create the FeatureCollection wrapper standard to GeoJson
        geojson_dict = {"type": "FeatureCollection", "features": feauture_list}   
        #Write to file
        with open(file_decomp, "w") as file:
            json.dump(geojson_dict, file)
        
        self.file_decomp = file_decomp

        
    def geometry_count(self):
        return len(open(self.file_comp).readlines('\n'))


# ---- UNARY ---- #
    def vertices(self, idx):
        line_at_idx = linecache.getline(self.file_comp, idx)
        geometry = shapely.wkt.loads(line_at_idx)
        return list(geometry.coords)

    def type(self, idx):
        line_at_idx = linecache.getline(self.file_comp, idx)
        geometry = shapely.wkt.loads(line_at_idx)
        return geometry.geom_type

    def vertex_count(self, idx):
        return len(self.vertices(idx))
            
    # For Polygon
    def area(self, idx):
        line_at_idx = linecache.getline(self.file_comp, idx)
        geometry = shapely.wkt.loads(line_at_idx)
        return geometry.area

    def length(self, idx):
        line_at_idx = linecache.getline(self.file_comp, idx)
        geometry = shapely.wkt.loads(line_at_idx)
        return geometry.length
    


def main():
        x = Geo2Wkb()
        x.compress('data/lund_building_highway.json', 'data/testbench_compressed')
        x.decompress('data/testbench_compressed', 'data/testbench_decompressed')

        
if __name__ == "__main__":
        main()
