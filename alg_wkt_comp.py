import json
import shutil
from base import CompressionAlgorithm
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape
import shapely.wkt
import shapely
import linecache
import gzip
import shutil
import os



class WktComp(CompressionAlgorithm):
    # Fields for saving most recent compressed/decompressed file
    file_comp = None
    file_decomp = None

    def compress(self, file_uncomp, file_comp):
        # Extract the nested feature attribute of the geo_json file containing the geometries
        with open(file_uncomp, 'r') as f:
            data = json.loads(f.read())
        file_df: pd.DataFrame = pd.json_normalize(
            data, record_path=['features'])
        df = pd.DataFrame(
            {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

        # Fill an array of all the geometries. Done for only writing once to the list
        rows = []
        for _, row in df.iterrows():
            rows.append(str(shape(row).wkt) + '\n')

        # Write the result to file
        f = open(file_comp + '_temp', 'w')
        f.writelines(rows)
        f.close()

        with open(file_comp + '_temp', 'rb') as f_in:
            with gzip.open(file_comp, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_comp + '_temp')

        self.file_comp = file_comp

    def decompress(self, file_comp, file_decomp):
        # Extract the file in WKT format where each line is a geometry
        f = gzip.open(file_comp, 'rb')
        lines = f.readlines()
        # For each geometry convert to corresponding json format
        feauture_list = []
        for line in lines:
            feauture_list.append({'type': 'Feature', 'properties': {
            }, 'geometry': shapely.geometry.mapping(shapely.wkt.loads(line))})

        # Create the FeatureCollection wrapper standard to GeoJson
        geojson_dict = {"type": "FeatureCollection", "features": feauture_list}
        # Write to file
        with open(file_decomp, "w") as file:
            json.dump(geojson_dict, file)

        self.file_decomp = file_decomp

    def geometry_count(self):
        f = gzip.open(self.file_comp, 'rb')
        lines = f.readlines()
        return len(lines)


# ---- UNARY ---- #
    # Note: getline is not zero-index based

    def vertices(self, idx):
        return self.uncompress_and_apply(idx, self.vertices_helper)


    def type(self, idx):
        return self.uncompress_and_apply(idx, self.type_helper)


    def vertex_count(self, idx):
        return self.uncompress_and_apply(idx, self.vertex_count_helper)


    # For Polygon
    def area(self, idx):
        return self.uncompress_and_apply(idx, self.area_helper)


    def length(self, idx):
        return self.uncompress_and_apply(idx, self.length_helper)


# ---- HELPER METHODS ---- #
    def vertices_helper(self, idx, file):
        line_at_idx = linecache.getline(file, idx + 1)
        geometry = shapely.wkt.loads(line_at_idx)
        return shapely.get_coordinates(geometry)

    def type_helper(self, idx, file):
        
        line_at_idx = linecache.getline(file, idx + 1)
        geometry = shapely.wkt.loads(line_at_idx)
        return geometry.geom_type

    def vertex_count_helper(self, idx,file):
        line_at_idx = linecache.getline(file, idx + 1)
        geometry = shapely.wkt.loads(line_at_idx)
        return shapely.count_coordinates(geometry)

    # For Polygon
    def area_helper(self, idx, file):
        line_at_idx = linecache.getline(file, idx + 1)
        geometry = shapely.wkt.loads(line_at_idx)
        return shapely.area(geometry)

    def length_helper(self, idx, file):
        line_at_idx = linecache.getline(file, idx + 1)
        geometry = shapely.wkt.loads(line_at_idx)
        return shapely.length(geometry)


    def uncompress_and_apply(self, idx, task):
        with gzip.open(self.file_comp, 'rb') as f_in:
            with open(self.file_comp + '_temp', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        res = task(idx, self.file_comp + '_temp')
        os.remove(self.file_comp + '_temp')
        return res



def main():
    x = WktComp()
    x.compress('data/lund_building_highway.json', 'data/testbench_compressed')
    x.decompress('data/testbench_compressed', 'data/testbench_decompressed')
    print(x.area(1))

if __name__ == "__main__":
    main()
