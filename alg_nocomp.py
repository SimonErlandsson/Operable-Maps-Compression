from base import CompressionAlgorithm
import shutil


import json
import geopandas as gpd


class NoCompression(CompressionAlgorithm):

    def compress(self, file_uncomp, file_comp):
        shutil.copyfile(file_uncomp, file_comp)
        self.file_comp = file_comp

    def decompress(self, file_comp, file_decomp):
        shutil.copyfile(file_comp, file_decomp)

    def geometry_count(self):
        with open(self.file_comp) as f:
            data = json.load(f)['features']
        return len(data)

# ---- UNARY ---- #
    def vertices(self, idx):
        with open(self.file_comp) as f:
            data = json.load(f)['features'][idx]['geometry']['coordinates']
        return data

    def type(self, idx):
        with open(self.file_comp) as f:
            data = json.load(f)['features'][idx]['geometry']['type']
        return data

    def vertex_count(self, idx):
        with open(self.file_comp) as f:
            data = json.load(f)['features'][idx]['geometry']['coordinates']
        return len(data)

    # For Polygon
    def area(self, idx):
        # with open(self.file_comp) as f:
        #    data = json.load(f)['features'][idx]
        # return gpd.GeoSeries.from_xy(data['geometry'])
        pass

    def length(self, idx):
        # with open(self.file_comp) as f:
        #    data = json.load(f)['features'][idx]['geometry']
        # return gpd.GeoSeries(geometry=data).length()
        pass
