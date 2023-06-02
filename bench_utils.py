import shapely
import json
import random
import pandas as pd
import geopandas as gpd
import glob
import tqdm
import pickle
import json

def load_shp_files(base_loc):
    """
        Can be used to read multiple .shp-files into df.
    """
    df = gpd.GeoDataFrame()
    files = glob.glob(base_loc + '/*.shp')
    for i in tqdm.tqdm(range(len(files))):
        f = files[i]
        print(i + 1, ':', f)
        to_append = gpd.read_file(f)
        print(i + 1, 'contains', len(to_append), 'geometries')
        df = gpd.pd.concat([df, to_append], copy=False)
    return df, list(range(len(df)))

def read_dataset(DATASET_PATH = "data/lund_building_highway.json", NBR_ITER = -1):
    """
        Path to dataset is either .json file or .shp file.
        Avoid loading large .shp files here. Instead use 'load_shp_files'.
    """
    #DATASET_PATH = "data/world.json"
    # Extract the nested feature attribute of the geo_json file containing the geometries
    if DATASET_PATH.endswith(".shp"):
        file_df = gpd.read_file(DATASET_PATH) # Parse .shp to GeoJson format
        geo = gpd.GeoSeries(file_df.geometry)
        geo_json = geo.to_json()
        file_df: pd.DataFrame = pd.json_normalize(json.loads(geo_json), record_path=['features'])
    else:
        with open(DATASET_PATH, 'r') as f:
            data = json.loads(f.read())
        file_df: pd.DataFrame = pd.json_normalize(data, record_path=['features'])
    if NBR_ITER != -1:
        file_df = file_df.head(NBR_ITER)
    # Create a dataframe suitable for the WKT format for easy convertion to shapely objects
    df = pd.DataFrame(
        {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

    if NBR_ITER != -1:
        max_idx = len(df) - 1
        unary_idxs = [random.randint(0, max_idx) for i in range(NBR_ITER)] # Generate list of indexes to query on
    else:
        unary_idxs = list(range(len(df)))
    return df, unary_idxs

def parse_intersection_data(file_name, max_shps=999999999, strip_precision=False):
    geom_pairs = []
    geom_stats = []
    if file_name.endswith('.json'):
        with open(f'data/intersection/{file_name}','r') as file:
            data = json.loads(file.read())
        for i in range(0, min(len(data), max_shps * 4), 4): # Don't include a new line in end of file
            type = data[i]
            p1 = shapely.from_wkt(data[i + 1])
            p2 = shapely.from_wkt(data[i + 2])
            if strip_precision:
                p1 = shapely.from_wkt(shapely.to_wkt(p1, rounding_precision=7))
                p2 = shapely.from_wkt(shapely.to_wkt(p2, rounding_precision=7))
            intersects = data[i + 3]
            geom_pairs.append((p1, p2))
            geom_stats.append((type, intersects))
    elif file_name.endswith('.pkl'):
        with open(f'data/intersection/{file_name}', 'rb') as f:
            intersections = pickle.load(f)
        for type, p1_wkb, p2_wkb, intersects in intersections: # Don't include a new line in end of file
            p1 = shapely.from_wkb(p1_wkb)
            p2 = shapely.from_wkb(p2_wkb)
            if strip_precision:
                p1 = shapely.from_wkt(shapely.to_wkt(p1, rounding_precision=7))
                p2 = shapely.from_wkt(shapely.to_wkt(p2, rounding_precision=7))
            geom_pairs.append((p1, p2))
            geom_stats.append((type, intersects))
    else:
        file = open(f'data/intersection/{file_name}', 'r')
        lines = file.read().splitlines()
        for i in range(0, min(len(lines), max_shps), 3): # Don't include a new line in end of file
            p1 = shapely.from_wkt(lines[i])
            p2 = shapely.from_wkt(lines[i + 1])
            geom_pairs.append((p1, p2))
    return geom_pairs, geom_stats