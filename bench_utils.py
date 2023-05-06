import shapely
import json
import random
import pandas as pd
import json

def read_dataset(DATASET_PATH = "data/lund_building_highway.json", NBR_ITER = -1):
    #DATASET_PATH = "data/world.json"
    # Extract the nested feature attribute of the geo_json file containing the geometries
    with open(DATASET_PATH, 'r') as f:
        data = json.loads(f.read())
    file_df: pd.DataFrame = pd.json_normalize(data, record_path=['features'])
    # Create a dataframe suitable for the WKT format for easy convertion to shapely objects
    df = pd.DataFrame(
        {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

    if NBR_ITER != -1:
        max_idx = len(df) - 1
        unary_idxs = [random.randint(0, max_idx) for i in range(NBR_ITER)] # Generate list of indexes to query on
    else:
        unary_idxs = list(range(len(df)))
    return df, unary_idxs

def parse_intersection_data(file_name, max_shps=999999999):
    geom_pairs = []
    geom_stats = []
    if file_name.endswith('.json'):
        with open(f'data/intersection/{file_name}','r') as file:
            data = json.loads(file.read())
        for i in range(0, min(len(data), max_shps * 4), 4): # Don't include a new line in end of file
            type = data[i]
            p1 = shapely.from_wkt(data[i + 1])
            p2 = shapely.from_wkt(data[i + 2])
            intersects = data[i + 3]
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