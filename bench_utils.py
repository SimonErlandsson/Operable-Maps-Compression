def read_dataset(DATASET_PATH = "data/lund_building_highway.json"):
    import random
    import pandas as pd
    import json
    #DATASET_PATH = "data/world.json"
    NBR_ITER = 16000


    SEED = 123 # If we want to enforce the same ordering and indexes for multiple runs, else None
    random.seed(SEED) # Init random

    # Extract the nested feature attribute of the geo_json file containing the geometries
    with open(DATASET_PATH, 'r') as f:
        data = json.loads(f.read())
    file_df: pd.DataFrame = pd.json_normalize(data, record_path=['features'])
    # Create a dataframe suitable for the WKT format for easy convertion to shapely objects
    df = pd.DataFrame(
        {'type': file_df['geometry.type'], 'coordinates': file_df['geometry.coordinates']})

    max_idx = len(df) - 1
    unary_idxs = [random.randint(0, max_idx) for i in range(NBR_ITER)] # Generate list of indexes to query on
    random.seed(SEED) # Reset random
    return df, unary_idxs

