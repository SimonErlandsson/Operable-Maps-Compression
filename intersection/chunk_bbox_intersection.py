import shapely
import numpy as np
from algos.alg_fpd_extended import FpdExtended

fpd = FpdExtended()

def common_bbox(bins): # Returns lower left corner, upper right corner in 1D array
    x_l_1, y_b_1, x_r_1, y_t_1 = fpd.bounding_box(bins[0])[1] # Left, bottom, right, top
    x_l_2, y_b_2, x_r_2, y_t_2 = fpd.bounding_box(bins[1])[1]
    bbox = [max(x_l_1, x_l_2), max(y_b_1, y_b_2), min(x_r_1, x_r_2), min(y_t_1, y_t_2)]
    x_l, y_b, x_r, y_t = bbox
    return bbox if x_r >= x_l and y_t >= y_b else None

def get_chunk_bounds(bin):
    pass

def is_intersecting(bins):
    bbox = common_bbox(bins)
    if bbox == None:
        return False
    
    # Bounding boxes intersect. Assume intersection, ensure that no intersection is in fact occuring: 
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box
    return True