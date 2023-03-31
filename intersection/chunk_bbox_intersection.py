import shapely
import numpy as np
from algos.alg_fpd_extended import FpdExtended
from intersection.plotting import *

fpd = FpdExtended()


def is_bbox_intersecting(bbox_1, bbox_2):
    x_l_1, y_b_1, x_r_1, y_t_1 = bbox_1  # Left, bottom, right, top
    x_l_2, y_b_2, x_r_2, y_t_2 = bbox_2
    bbox = [max(x_l_1, x_l_2), max(y_b_1, y_b_2), min(x_r_1, x_r_2), min(y_t_1, y_t_2)]
    x_l, y_b, x_r, y_t = bbox
    return True if x_r >= x_l and y_t >= y_b else False


def common_bbox(bins):  # Returns lower left corner, upper right corner in 1D array
    x_l_1, y_b_1, x_r_1, y_t_1 = fpd.bounding_box(bins[0])[1]  # Left, bottom, right, top
    x_l_2, y_b_2, x_r_2, y_t_2 = fpd.bounding_box(bins[1])[1]
    bbox = [max(x_l_1, x_l_2), max(y_b_1, y_b_2), min(x_r_1, x_r_2), min(y_t_1, y_t_2)]
    x_l, y_b, x_r, y_t = bbox
    return bbox if x_r >= x_l and y_t >= y_b else None


# --------- METHODS REQUIRING IMPLEMENTATIONS IN FPDE ---------------
def get_chunks_idxs_within_bounds(bin, bbox):
    chunks_bounds, _ = calculate_chunks_bounds(bin)
    chunks_idxs = [i for i in range(len(chunks_bounds)) if is_bbox_intersecting(bbox, chunks_bounds[i])]
    return chunks_idxs


def get_chunk(bin, idx, include_next=True):
    chunks = fpd.get_chunks(bin)
    vertices = chunks[idx] + [chunks[(idx + 1) % len(chunks)][0]]
    return vertices


# --------- /END/ METHODS REQUIRING IMPLEMENTATIONS IN FPDE ---------------


def is_intersecting(bins):
    bbox = common_bbox(bins)
    if bbox == None:
        return False

    # Bounding boxes intersect. Assume intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    chks = []
    segments = [[], []]
    for i in range(2):
        chks.append(get_chunks_idxs_within_bounds(bins[i], bbox))
        # Create list of segments for each chunk
        chunks_segments = [shapely.LineString(get_chunk(bins[i], c)) for c in chks[i]]
        segments[i] += chunks_segments
        for cs in chunks_segments:
            plot_geometry(cs)

    intersecting_points = []
    for i in segments[0]:
        for j in segments[1]:
            if i.intersects(j):
                pts = i.intersection(j)
                intersecting_points += list(shapely.get_coordinates(pts))

    plot_intersecting_points(intersecting_points)

    plot_chunks_bounds(bins[0], include_next_chunk_start=True, avoid_show=True, idxs=chks[0])
    plot_chunks_bounds(bins[1], include_next_chunk_start=True, avoid_create_frame=True, idxs=chks[1])

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box
    return True
