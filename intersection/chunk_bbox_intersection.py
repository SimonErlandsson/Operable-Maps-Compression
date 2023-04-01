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

    if x_r < x_l or y_t < y_b:
        type = 'None'
    elif x_l_1 == x_l_2 and y_b_1 == y_b_2 and x_r_1 == x_r_2 and y_t_1 == y_t_2:
        type = 'Equal'
    elif x_l_1 > x_l_2 and y_b_1 > y_b_2 and x_r_1 < x_r_2 and y_t_1 < y_t_2:
        type = '1 in 2'
    elif x_l_1 < x_l_2 and y_b_1 < y_b_2 and x_r_1 > x_r_2 and y_t_1 > y_t_2:
        type = '2 in 1'
    else:
        type = 'Partial'

    return (bbox, type) if type != 'None' else (None, type)


# --------- METHODS REQUIRING IMPLEMENTATIONS IN FPDE ---------------
def get_chunks_idxs_within_bounds(bin, bbox):
    chunks_bounds, _ = calculate_chunks_bounds(bin)
    chunks_idxs = [i for i in range(len(chunks_bounds)) if is_bbox_intersecting(bbox, chunks_bounds[i])]
    return chunks_idxs


def get_chunk(bin, idx, include_next=True):
    chunks, is_last_chunk_ring = fpd.get_chunks(bin)
    vertices = chunks[idx]
    chk_idx_ring_start = 0
    if is_last_chunk_ring[idx]:
        vertices += [chunks[chk_idx_ring_start][0]]
        chk_idx_ring_start = idx + 1
    else:
        vertices += [chunks[idx + 1][0]]
    return vertices


# --------- /END/ METHODS REQUIRING IMPLEMENTATIONS IN FPDE ---------------

PLOT_ONLY_INCORRECT = True


def is_intersecting(bins, debug_correct_ans):
    bbox, overlap_type = common_bbox(bins)
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

    intersecting_points = []
    for i in segments[0]:
        for j in segments[1]:
            if i.intersects(j):
                pts = i.intersection(j)
                intersecting_points += list(shapely.get_coordinates(pts))
                # return True # Optimize, disable for plotting

    if not PLOT_ONLY_INCORRECT or debug_correct_ans != None and debug_correct_ans != (len(intersecting_points) > 0):
        for cs in segments[0] + segments[1]:
            plot_geometry(cs)
        plot_intersecting_points(intersecting_points)
        plot_chunks_bounds(bins[0], include_next_chunk_start=True, avoid_show=True, idxs=chks[0])
        plot_chunks_bounds(bins[1], include_next_chunk_start=True, avoid_create_frame=True, idxs=chks[1])

    if len(intersecting_points) > 0:
        return True
    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box

    def is_contained_within(containee, container):
        origin = fpd.access_vertex(containee, 0)[0]
        # TODO: Create ray and check
        return True

    if overlap_type == '1 in 2':
        return is_contained_within(bins[0], bins[1])
    elif overlap_type == '2 in 1':
        return is_contained_within(bins[1], bins[0])
    return False
