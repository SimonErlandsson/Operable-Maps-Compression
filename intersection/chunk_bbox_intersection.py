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
    #chunks_bounds, _ = calculate_chunks_bounds(bin)
    chunks_bounds = fpd.get_chunk_bounds(bin)
    chunks_idxs = [i for i in range(len(chunks_bounds)) if is_bbox_intersecting(bbox, chunks_bounds[i])]
    return chunks_idxs


def get_chunk(bin, idx, include_next=True):
    chunks, is_last_chunk_ring = fpd.get_chunks(bin)
    _, type = fpd.type(bin)
    vertices = chunks[idx]
    if not is_last_chunk_ring[idx] and type != 'LineString':
        vertices += [chunks[idx + 1][0]]
    return vertices


# --------- /END/ METHODS REQUIRING IMPLEMENTATIONS IN FPDE ---------------
def chunk_to_shape(chk): return shapely.Point(chk[0]) if len(chk) == 1 else shapely.LineString(chk)


def is_contained_within(containee, container, debug_correct_ans, plot_all=False):
    if fpd.type(container)[1] == 'LineString':
        return False

    # Origin
    x, y = fpd.access_vertex(containee, 0)[0]
    # Outher bounding box
    x_l, y_b, x_r, y_t = fpd.bounding_box(container)[1]
    distances = [(x_l - x, 0), (0, y_b - y), (x_r - x, 0), (0, y_t - y)]
    ray_end = min(distances, key=lambda x: np.linalg.norm(x))
    ray_end = (x + ray_end[0], y + ray_end[1])
    ray = shapely.LineString([(x, y), ray_end])

    # Create list of segments for other shape, if the chunk collides with LineString
    chks = get_chunks_idxs_within_bounds(container, ray.bounds)
    segments = [chunk_to_shape(get_chunk(container, c)) for c in chks]
    intersecting_points = []
    for i in segments:
        if i.intersects(ray):
            intersecting_points += list(shapely.get_coordinates(i.intersection(ray)))
    # DEBUG
    if plot_all or debug_correct_ans != None and debug_correct_ans != (len(intersecting_points) % 2 == 1):
        print(fpd.type(container), fpd.type(containee))
        plot_chunks_bounds(container, include_next_chunk_start=True, avoid_show=True)
        plot_chunks_bounds(containee, include_next_chunk_start=True, avoid_create_frame=True)
        plot_chunks_bounds(container, include_next_chunk_start=True, avoid_show=True, idxs=chks)
        for cs in segments:
            plot_geometry(cs)
        plot_geometry(ray, solid=False)
        plot_intersecting_points(intersecting_points)
        plot_chunks_bounds(containee, include_next_chunk_start=False, avoid_create_frame=True, idxs=[], txt=f" : was {len(intersecting_points) % 2 == 1} expected {debug_correct_ans}")
    # END DEBUG
    return len(intersecting_points) % 2 == 1

# Based on the common bbox, extracts the chunks for both geometries within the bbox,
# and performs intersection testing between the line segments.


def has_line_intersection(bins, bbox, debug_correct_ans, plot_all=False):
    chks = [[], []]
    segments = [[], []]
    for i in range(2):
        chks[i] = get_chunks_idxs_within_bounds(bins[i], bbox)
        # Create list of segments for each chunk
        segments[i] = [chunk_to_shape(get_chunk(bins[i], c)) for c in chks[i]]

    for i in segments[0]:
        for j in segments[1]:
            if i.intersects(j):
                intersecting_points = list(shapely.get_coordinates(i.intersection(j)))
                # DEBUG ------------------
                if plot_all or debug_correct_ans != None and debug_correct_ans != (len(intersecting_points) > 0):
                    for cs in segments[0] + segments[1]:
                        plot_geometry(cs)
                    plot_intersecting_points(intersecting_points)
                    plot_chunks_bounds(bins[0], include_next_chunk_start=True, avoid_show=True, idxs=chks[0])
                    plot_chunks_bounds(bins[1], include_next_chunk_start=True, avoid_create_frame=True, idxs=chks[1], txt=f" : was {len(intersecting_points) > 0} expected {debug_correct_ans}")
                # END ----------------------
                return True
    return False


def is_intersecting(bins, debug_correct_ans, plot_all=False):
    bbox, overlap_type = common_bbox(bins)
    if bbox == None:
        return False

    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    if has_line_intersection(bins, bbox, debug_correct_ans, plot_all):
        return True

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box

    if overlap_type == '1 in 2':
        return is_contained_within(bins[0], bins[1], debug_correct_ans, plot_all)
    elif overlap_type == '2 in 1':
        return is_contained_within(bins[1], bins[0], debug_correct_ans, plot_all)
    return False
