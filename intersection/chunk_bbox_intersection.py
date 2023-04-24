from collections import defaultdict
import shapely
import numpy as np
from algos.alg_fpd_extended import FpdExtended
from intersection.plotting import plot_chunks_bounds, plot_geometry, plot_intersecting_points, create_canvas
import matplotlib.pyplot as plt
import math

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
# def slow_get_chunk(bin, idx, include_next=True):
#     chunks, is_last_chunk_ring = fpd.get_chunks(bin)
#     _, type = fpd.type(bin)
#     vertices = chunks[idx]
#     if not is_last_chunk_ring[idx] and type != 'LineString':
#         vertices += [chunks[idx + 1][0]]
#     return vertices

#chunks_bounds, _ = calculate_chunks_bounds(bin)
# --------- /END/ METHODS REQUIRING IMPLEMENTATIONS IN FPDE ---------------
def get_chunks_idxs_within_bounds(bin, bbox):
    chunks_bounds = fpd.get_chunk_bounds(bin)
    chunks_idxs = [i for i in range(len(chunks_bounds)) if is_bbox_intersecting(bbox, chunks_bounds[i])]
    return chunks_idxs

get_chunk = lambda bin, idx: fpd.get_chunk(bin, idx)[0] # Can also use slow above for debugging

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
    # if plot_all or debug_correct_ans != None and debug_correct_ans != (len(intersecting_points) % 2 == 1):
    #     print(fpd.type(container)[1], fpd.type(containee)[1])
    #     plot_chunks_bounds(container, include_next_chunk_start=True, avoid_show=True)
    #     plot_chunks_bounds(containee, include_next_chunk_start=True, avoid_create_frame=True)
    #     plot_chunks_bounds(container, include_next_chunk_start=True, avoid_show=True, idxs=chks)
    #     for cs in segments:
    #         plot_geometry(cs)
    #     plot_geometry(ray, solid=False)
    #     plot_intersecting_points(intersecting_points)
    #     plot_chunks_bounds(containee, include_next_chunk_start=False, avoid_create_frame=True, idxs=[], txt=f" : was {len(intersecting_points) % 2 == 1} expected {debug_correct_ans}")
    # END DEBUG
    return len(intersecting_points) % 2 == 1

def is_point_on_segment(seg, pt):
    seg_pt_1, seg_pt_2, pt = (np.array(seg[0]), np.array(seg[1]), np.array(pt))
    return abs(np.linalg.norm(seg_pt_1 - pt) + np.linalg.norm(seg_pt_2 - pt) - np.linalg.norm(seg_pt_1 - seg_pt_2)) < 1e-15

# Based on the common bbox, extracts the chunks for both geometries within the bbox,
# and performs intersection testing between the line segments.
def line_intersection(bins, bbox, debug_correct_ans, res_list=None, plot_all=False):
    chk_idxs = [[], []]
    chks = [[], []]
    polylines = [[], []]
    for i in range(2):
        chk_idxs[i] = get_chunks_idxs_within_bounds(bins[i], bbox)
        chks[i] = [get_chunk(bins[i], c_i) for c_i in chk_idxs[i]]
        # Each chunk becomes a polyline
        polylines[i] = [chunk_to_shape(c) for c in chks[i]]

    intersecting_points = []

    for i in polylines[0]:
        for j in polylines[1]:
            if i.intersects(j):
                # DEBUG ------------------
                # if plot_all or debug_correct_ans != None and debug_correct_ans != True:
                #     for cs in polylines[0] + polylines[1]:
                #         plot_geometry(cs)
                #     plot_intersecting_points(list(shapely.get_coordinates(i.intersection(j))))
                #     plot_chunks_bounds(bins[0], include_next_chunk_start=True, avoid_show=True, idxs=chk_idxs[0])
                #     plot_chunks_bounds(bins[1], include_next_chunk_start=True, avoid_create_frame=True, idxs=chk_idxs[1], txt=f" : was True expected {debug_correct_ans}")
                # END ----------------------
                if res_list == None:
                    return True
                intersecting_points += list(shapely.get_coordinates(i.intersection(j)))

    if len(intersecting_points) == 0:
        return False

    ## Append to res_list
    res_list.append(intersecting_points)
    segments = [[], []]
    seg_to_cross = [defaultdict(list), defaultdict(list)]
    cross_to_seg = [[[], []] for _ in range(len(intersecting_points))]
    for s in range(2):
        segments[s] = [[c[i], c[i+1]] for c in chks[s] for i in range(len(c) - 1)]


        # Fix check restart
        for seg_idx, seg in enumerate(segments[s]):
            for p_idx, p in enumerate(intersecting_points):
                #plot_intersecting_points([p])
                if is_point_on_segment(seg, p):
                    #plot_geometry(shapely.LineString(seg))
                    seg_to_cross[s][seg_idx].append(p_idx)
                    cross_to_seg[p_idx][s].append(seg_idx)
                # else:
                #     plot_geometry(shapely.LineString(seg), solid=False)
            seg_to_cross[s][seg_idx].sort(key=lambda x: x[0])
    res_list.append(segments)
    res_list.append(seg_to_cross)
    res_list.append(cross_to_seg)


def is_intersecting(bins, debug_correct_ans=None, plot_all=False):
    bbox, overlap_type = common_bbox(bins)
    if bbox == None:
        return False

    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    if line_intersection(bins, bbox, debug_correct_ans, plot_all=plot_all):
        return True

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box

    if overlap_type == '1 in 2':
        return is_contained_within(bins[0], bins[1], debug_correct_ans, plot_all)
    elif overlap_type == '2 in 1':
        return is_contained_within(bins[1], bins[0], debug_correct_ans, plot_all)
    return False

def intersection(bins, debug_correct_ans=None, plot_all=False):
    bbox, overlap_type = common_bbox(bins)
    if bbox == None:
        return shapely.Polygon(None)

    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    line_data = []
    line_intersection(bins, bbox, debug_correct_ans, line_data, plot_all)

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box
    if len(line_data) == 0:
        if overlap_type == '1 in 2' and is_contained_within(bins[0], bins[1], debug_correct_ans, plot_all):
            return fpd.decompress(bins[0])[1]
        elif overlap_type == '2 in 1' and is_contained_within(bins[1], bins[0], debug_correct_ans, plot_all):
            return fpd.decompress(bins[1])[1]
        return shapely.Polygon(None)

    # Have intersecting points, construct resulting polygon
    # 1. Create set of intersection points, mapping: intersection point <-> line segments
    # 2. Take random intersection point from set, follow path inside both shapes
    # 3. Continue until encountering intersection point or already visited point

    # Returns the next unvisited shapes inside both shapes connected to c_i
    def next_point(segments, seg_to_cross, cross_to_seg, c_i):
        seg_idxs = cross_to_seg[c_i][s]

    intersecting_points, segments, seg_to_cross, cross_to_seg = line_data
    cross_left = set(range(len(intersecting_points)))
    while len(cross_left) > 0:
        c_i = cross_left.pop()
        for segs in cross_to_seg[c_i]:
            next_point = segs

    create_canvas()
    plot_intersecting_points(line_data[0])
    plot_intersecting_points(line_data[1][0])
    plot_intersecting_points(line_data[1][1])
    plt.show()
    return line_data
