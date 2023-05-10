from collections import defaultdict
import shapely
import numpy as np
from algos.alg_fpd_extended import FpdExtended
from intersection.plotting import plot_chunks_bounds, plot_geometry, plot_intersecting_points, create_canvas, plot_line
from algos.fpd_extended_lib.cfg import *
import matplotlib.pyplot as plt
import math
from collections import deque 
import itertools
import bisect 



fpd = FpdExtended()
@profile
def is_bbox_intersecting(bbox_1, bbox_2):
    bbox = [max(bbox_1[0], bbox_2[0]), max(bbox_1[1], bbox_2[1]), min(bbox_1[2], bbox_2[2]), min(bbox_1[3], bbox_2[3])]
    return bbox[2] >= bbox[0] and bbox[3] >= bbox[1]

@profile
def common_bbox(bins):  # Returns lower left corner, upper right corner in 1D array
    x_l_1, y_b_1, x_r_1, y_t_1 = fpd.bounding_box(bins[0])[1]  # Left, bottom, right, top
    x_l_2, y_b_2, x_r_2, y_t_2 = fpd.bounding_box(bins[1])[1]
    bbox = [max(x_l_1, x_l_2), max(y_b_1, y_b_2), min(x_r_1, x_r_2), min(y_t_1, y_t_2)]
    x_l, y_b, x_r, y_t = max(x_l_1, x_l_2), max(y_b_1, y_b_2), min(x_r_1, x_r_2), min(y_t_1, y_t_2) # Common bounding box

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
@profile
def get_chunks_idxs_within_bounds(bin, bbox, get_geom_bounds=False):
    chunks_bounds = fpd.get_chunk_bounds(bin)
    chk_idxs = [i for i in range(len(chunks_bounds)) if is_bbox_intersecting(bbox, chunks_bounds[i])]
    return chk_idxs if not get_geom_bounds else chk_idxs, dict(zip(range(len(chunks_bounds)),chunks_bounds))


get_chunk = lambda bin, idx: fpd.get_chunk(bin, idx)[0] # Can also use slow above for debugging

@profile
def chunk_to_shape(chk): return shapely.Point(chk[0]) if len(chk) == 1 else shapely.LineString(chk)

@profile
def is_contained_within(containee, container, container_bounds, debug_correct_ans=None, plot_all=False):
    '''
    Containee is either FPDE binary object or a tuple of coordinates.
    '''
    if fpd.type(container)[1] == 'LineString':
        return False

    # Origin
    if type(containee) == bytes:
        x, y = fpd.access_vertex(containee, 0)[0]
    else:
        x, y = containee

    # Outher bounding box
    x_l, y_b, x_r, y_t = fpd.bounding_box(container)[1]
    distances = [(x_l - x, 0), (0, y_b - y), (x_r - x, 0), (0, y_t - y)]
    # Is point outside bounding box?
    if distances[0][0] > 0 or distances[1][1] > 0 or distances[2][0] < 0 or distances[3][1] < 0:
        return False
    ray_end = min(distances, key=lambda x: np.linalg.norm(x))
    ray_end = (x + ray_end[0], y + ray_end[1])
    ray = shapely.LineString([(x, y), ray_end])
    ray_bound = ray.bounds

    if OPTIMIZED_INTERSECTION:
        chks = [idx for idx, bound in container_bounds.items() if is_bbox_intersecting(ray_bound, bound)]
    else:
        chks = [i for i in range(len(container_bounds)) if is_bbox_intersecting(ray_bound, container_bounds[i])]

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
@profile
def is_point_on_segment(seg, pt):
    if  OPTIMIZED_INTERSECTION:
        x, y = pt
        x1, y1 = seg[0]
        x2, y2 = seg[1]

        cross = (y - y1) * (x2 - x1) - (x - x1) * (y2 - y1)

        # Check if the point is collinear with the segment
        if abs(cross) > POINT_ERR_TOL:
            return False

        dot_product = (x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)
        segment_length = (x2 - x1) ** 2 + (y2 - y1) ** 2

        if dot_product < -POINT_ERR_TOL or dot_product > segment_length + POINT_ERR_TOL:
            return False
        return True
    else:
        seg_pt_1, seg_pt_2, pt = (np.array(seg[0]), np.array(seg[1]), np.array(pt))
        return abs(np.linalg.norm(seg_pt_1 - pt) + np.linalg.norm(seg_pt_2 - pt) - np.linalg.norm(seg_pt_1 - seg_pt_2)) < 1e-15

   

# Based on the common bbox, extracts the chunks for both geometries within the bbox,
# and performs intersection testing between the line segments.
@profile
def line_intersection(bins, bbox, debug_correct_ans, res_list=None, plot_all=False):
    #General variables for whole geometries
    bounds =  [[], []]

    chks =  [[], []]
    chk_idxs =  [[], []]
    polylines = [[], []]
    segments =  [[], []]
    seg_idxs =  [[], []]
    #For not needing to recalculate boudning boxes and polylines in common bb
    chk_coords =    [[], []]
    chk_polylines = [[], []]
    chk_segments =  [[], []]

    #Variable for intersection
    intersecting_points = []


    for i in range(2):
        #Get all bounds for geometry as well as idxes which overlap in common bounding box
        chk_idxs[i], bounds[i] = get_chunks_idxs_within_bounds(bins[i], bbox, get_geom_bounds=True)
        chk_coords[i] = {c_i: get_chunk(bins[i], c_i) for c_i in chk_idxs[i]} # Get chunk to coords_data
        chk_polylines[i] = {c_i: chunk_to_shape(coords) for c_i, coords in chk_coords[i].items()} #Chunk -> polylines used in isintersection checks
        chk_segments[i] = {c_i: [(coords[i], coords[i+1]) for i in range(len(coords) - 1)] for c_i, coords in chk_coords[i].items()} # Get chunk -> segments data
        # Each chunk becomes a polyline
        chks[i] = list(chk_coords[i].values()) # Get chunk coords
        polylines[i] = list(chk_polylines[i].values()) # Transform each chunk to polyline
        segments[i] = list(itertools.chain(*chk_segments[i].values())) #Get segment data for geometry in total
        seg_idxs[i] = {segments[i][idx]: idx for idx in range(len(segments[i]))} # Get chunk -> segments data

    if OPTIMIZED_INTERSECTION:
        intersecting_points_idxs = dict() #To not have duplicates

        #Variables to save segment -> crossings and crossings -> segments
        seg_to_cross = [defaultdict(list), defaultdict(list)]
        cross_to_seg = []
        current_p_idx = 0
        #Loops through both geometrie's chunks, checks if polylines for those chunks intersect and if, save segment <-> crossings
        for chk_idx1, polylines1 in chk_polylines[0].items():
            for chk_idx2, polylines2 in chk_polylines[1].items():
                if polylines1.intersects(polylines2):
                    if res_list == None: # If doing predicate, return here. Else save point
                        return True, bounds
                    
                    curr_intersect_pnts = set((*coord,) for coord in list(shapely.get_coordinates(polylines1.intersection(polylines2)))) # Get all intersection points
                    # If intersection, save crossing <-> segment_idx
                    for p in curr_intersect_pnts: #For each intersection point found
                        #If point has not been extracted before get its index
                        if p not in intersecting_points_idxs:
                            cross_to_seg.append(([], []))
                            p_idx = current_p_idx
                            current_p_idx += 1
                            intersecting_points.append(p) #Point have not been discovered before so add to intersection points list (No duplicates allowed)
                            intersecting_points_idxs[p] = p_idx
                        else:
                            p_idx = intersecting_points_idxs[p]
                        for s, chk_idx in [(0,chk_idx1), (1, chk_idx2)]: #Go through both chunks for both geometries
                            chk_segs = chk_segments[s][chk_idx] #Extract the relevant segment
                            for seg in chk_segs: #Go through each of those segments
                                seg_idx = seg_idxs[s][seg] #Get the correct segment index !O(n)!
                                if not seg_idx in cross_to_seg[p_idx][s] and is_point_on_segment(seg, p):
                                #Checking if a segment intersect with intersection point or line
                                    seg_to_cross[s][seg_idx].append(p_idx)
                                    cross_to_seg[p_idx][s].append(seg_idx)
            for s in range(2):
                for seg_idx in range(len(seg_to_cross[s])):    
                    seg_to_cross[s][seg_idx].sort(key=lambda x: math.dist(segments[s][seg_idx][0],intersecting_points[x])) # Sort ordered by distance from p[0]
        if len(intersecting_points) == 0:
            return False, bounds 
    else:
        intersecting_points_set = set() #To not have duplicates
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
                    if res_list == None: # If doing predicate, return here. Else save point
                        return True, bounds
                    intersect_points = list(shapely.get_coordinates(i.intersection(j)))
                    for point in intersect_points: 
                        if str(point) not in intersecting_points_set:
                            intersecting_points.append(point)
                            intersecting_points_set.add(str(point))
            
        if len(intersecting_points) == 0:
            return False, bounds
        ## Append to res_list
        seg_to_cross = [defaultdict(list), defaultdict(list)]
        cross_to_seg = [[[], []] for _ in range(len(intersecting_points))]
        for s in range(2):
            # Find which segments contain intersecting points
            for seg_idx, seg in enumerate(segments[s]):
                for p_idx, p in enumerate(intersecting_points):
                    #plot_intersecting_points([p])
                    if is_point_on_segment(seg, p):
                        #plot_geometry(shapely.LineString(seg))
                        seg_to_cross[s][seg_idx].append(p_idx)
                        cross_to_seg[p_idx][s].append(seg_idx)
                    # else:
                    #     plot_geometry(shapely.LineString(seg), solid=False)
                seg_to_cross[s][seg_idx].sort(key=lambda x: np.linalg.norm(segments[s][seg_idx][0] - intersecting_points[x])) # Sort ordered by distance from p[0]
    res_list.append(intersecting_points)
    res_list.append(segments)
    res_list.append(seg_to_cross)
    res_list.append(cross_to_seg)
    res_list.append(bounds)



@profile
def is_intersecting(bins, debug_correct_ans=None, plot_all=False):
    bbox, overlap_type = common_bbox(bins)
    if bbox == None:
        return False

    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    intersects, bounds = line_intersection(bins, bbox, debug_correct_ans, plot_all=plot_all)
    if intersects:
        return True

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box

    if overlap_type == '1 in 2':
        return is_contained_within(bins[0], bins[1], bounds[1], debug_correct_ans=debug_correct_ans, plot_all=plot_all)
    elif overlap_type == '2 in 1':
        return is_contained_within(bins[1], bins[0], bounds[0], debug_correct_ans=debug_correct_ans, plot_all=plot_all)
    return False

@profile
def intersection(bins, debug_correct_ans=None, plot_all=False):
    bbox, overlap_type = common_bbox(bins)
    if bbox == None:
        return shapely.Polygon(None)

    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    line_data = []
    result = line_intersection(bins, bbox, debug_correct_ans, line_data, plot_all)

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box
    if len(line_data) == 0:
        _, bounds = result
        if overlap_type == '1 in 2' and is_contained_within(bins[0], bins[1], bounds[1], debug_correct_ans=debug_correct_ans, plot_all=plot_all):
            return fpd.decompress(bins[0])[1] # Return whole smaller shape
        elif overlap_type == '2 in 1' and is_contained_within(bins[1], bins[0], bounds[0], debug_correct_ans=debug_correct_ans, plot_all=plot_all):
            return fpd.decompress(bins[1])[1]
        return shapely.Polygon(None)

    # Have intersecting points, construct resulting polygon
    # 1. Create set of intersection points, mapping: intersection point <-> line segments
    # 2. Take random intersection point from set, follow path inside both shapes
    # 3. Continue until encountering intersection point or already visited segment

    # Takes a segment and vertex_index (can be > 2 if cross-point, idx 0 is beg, idx 1 is end)
    def seg_to_point(s, seg_idx, v_idx):
        if v_idx > 1:
            return intersecting_points[seg_to_cross[s][seg_idx][v_idx - 2]]
        else:
            return segments[s][seg_idx][v_idx]

    # Find the point in the middle of the segment, taking into account cross points
    def seg_to_middle_point(s, seg_idx, v_idxs):
        l = seg_to_point(s, seg_idx, v_idxs[0])
        r = seg_to_point(s, seg_idx, v_idxs[1])
        return [(l[0] + r[0]) / 2, (l[1] + r[1]) / 2]


    def DEBUG_print_paths(paths, c_i=None):
        if c_i != None:
            print(f"--: {intersecting_points[c_i]} :--")
        #print("Raw:", paths)
        str = ""
        for p in paths:
            s, seg_idx, v_idxs, p_dir = p
            p_dir_sym = '+' if p_dir == 1 else '-'
            str += f'[{seg_to_point(s, seg_idx, v_idxs[0])} -> {seg_to_point(s, seg_idx, v_idxs[1])} ({p_dir_sym})], '
        print(str)


    # Returns the possible paths (directed segment) from an intersection point. Also checks that it is within both shapes.
    @profile
    def possible_paths(c_i, bounds):
        possible_paths = []
        seg_idxs = cross_to_seg[c_i]
        for s in range(2): # Both shapes
            for seg_idx in seg_idxs[s]: # Enumerate the segments containing crosspoint with id c_i
                seg_cross_cnt = len(seg_to_cross[s][seg_idx])
                # Get possible successor points
                # Where is the current cross?
                c_i_in_seg = seg_to_cross[s][seg_idx].index(c_i) + 2 # Index of vertex within segment

                # Get previous cross point if exists
                start_v = c_i_in_seg - 1 if c_i_in_seg != 2 else 0 # V:0 left of segment, V:1 right in segment, V:2+ intersection points
                if not np.array_equal(seg_to_point(s, seg_idx, start_v), seg_to_point(s, seg_idx, c_i_in_seg)): # Dont add line segments consisting of one point
                    possible_paths.append((s, seg_idx, (start_v, c_i_in_seg), -1)) # Vertex 2 being the first cross, 0 is first in shape order

                # Has cross point after?
                end_v = c_i_in_seg + 1 if c_i_in_seg - 1 != seg_cross_cnt else 1
                if not np.array_equal(seg_to_point(s, seg_idx, c_i_in_seg), seg_to_point(s, seg_idx, end_v)):
                    possible_paths.append((s, seg_idx, (c_i_in_seg, end_v), 1)) # Vertex 2 being the first cross, 0 is first in shape order

        paths_to_save = set()
        for i in range(len(possible_paths)):
            for j in range(len(possible_paths)):
                if j != i:
                    s1, seg_idx1, v_idxs1, _ = possible_paths[i]
                    s2, seg_idx2, v_idxs2, _ = possible_paths[j]
                    path1 = shapely.LineString([seg_to_point(s1, seg_idx1, v_idxs1[0]), seg_to_point(s1, seg_idx1, v_idxs1[1])])
                    path2 = shapely.LineString([seg_to_point(s2, seg_idx2, v_idxs2[0]), seg_to_point(s2, seg_idx2, v_idxs2[1])])
                    if path1.intersection(path2).geom_type == "LineString":
                        paths_to_save.add(possible_paths[i])
                        paths_to_save.add(possible_paths[j])
                        
        possible_paths = list(filter(lambda p: is_contained_within(seg_to_middle_point(*p[0:3]), bins[(p[0] + 1) % 2], bounds[(p[0] + 1) % 2]), possible_paths))
        paths_to_save.update(possible_paths)
        
        possible_paths = list(paths_to_save)
        
        # print("")
        # DEBUG_print_paths(list(paths_to_save))
        # DEBUG_print_paths(possible_paths, c_i)
        # #Make sure points are within both shapes
        # print("")
        # DEBUG_print_paths(possible_paths)
        return possible_paths


    intersecting_points, segments, seg_to_cross, cross_to_seg, bounds = line_data
    #print("Intersecting Points:", intersecting_points)
    cross_left = set(range(len(intersecting_points))) # Ids of unprocessed intersection points
    processed_ways = [[set(), set()] for _ in range(len(intersecting_points))] # Avoid processing visited segments
    res_segs = deque() # Segments which are part of the resulting shape
    res_points = []
    visited_edges = set()
    while len(cross_left) > 0:
        c_i = cross_left.pop() # Take one cross-point
        paths = possible_paths(c_i, bounds) # Find possible paths from the cross-point
        if len(paths) == 0:
            res_points.append(intersecting_points[c_i])
        while len(paths) > 0:
            path = paths.pop() # Process one path

            s, seg_idx, v_idxs, p_dir = path # V_idxs contains: start_idx, end_index. Can also be cross-points!
            start_idx = 0 if p_dir == 1 else 1 # Start index. Avoids creating segments which are flipped
            end_idx = 1 if p_dir == 1 else 0

            if seg_idx in processed_ways[c_i][s]: # Skip processed
                continue
            
            path_segs = deque()
            path_append = path_segs.append if p_dir == 1 else lambda x: path_segs.insert(0, x)
            while True: # While no cross point

                v1, v2 = seg_to_point(s, seg_idx, v_idxs[start_idx]), seg_to_point(s, seg_idx, v_idxs[end_idx])
                
                #Add front and back edges to visited edges set
                if OPTIMIZED_INTERSECTION and (v1, v2) not in visited_edges and (v2, v1) not in visited_edges:
                    path_append([v1, v2]) # Add segment to resulting shape a
                    visited_edges.add((v1, v2))
                    visited_edges.add((v2, v1))

                elif not OPTIMIZED_INTERSECTION and str([v1, v2]) not in visited_edges and str([v2, v1]) not in visited_edges:
                    path_append([seg_to_point(s, seg_idx, v_idxs[0]), seg_to_point(s, seg_idx, v_idxs[1])]) # Add segment to resulting shape
                    visited_edges.add(str([v1, v2]))
                    visited_edges.add(str([v2, v1]))

                e_v = v_idxs[end_idx] # Find index of actual end vertex (i.e. flip segment based on direction)

                next_seg_idx = seg_idx + p_dir
                if e_v > 1: # Is cross point?
                    encountered_c_idx = seg_to_cross[s][seg_idx][e_v - 2] # Find cross_idx of collided cross-point
                    processed_ways[encountered_c_idx][s].add(seg_idx) # Add to shape's list for the cross-point
                    break
                elif next_seg_idx == -1 or next_seg_idx == len(segments[s]) or not np.array_equal(seg_to_point(s, seg_idx, e_v), seg_to_point(s, next_seg_idx, start_idx)): #<- HERE ERROR
                    break # Break if no more segments, or if path formed by segments is not continuous
                else:
                    seg_idx = next_seg_idx # Next segment

                    seg_cross_cnt = len(seg_to_cross[s][seg_idx])
                    if p_dir == 1:
                        v_idxs = [0, 2 if seg_cross_cnt != 0 else 1] # Has next line-segment crosspoint? If so take the first cross point as segment end-point.
                    else:
                        v_idxs = [2 if seg_cross_cnt != 0 else 0, 1]

            # Append path-segment to total segments
            # TODO: Fix ordering in resulting shape
            if p_dir == 1:
                res_segs += path_segs
            else:
                path_segs += res_segs
                res_segs = path_segs 


    #Merge all segments into LineString or MultiLineString
    unfilt_line_strs = shapely.ops.linemerge(res_segs)
    type_unfilt = unfilt_line_strs.geom_type
    unfilt_line_strs = [unfilt_line_strs] if type_unfilt == "LineString" else list(unfilt_line_strs.geoms) #For making MultiLineString and LineString be handleded similarly

    #List of all resulting linestrings and polygons
    line_strs, polygons = [], []

    #Check if LineString has the shape of a polygon, then convert it and divide polygons and line strings
    for line_str in unfilt_line_strs:
        if line_str.is_ring:
            polygons.append(shapely.Polygon(line_str.coords))
        else:
            line_strs.append(line_str)

    #Variables for different cases of geometries
    is_MultiLineString, is_MultiPolygon, is_MultiPoint = len(line_strs) > 1, len(polygons) > 1, len(res_points) > 1
    has_LineString, has_Polygon, has_Point = len(line_strs) > 0, len(polygons) > 0, len(res_points) > 0

    result = []
    if is_MultiLineString:
        result.append(shapely.MultiLineString(line_strs))
    if is_MultiPolygon:
        result.append(shapely.MultiPolygon(polygons))
    if is_MultiPoint:
        result.append(shapely.MultiPoint(res_points))
    if not is_MultiLineString and has_LineString:
        result.append(line_strs[0])
    if not is_MultiPolygon and has_Polygon:
        result.append(polygons[0])
    if not is_MultiPoint and has_Point:
        result.append(shapely.Point(res_points[0]))

    if len(result) > 1:
        return shapely.GeometryCollection(result)
    
    elif len(result) == 1:
        return result[0]

    return res_segs
