from collections import defaultdict
from shapely import LineString, Point, MultiPolygon, MultiLineString, MultiPoint, GeometryCollection, Polygon, ops, get_coordinates
import numpy as np
from algos.alg_fpd_extended import FpdExtended
from intersection.plotting import plot_chunks_bounds, plot_geometry, plot_intersecting_points, create_canvas, plot_line
from algos.fpd_extended_lib.cfg import *
import matplotlib.pyplot as plt
import math
from collections import deque 
import itertools
import time 
import algos.fpd_extended_lib.cfg as cfg

fpd = FpdExtended()

#Extract the bounding box based on settings
def bounding_box(bin, s):
    if not cfg.DISABLE_OPTIMIZED_INTERSECTION or 'glob_bounding_boxes' not in globals():
        return fpd.bounding_box(bin)[1]
    else:
        return glob_bounding_boxes[s]

#Calculates if two given bounding boxes overlap
def is_bboxs_intersecting(bbox_1, bbox_2):
    return min(bbox_1[2], bbox_2[2]) >= max(bbox_1[0], bbox_2[0]) and min(bbox_1[3], bbox_2[3]) >= max(bbox_1[1], bbox_2[1])


def calculate_final_stats(start_time):
        #Shifts the stats
        perf_counter[1] = perf_counter[1] + perf_counter[2]
        perf_counter[2] = perf_counter[3] + perf_counter[4]
        perf_counter[3] = time.perf_counter() - start_time        
        return perf_counter[:-1]

def select_return(get_stats, ret_value, start_time):
    if not get_stats:
        return ret_value
    else:
        return calculate_final_stats(start_time), ret_value
    

#Checks if segments are parallell to each other
def are_lines_parallel(seg1, seg2):
        (x1, y1), (x2, y2) = seg1[0], seg1[1]
        (x3, y3),(x4, y4) = seg2[0], seg2[1]

        return ((y2 - y1) / (x2 - x1) if (x2 - x1) != 0 else float('inf')) == ((y4 - y3) / (x4 - x3) if (x4 - x3) != 0 else float('inf'))
     

def get_chunk(bin, idx, offset_cache=None, glob_idx = None): 
    s = time.perf_counter()
    if not cfg.DISABLE_OPTIMIZED_INTERSECTION:
        perf_counter[1 + glob_idx] += 1
        res = fpd.get_chunk(bin, idx, offset_cache=offset_cache)[0] # Can also use slow above for debugging
    else: 
        res = [tuple(l) for l in glob_chunks[glob_idx][idx]]
    perf_counter[0] += time.perf_counter() - s
    return res


def chunk_to_shape(chk): return LineString(chk) if len(chk) != 1 else  Point(chk[0]) #Convert chunk segments to shapely object


#Calculate the global variables for chunks and bounding boxes when cfg.DISABLE_OPTIMIZED_INTERSECTION = True
def calculate_globals(bins):
    global perf_counter
    perf_counter = [0, 0, 0, 0, 0]
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        start_time = time.perf_counter()

        global glob_bounding_boxes
        global glob_chunk_bounds
        global glob_chunks

        glob_bounding_boxes = []
        glob_chunk_bounds = []
        glob_chunks = []

        for s in range(2):
            bin = bins[s]
            chunks = fpd.get_chunks(bin)[0]
            chunks_bounds = []
            #Global bbox
            g_x_l, g_y_b, g_x_r, g_y_t = 200, 200, -200, -200
            for chunk in chunks:
                #Chunk local bbox
                x_l, y_b, x_r, y_t = 200, 200, -200, -200
                for x, y in chunk:
                    x_l, y_b, x_r, y_t = min(x, x_l), min(y, y_b), max(x, x_r), max(y, y_t)
                    g_x_l, g_y_b, g_x_r, g_y_t = min(x, g_x_l), min(y, g_y_b), max(x, g_x_r), max(y, g_y_t)

                chunks_bounds.append([x_l, y_b, x_r, y_t])
            glob_chunk_bounds.append(chunks_bounds)
            glob_bounding_boxes.append([g_x_l, g_y_b, g_x_r, g_y_t])
            glob_chunks.append(chunks)
        perf_counter[0] += time.perf_counter() - start_time


#Calculate the common bounding box from two FPDE compressed binaries
def common_bbox(bins):  
    # Left, bottom, right, top of bounding boxes
    x_l_1, y_b_1, x_r_1, y_t_1 = bounding_box(bins[0], 0)
    x_l_2, y_b_2, x_r_2, y_t_2 = bounding_box(bins[1], 1)

    #Calculate common bounding box coordinates
    x_l, y_b, x_r, y_t = [max(x_l_1, x_l_2), max(y_b_1, y_b_2), min(x_r_1, x_r_2), min(y_t_1, y_t_2)]

    #Match common bounding box with context
    if x_r < x_l or y_t < y_b:
        return (None, 'None') #Return if common bounding box do is empty
    
    elif x_l_1 == x_l_2 and y_b_1 == y_b_2 and x_r_1 == x_r_2 and y_t_1 == y_t_2:
        type = 'Equal'
    elif x_l_1 > x_l_2 and y_b_1 > y_b_2 and x_r_1 < x_r_2 and y_t_1 < y_t_2:
        type = '1 in 2'
    elif x_l_1 < x_l_2 and y_b_1 < y_b_2 and x_r_1 > x_r_2 and y_t_1 > y_t_2:
        type = '2 in 1'
    else:
        type = 'Partial'

    return ([x_l, y_b, x_r, y_t], type)


#Get chunk indexes in FPDE binary containing a segment connected to the given bounding box.
# Also return the bounds of the geometry if get_geom_bounds=True
def get_chunks_idxs_within_bounds(bin, bbox, glob_idx, get_geom_bounds=False):
    s = time.perf_counter()
    if cfg.DISABLE_OPTIMIZED_INTERSECTION:
        chunks_bounds = glob_chunk_bounds[glob_idx]
        perf_counter[1 + glob_idx] += len(chunks_bounds)
    else:
        chunks_bounds = fpd.get_chunk_bounds(bin) #Get bounds from geometry

    chk_idxs = {i for i in range(len(chunks_bounds)) if is_bboxs_intersecting(bbox, chunks_bounds[i])} #Add chunk index to list if chunk overlap with bounding box
    res = chk_idxs if not get_geom_bounds else (chk_idxs, chunks_bounds)    
    
    perf_counter[0] += time.perf_counter() - s
    perf_counter[3 + glob_idx] += len(chunks_bounds)
   
    return res


def is_contained_within(containee, container, container_bounds, glob_idx = None, debug_correct_ans=None, plot_all=False, cache=None):
    '''
    Containee is either FPDE binary object or a tuple of coordinates. Uses Ray Casting algorithm to
    see if a point (containee) is inside container.
    PS. To solve for edge cases structure of algorithm makes method not to be called if containee is exactly on a container segment.
    '''

    #Return false since a a line string can not encapsulate another geometry
    if fpd.type(container)[1] == 'LineString':
        return False

    # Origin
    if type(containee) == bytes:
        x, y = fpd.access_vertex(containee, 0)[0]
    else:
        x, y = containee

    # Outher bounding box
    x_l, y_b, x_r, y_t = bounding_box(container, glob_idx)
    distances = [(x_l - x, 0), (0, y_b - y), (x_r - x, 0), (0, y_t - y)]
    # Is point outside bounding box?
    if distances[0][0] > 0 or distances[1][1] > 0 or distances[2][0] < 0 or distances[3][1] < 0:
        return False
    
    #Create ray
    ray_end = min(distances, key=lambda x: x[0]**2 + x[1]**2)
    ray_end = (x + ray_end[0], y + ray_end[1])
    ray = LineString([(x, y), ray_end])
    ray_bound = ray.bounds

    #Uses caching to not call get_chunk multiple times
    s = time.perf_counter()
    chks = [idx for idx, bound in enumerate(container_bounds) if is_bboxs_intersecting(ray_bound, bound)]
    perf_counter[0] += time.perf_counter() - s
    segments = []
    for c in chks:
        if c not in cache:
            #if not cfg.DISABLE_OPTIMIZED_INTERSECTION:
            chunk_shape = chunk_to_shape(get_chunk(container, c, glob_idx=glob_idx))
            cache[c] = chunk_shape
            segments.append(chunk_shape)
        else:
            segments.append(cache[c])

    #Count how many times the ray intersects container boundries
    intersecting_points = []
    for i in segments:
        if i.intersects(ray):
            intersecting_points += list(get_coordinates(i.intersection(ray)))

    # DEBUG
    # if plot_all or debug_correct_ans != None and debug_correct_ans != (len(intersecting_points) % 2 == 1):
    #     print(fpd.type(container)[1], fpd.type(containee)[1])
    #     create_canvas(no_frame=True, zoom=1.5)
    #     chunks_not_loaded = list(filter(lambda x: x not in chks, list(range(len(container_bounds.items())))))
    #     plot_chunks_bounds(container, include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, idxs=chunks_not_loaded, solid=False, alpha=0.3, fill_alpha=0.05, color='black', point_color='cyan', point_size=9)
    #     plot_chunks_bounds(containee, include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, alpha=1, color='orange', point_color='lime', point_size=14)
    #     plot_chunks_bounds(container, include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, idxs=chks, alpha=1, color='orange', point_color='red', point_size=12)
    #     ## for cs in segments:
    #     ##     plot_geometry(cs)
    #     plot_geometry(fpd.decompress(containee)[1], fill_alpha=0.1, hatch='\\', color='blue')
    #     plot_geometry(fpd.decompress(container)[1], fill_alpha=0.1, hatch='/', color='green')
    #     plot_geometry(fpd.decompress(containee)[1], fill_alpha=0.3, hatch='X', color='purple')
    #     plot_geometry(ray, solid=False)
    #     plot_intersecting_points(intersecting_points, color='orange', size=25)
    #     plt.show()
        ##plot_chunks_bounds(containee, include_next_chunk_start=False, avoid_create_frame=True, idxs=[], txt=f" : was {len(intersecting_points) % 2 == 1} expected {debug_correct_ans}")
    # END DEBUG

    return len(intersecting_points) % 2 == 1


# Based on the common bbox, extracts the chunks for both geometries within the bbox,
# and performs intersection testing between the line segments.
def line_intersection(bins, bbox, debug_correct_ans, res_list=None, plot_all=False, cache=None):
    offset_cache = [{}, {}]
    intersecting_points = []

    #General variables for whole geometries
    bounds =    [None, None]
    chks =      [None, None]
    chk_idxs =  [None, None]
    polylines = [None, None]
    segments =  [None, None]
    seg_idxs =  [None, None]

    #For not needing to recalculate boudning boxes and polylines in common bb
    chk_coords =    [None, None]
    chk_polylines = [None, None]
    chk_segments =  [None, None]
    

    for i in range(2):
        chk_idxs[i], bounds[i] = get_chunks_idxs_within_bounds(bins[i], bbox, glob_idx=i, get_geom_bounds=True) #Get all bounds for a geometry as well as idxes which overlap in common bounding box

    #Only include chunks which have any intersection with any other bbox chunk of the other geometry
    chks_filt =  [set(), set(), set()] #Last index is for saving chunks that overlap to minimize line intersection calls
    for chk1 in chk_idxs[0]:
        for chk2 in chk_idxs[1]:
            if is_bboxs_intersecting(bounds[0][chk1], bounds[1][chk2]):
                chks_filt[2].add((chk1, chk2))
                if not cfg.DISABLE_OPTIMIZED_INTERSECTION or res_list == None:
                    chks_filt[0].add(chk1)
                    chks_filt[1].add(chk2)
    
    #If Is_intersection: Do not create mock segments for "not yet" unfolded chunk
    if res_list == None:
        chk_idxs = chks_filt    
    
    #Extract data for the relevant chunks
    for i in range(2):
        if not DISABLE_OPTIMIZED_UNPACKING:
            #Create mock segment if chunk from geometry 1 do not intersecting with another chunk from geometry 2. For cfg.DISABLE_OPTIMIZED_INTERSECTION, it is unnessesary since  chunks is  already unpacked fully
            chk_coords[i] = {c_i: [(None,c_i), (None,c_i)] if c_i not in chks_filt[i] and not cfg.DISABLE_OPTIMIZED_INTERSECTION else get_chunk(bins[i], c_i, offset_cache=offset_cache[i], glob_idx=i) for c_i in chk_idxs[i]} # Get chunk -> coordinates for those inside common boudning box
            chk_polylines[i] = {c_i: chunk_to_shape(coords) for c_i, coords in chk_coords[i].items() if c_i in chks_filt[i] or cfg.DISABLE_OPTIMIZED_INTERSECTION} #Chunk -> polylines used in isintersection checks
        else:    
            chk_coords[i] = {c_i: get_chunk(bins[i], c_i, offset_cache=offset_cache[i], glob_idx=i) for c_i in chk_idxs[i]} # Get chunk -> coordinates for those inside common boudning box
            chk_polylines[i] = {c_i: chunk_to_shape(coords) for c_i, coords in chk_coords[i].items()} #Chunk -> polylines used in isintersection checks
        
        
        chks[i] = list(chk_coords[i].values()) # Get chunk coords
        polylines[i] = list(chk_polylines[i].values()) # Transform each chunk to polyline
        cache[i] = chk_polylines[i]

        
        #Avoid declaring unused variables if predicate intersection
        if res_list != None:
            chk_segments[i] = {c_i: [(coords[i], coords[i+1]) for i in range(len(coords) - 1)] for c_i, coords in chk_coords[i].items()} # Get chunk -> segments data
            segments[i] = list(itertools.chain(*chk_segments[i].values())) #Get segment data for geometry in total
            seg_idxs[i] = {segments[i][idx]: idx for idx in range(len(segments[i]))} # Get chunk -> segments data

    intersecting_points_idxs = dict() #To not have duplicates

    #Variables to save segment -> crossings and crossings -> segments
    seg_to_cross = [defaultdict(list), defaultdict(list)]
    cross_to_seg = []
    current_p_idx = 0
    #Loops through both geometrie's chunks, checks if polylines for those chunks intersect and if, save segment <-> crossings
    # for chk_idx1, polylines1 in chk_polylines[0].items():
    #     for chk_idx2, polylines2 in chk_polylines[1].items():
    for (chk_idx1, chk_idx2) in chks_filt[2]:
            polylines1, polylines2 = chk_polylines[0][chk_idx1], chk_polylines[1][chk_idx2]
            if polylines1.intersects(polylines2):
                
                if res_list == None: # If doing predicate, return here. Else save point
                    return True, bounds
                
                curr_intersect_pnts = set((*coord,) for coord in list(get_coordinates(polylines1.intersection(polylines2)))) # Get all intersection points
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
                        found = False #Each cross point can only have 2 segments from one geometry (non self intersecting)
                        chk_segs = chk_segments[s][chk_idx] #Extract the relevant segment
                        for seg in chk_segs: #Go through each of those segments   
                            seg_idx = seg_idxs[s][seg] #Get the correct segment index !O(n)!
                            if not seg_idx in cross_to_seg[p_idx][s]:
                                
                                #Inline for  checking if a point "p" is on a segment "seg"
                                if abs(math.hypot(seg[1][0] - seg[0][0], seg[1][1] - seg[0][1]) - (math.hypot(p[0] - seg[0][0], p[1] - seg[0][1]) +  math.hypot(seg[1][0] - p[0], seg[1][1] - p[1]))) < 1e-13:                               
                                    seg_to_cross[s][seg_idx].append(p_idx)
                                    cross_to_seg[p_idx][s].append(seg_idx)
                                    if found:
                                        break
                                    found = True 
       
    if res_list != None:
            for s in range(2):
                for seg_idx in seg_to_cross[s]:    
                    seg_to_cross[s][seg_idx].sort(key=lambda x: math.dist(segments[s][seg_idx][0],intersecting_points[x])) # Sort ordered by distance from p[0]

    # create_canvas(no_frame=True, zoom=1.5)
    # for s in range(2):
    #     chunks_not_loaded = list(filter(lambda x: x not in chk_idxs[s], list(range(len(bounds[s])))))
    #     plot_chunks_bounds(bins[s], include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, idxs=chunks_not_loaded, solid=False, alpha=0.3, fill_alpha=0.05, color='black', point_color='cyan', point_size=9)
    #     plot_chunks_bounds(bins[s], include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, idxs=chk_idxs[s], alpha=1.0, color='orange', point_color='red', point_size=14)
    #     plot_geometry(fpd.decompress(bins[s])[1], fill_alpha=0.1, hatch=('\\' if s == 0 else '/'), color=('blue' if s == 0 else 'green'))
    # from shapely import intersection as s_inter
    # inter_shape = s_inter(fpd.decompress(bins[0])[1], fpd.decompress(bins[1])[1])
    # plot_geometry(inter_shape, fill_alpha=0.3, hatch='X', color='purple')
    # plot_intersecting_points(get_coordinates(inter_shape), color='lime', zorder=100, size=14)
    # plt.show()


    if len(intersecting_points) == 0:
        return False, bounds 

    #Update results list (return value)
    res_list.append(intersecting_points)
    res_list.append(segments)
    res_list.append(seg_to_cross)
    res_list.append(cross_to_seg)
    res_list.append(bounds)
    res_list.append(offset_cache)


# Returns the possible paths (directed segment) from an intersection point. Also checks that it is within both shapes.
def possible_paths(c_i, bounds, cross_to_seg, seg_to_cross, seg_to_point, seg_to_middle_point, bins, cache):
    possible_paths, removed = [], False
    seg_idxs = cross_to_seg[c_i]
    for s in range(2): # Both shapes
        for seg_idx in seg_idxs[s]: # Enumerate the segments containing crosspoint with id c_i
            seg_cross_cnt = len(seg_to_cross[s][seg_idx])
            # Get possible successor points
            # Where is the current cross?
            c_i_in_seg = seg_to_cross[s][seg_idx].index(c_i) + 2 # Index of vertex within segment

            # Get previous cross point if exists
            start_v = c_i_in_seg - 1 if c_i_in_seg != 2 else 0 # V:0 left of segment, V:1 right in segment, V:2+ intersection points
            
            refactored_path = seg_to_point(s, seg_idx, c_i_in_seg)
            if not seg_to_point(s, seg_idx, start_v) == refactored_path: # Dont add line segments consisting of one point
                possible_paths.append((s, seg_idx, (start_v, c_i_in_seg), -1)) # Vertex 2 being the first cross, 0 is first in shape order

            # Has cross point after?
            end_v = c_i_in_seg + 1 if c_i_in_seg - 1 != seg_cross_cnt else 1
            if not refactored_path == seg_to_point(s, seg_idx, end_v):
                possible_paths.append((s, seg_idx, (c_i_in_seg, end_v), 1)) # Vertex 2 being the first cross, 0 is first in shape order
    

    nbr_paths = len(possible_paths)
    paths_to_check = set(range(nbr_paths))
    valid_paths = set()
    for i in range(nbr_paths):
        for j in range(i + 1, nbr_paths):
            if j in paths_to_check and i in paths_to_check:
                s1, seg_idx1, v_idxs1, _ = possible_paths[i]
                s2, seg_idx2, v_idxs2, _ = possible_paths[j]
                if s1 != s2: 
                    if are_lines_parallel([seg_to_point(s1, seg_idx1, v_idxs1[0]), seg_to_point(s1, seg_idx1, v_idxs1[1])], [seg_to_point(s2, seg_idx2, v_idxs2[0]), seg_to_point(s2, seg_idx2, v_idxs2[1])]):
                        valid_paths.update({possible_paths[i], possible_paths[j]})
                        paths_to_check -= {i, j}
                        break

    possible_paths = list(filter(lambda p: is_contained_within(seg_to_middle_point(*p[0:3]), 
                                                               bins[(p[0] + 1) % 2], 
                                                               bounds[(p[0] + 1) % 2], 
                                                               (p[0] + 1) % 2,
                                                               cache=cache[(p[0] + 1) % 2]), 
                                                               [possible_paths[i] for i in paths_to_check]))
    valid_paths.update(possible_paths)
    possible_paths = list(valid_paths)
    
    # print("")
    # DEBUG_print_paths(list(paths_to_save))
    # DEBUG_print_paths(possible_paths, c_i)
    # #Make sure points are within both shapes
    # print("")
    # DEBUG_print_paths(possible_paths)
    return possible_paths, removed


def is_intersecting(bins, debug_correct_ans=None, plot_all=False, get_stats=False):
    start_time = time.perf_counter()
    calculate_globals(bins)
    cache = [{},{}]
    bbox, overlap_type = common_bbox(bins)
    
    if bbox == None:
        return select_return(get_stats, False, start_time)


    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    intersects, bounds = line_intersection(bins, bbox, debug_correct_ans, plot_all=plot_all, cache=cache)
    if intersects:
        return select_return(get_stats, True, start_time)

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box

    if overlap_type == '1 in 2':
        res = is_contained_within(bins[0], bins[1], bounds[1], glob_idx=1, debug_correct_ans=debug_correct_ans, plot_all=plot_all,cache=cache[1])
    elif overlap_type == '2 in 1':
        res = is_contained_within(bins[1], bins[0], bounds[0], glob_idx=0, debug_correct_ans=debug_correct_ans, plot_all=plot_all, cache=cache[0])
    else:
        res = False

    return select_return(get_stats, res, start_time)


def intersection(bins, debug_correct_ans=None, plot_all=False, get_stats=False):
    start_time = time.perf_counter()
    calculate_globals(bins)
    cache = [{},{}]
    bbox, overlap_type = common_bbox(bins)
    
    if bbox == None:
        return select_return(get_stats, Polygon(None), start_time)


    # Bounding boxes intersect. Assume no intersection, ensure that no intersection is in fact occuring:
    # 1. Find all chunks which are inside the common bounding box
    #    Construct LineStrings and check for intersections
    line_data = []
    result = line_intersection(bins, bbox, debug_correct_ans, line_data, plot_all, cache=cache)

    # 2. Ensure that the polygon is not fully contained
    #    Send ray and verify that it hits other polygon zero or even amount of times
    #    - Possibly pick point closest to other polygon's bounding box
    if len(line_data) == 0:
        _, bounds = result
        if overlap_type == '1 in 2' and is_contained_within(bins[0], bins[1], bounds[1], 1, debug_correct_ans=debug_correct_ans, plot_all=plot_all, cache=cache[1]):
            s = time.perf_counter()
            res = fpd.decompress(bins[0])[1] # Return whole smaller shape
            perf_counter[1] = perf_counter[3]
            perf_counter[0] += time.perf_counter() - s

        elif overlap_type == '2 in 1' and is_contained_within(bins[1], bins[0], bounds[0], 0, debug_correct_ans=debug_correct_ans, plot_all=plot_all, cache=cache[0]):
            s = time.perf_counter()
            res = fpd.decompress(bins[1])[1] # Return whole smaller shape
            perf_counter[2] = perf_counter[4]
            perf_counter[0] += time.perf_counter() - s
        else:
            res = Polygon(None)
        
        return select_return(get_stats, res, start_time)


    # Have intersecting points, construct resulting polygon
    # 1. Create set of intersection points, mapping: intersection point <-> line segments
    # 2. Take random intersection point from set, follow path inside both shapes
    # 3. Continue until encountering intersection point or already visited segment

    # Takes a segment and vertex_index (can be > 2 if cross-point, idx 0 is beg, idx 1 is end)
    def seg_to_point(s, seg_idx, v_idx):
        if v_idx <= 1:
            return segments[s][seg_idx][v_idx]
        else:
            return intersecting_points[seg_to_cross[s][seg_idx][v_idx - 2]]


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



    intersecting_points, segments, seg_to_cross, cross_to_seg, bounds, offset_cache = line_data
    #print("Intersecting Points:", intersecting_points)
    cross_left = set(range(len(intersecting_points))) # Ids of unprocessed intersection points
    res_segs = deque() # Segments which are part of the resulting shape
    res_points = []
    visited_edges = set()
    while len(cross_left) > 0:
        c_i = cross_left.pop() # Take one cross-point
        paths, removed = possible_paths(c_i, bounds, cross_to_seg, seg_to_cross, seg_to_point, seg_to_middle_point, bins, cache)
        if len(paths) == 0 and not removed:
            res_points.append(intersecting_points[c_i])
        while len(paths) > 0:
            path = paths.pop() # Process one path

            s, seg_idx, v_idxs, p_dir = path # V_idxs contains: start_idx, end_index. Can also be cross-points!
            
            
            path_segs = deque()
            path_append = path_segs.append
            while True: # While no cross point

                v1, v2 = seg_to_point(s, seg_idx, v_idxs[0]), seg_to_point(s, seg_idx, v_idxs[1])
                
                #Add front and back edges to visited edges set
                if (v1, v2) not in visited_edges:
                    path_append([v1, v2]) # Add segment to resulting shape a

                    # create_canvas(zoom=0.8, no_frame=True)
                    # #plot_intersecting_points()
                    # for p1, p2 in res_segs:
                    #     plot_intersecting_points([p1, p2])
                    #     plot_line(p1, p2, color='red')
                    # for p1, p2 in path_segs:
                    #     plot_intersecting_points([p1, p2])
                    #     plot_line(p1, p2, color='green')
                    # plot_geometry(fpd.decompress(bins[0])[1], solid=False, alpha=0.1)
                    # plot_geometry(fpd.decompress(bins[1])[1], solid=False, alpha=0.1)

                    visited_edges.update({(v1, v2), (v2, v1)})
                else:
                    break
                
                e_v = v_idxs[1 if p_dir == 1 else 0] # Find index of actual end vertex (i.e. flip segment based on direction)
                end_vertex = seg_to_point(s, seg_idx, e_v)
                next_seg_idx = seg_idx + p_dir
                
                
                if e_v > 1: # Is cross point?
                    break

                #For checking if next chunk has not been unfolded yet
                if next_seg_idx > -1 and next_seg_idx < len(segments[s]) and segments[s][next_seg_idx][0][0] == None:
                    #While next chunk is an non-unfolded one
                    while(next_seg_idx > -1 and next_seg_idx < len(segments[s]) and segments[s][next_seg_idx][0][0] == None):
                        follow_chunk = get_chunk(bins[s], segments[s][next_seg_idx][0][1], offset_cache=offset_cache[s], glob_idx=s) #compress that chunk
                        follow_segments = [(follow_chunk[i], follow_chunk[i + 1]) for i in range(len(follow_chunk) - 1)]  #Create  the segments
                        
                        #If there is no segments to  traverse
                        if len(follow_segments) == 0:
                            break
            
                        curr_start = follow_segments[0 if p_dir == 1 else -1][0 if p_dir == 1 else -1] #Create the start
                        if curr_start == end_vertex: #Check if current chunk is continuous to the previous one
                            end_vertex = follow_segments[-1 if p_dir == 1 else 0][-1 if p_dir == 1 else 0] #Update previous index to be this one for the next chunk/segment
                            for v1, v2 in follow_segments: #Add all segments to resulting shape
                                path_append([v1, v2])
                                visited_edges.update({(v1, v2), (v2, v1)})
                            next_seg_idx = next_seg_idx + p_dir
                        else:
                            break #Break if non continuous (Similar to the if below)

                if next_seg_idx == -1 or next_seg_idx == len(segments[s]) or not end_vertex == seg_to_point(s, next_seg_idx, 0 if p_dir == 1 else 1): #<- HERE ERROR
                    break # Break if no more segments, or if path formed by segments is not continuous
                else:
                    seg_idx = next_seg_idx # Next segment
                    

                    seg_cross_cnt = len(seg_to_cross[s][seg_idx])
                    if p_dir == 1:
                        v_idxs = [0, 2 if seg_cross_cnt != 0 else 1] # Has next line-segment crosspoint? If so take the first cross point as segment end-point.
                    else:
                        v_idxs = [seg_cross_cnt + 1 if seg_cross_cnt != 0 else 0, 1]

            # if len(path_segs) > 0:
            #     pt_size = 38
            #     create_canvas(zoom=0.8, no_frame=True)
            #     plot_geometry(fpd.decompress(bins[0])[1], solid=False, alpha=0.3, fill_alpha=0.1, hatch='/')
            #     plot_geometry(fpd.decompress(bins[1])[1], solid=False, alpha=0.3, fill_alpha=0.1, hatch='\\')
            #     for p1, p2 in res_segs:
            #         plot_intersecting_points([p1, p2], color='red', size=pt_size)
            #         plot_line(p1, p2, color='red', zorder=10)
            #     for idx, pts in enumerate(path_segs):
            #         p1, p2 = pts
                    
            #         col = (0, ((idx + 1) / len(path_segs) * 0.5 + 0.5), 0)
            #         plot_intersecting_points([p1, p2], color=col, size=pt_size)
            #         plot_line(p1, p2, color=col, zorder=10)
            #     plot_intersecting_points([intersecting_points[c_i]], color='black', size=pt_size)
            #     plt.savefig("animate/" + str(intersecting_points[0][0]) + " " + str(len(res_segs)) + ".png", bbox_inches='tight')
            #     plt.show()
            # Append path-segment to total segments
            res_segs += path_segs 

    # create_canvas(zoom=0.8, no_frame=True)
    # plot_geometry(fpd.decompress(bins[0])[1], solid=False, alpha=0.3, fill_alpha=0.1, hatch='/')
    # plot_geometry(fpd.decompress(bins[1])[1], solid=False, alpha=0.3, fill_alpha=0.1, hatch='\\')
    # for p1, p2 in res_segs:
    #     plot_intersecting_points([p1, p2], color='red', size=pt_size)
    #     plot_line(p1, p2, color='red', zorder=10)
    # plot_intersecting_points(intersecting_points, color='black', size=pt_size)
    # plt.savefig("animate/" + str(intersecting_points[0][0]) + " " + str(len(res_segs)) + ".png", bbox_inches='tight')
    # plt.show()

    #Merge all segments into LineString or MultiLineString
    unfilt_line_strs = ops.linemerge(res_segs)
    type_unfilt = unfilt_line_strs.geom_type
    unfilt_line_strs = [unfilt_line_strs] if type_unfilt == "LineString" else list(unfilt_line_strs.geoms) #For making MultiLineString and LineString be handleded similarly

    #List of all resulting linestrings and polygons
    line_strs, polygons = [], []

    #Check if LineString has the shape of a polygon, then convert it and divide polygons and line strings
    for line_str in unfilt_line_strs:
        if line_str.is_ring:
            polygons.append(Polygon(line_str.coords))
        else:
            line_strs.append(line_str)

    #Variables for different cases of geometries
    is_MultiLineString, is_MultiPolygon, is_MultiPoint = len(line_strs) > 1, len(polygons) > 1, len(res_points) > 1
    has_LineString, has_Polygon, has_Point = len(line_strs) > 0, len(polygons) > 0, len(res_points) > 0

    result = []
    if is_MultiLineString:
        result.append(MultiLineString(line_strs))
    if is_MultiPolygon:
        result.append(MultiPolygon(polygons))
    if is_MultiPoint:
        result.append(MultiPoint(res_points))
    if not is_MultiLineString and has_LineString:
        result.append(line_strs[0])
    if not is_MultiPolygon and has_Polygon:
        result.append(polygons[0])
    if not is_MultiPoint and has_Point:
        result.append(Point(res_points[0]))

    if len(result) > 1:
        res = GeometryCollection(result)
    
    elif len(result) == 1:
        res = result[0]

    return select_return(get_stats, res, start_time)