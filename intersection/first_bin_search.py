#Helpers
import numpy as np
import shapely
import matplotlib.pyplot as plt
import math

from intersection.plotting import create_canvas, plot_common_bbox, plot_coordinates, plot_geometry, plot_geometry_bbox, plot_intersecting_points

'''

    NOTE: This method was abandoned for reasons outlined in the report.

'''

def get_bounds_intersect(bound_min1, bound_min2, bound_max1, bound_max2):
    if bound_min1 <= bound_min2 and bound_max1 <= bound_max2 and bound_max1 >= bound_min2:
        return (bound_min2, bound_max1)
    elif bound_min1 >= bound_min2 and bound_max1 <= bound_max2:
        return (bound_min1, bound_max1)
    if bound_min2 <= bound_min1 and bound_max2 <= bound_max1 and bound_max2 >= bound_min1:
        return (bound_min1, bound_max2)
    elif bound_min2 >= bound_min1 and bound_max2 <= bound_max1:
        return (bound_min2, bound_max2)
    else:
        return (None, None)

def binary_search_min_bound(sorted_indicies, point_list, bound_right):
    hi, lo = len(sorted_indicies) - 1, 0
    if point_list[sorted_indicies[hi]] < bound_right:
        return None
    while lo < hi:
        mid = (lo+hi)//2
        if bound_right <= point_list[sorted_indicies[mid]]: hi = mid
        else: lo = mid+1
    return lo

def binary_search_max_bound(sorted_indicies, point_list, bound_left):
    hi, lo = len(sorted_indicies) - 1, 0
    if point_list[sorted_indicies[lo]] > bound_left:
        return None
    while lo < hi:
        mid = math.ceil((lo+hi)/2)
        if point_list[sorted_indicies[mid]] > bound_left: hi = mid - 1
        else: lo = mid
    return lo

def divide_xs_ys(geom):
    xs, ys = [], []
    coords = shapely.get_coordinates(geom)[:-1]
    for coord in coords:
        xs.append(coord[0])
        ys.append(coord[1])
    return coords, xs, ys

def get_min_max_bounds(argsorted_xs, argsorted_ys, xs, ys, bounds):
    xs_min_idx = binary_search_min_bound(argsorted_xs, xs, bounds[0])
    xs_max_idx = binary_search_max_bound(argsorted_xs, xs, bounds[1])

    ys_min_idx = binary_search_min_bound(argsorted_ys, ys, bounds[2])
    ys_max_idx = binary_search_max_bound(argsorted_ys, ys, bounds[3])
    return xs_min_idx, xs_max_idx, ys_min_idx, ys_max_idx

def create_linesegments(in_bounds_idxs, coords):
    coords_linestrings = set()
    for i in in_bounds_idxs:
        if i == 0:
            coords_linestrings.add(shapely.LineString([coords[0], coords[1]]))
        elif i == len(coords) - 1:
            coords_linestrings.add(shapely.LineString([coords[len(coords) - 1], coords[0]]))
        else:
            coords_linestrings.add(shapely.LineString([coords[i - 1], coords[i]]))
            coords_linestrings.add(shapely.LineString([coords[i], coords[i + 1]]))
    return coords_linestrings


def binary_intersection(geom1, geom2):
    coords1, xs1, ys1 = divide_xs_ys(geom1)
    geom1_bounds = geom1.bounds

    coords2, xs2, ys2 = divide_xs_ys(geom2)
    geom2_bounds = geom2.bounds

    argsorted_xs1, argsorted_ys1 = np.argsort(xs1), np.argsort(ys1)
    argsorted_xs2, argsorted_ys2 = np.argsort(xs2), np.argsort(ys2)

    bound_xmin, bound_xmax = get_bounds_intersect(geom1_bounds[0], geom2_bounds[0], geom1_bounds[2], geom2_bounds[2])
    bound_ymin, bound_ymax = get_bounds_intersect(geom1_bounds[1], geom2_bounds[1], geom1_bounds[3], geom2_bounds[3])
    bounds = (bound_xmin, bound_xmax, bound_ymin, bound_ymax)

    if bound_xmin == None or bound_ymin == None:
        print("Geometries cant intersect since bounding boxes are non intersecting")
        return False, []
    else:
        xs1_min_idx, xs1_max_idx, ys1_min_idx, ys1_max_idx = get_min_max_bounds(argsorted_xs1, argsorted_ys1, xs1, ys1, bounds)
        xs2_min_idx, xs2_max_idx, ys2_min_idx, ys2_max_idx = get_min_max_bounds(argsorted_xs2, argsorted_ys2, xs2, ys2, bounds)

        in_bounds_idxs1 = set(argsorted_xs1[xs1_min_idx:xs1_max_idx + 1]).intersection(set(argsorted_ys1[ys1_min_idx:ys1_max_idx + 1]))
        in_bounds_idxs2 = set(argsorted_xs2[xs2_min_idx:xs2_max_idx + 1]).intersection(set(argsorted_ys2[ys2_min_idx:ys2_max_idx + 1]))

        ## Add left and right segs
        for i in in_bounds_idxs1.copy():
            left = (i - 1) % len(coords1)
            right = (i + 1) % len(coords1)
            in_bounds_idxs1.add(left)
            in_bounds_idxs1.add(right)
        for i in in_bounds_idxs2.copy():
            left = (i - 1) % len(coords2)
            right = (i + 1) % len(coords2)
            in_bounds_idxs2.add(left)
            in_bounds_idxs2.add(right)

        coords1_linestrings = create_linesegments(in_bounds_idxs1, coords1)  
        coords2_linestrings = create_linesegments(in_bounds_idxs2, coords2)

        # for i in in_bounds_idxs2:
        #     if i == 0:
        #         coords2_linestrings.add(shapely.LineString([coords2[0], coords2[1]]))
        #     elif i == len(coords2) - 1:
        #         coords2_linestrings.add(shapely.LineString([coords2[len(coords2) - 2], coords2[len(coords2) - 1]]))
        #     else:
        #         coords2_linestrings.add(shapely.LineString([coords2[i - 1], coords2[i]]))
        #         coords2_linestrings.add(shapely.LineString([coords2[i], coords2[i + 1]]))
        
        intersection_points_x = []
        intersection_points_y = []
        for i in coords1_linestrings:
            for j in coords2_linestrings:
                if i.intersects(j):
                    int_point = i.intersection(j)
                    intersection_points_x.append(int_point.x)
                    intersection_points_y.append(int_point.y)


        # points_in_bounds_idxs1 = [coords1[i] for i in in_bounds_idxs1]
        # points_in_bounds_idxs2 = [coords2[i] for i in in_bounds_idxs2]
        # points_not_in_bounds_idxs1 = [c for idx, c in enumerate(coords1) if idx not in in_bounds_idxs1]
        # points_not_in_bounds_idxs2 = [c for idx, c in enumerate(coords2) if idx not in in_bounds_idxs2]
        # create_canvas(no_frame=True, zoom=1.25)
        # for s, g in enumerate([geom1, geom2]):
        #     #plot_chunks_bounds(bins[s], include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, idxs=chunks_not_loaded, solid=False, alpha=0.3, fill_alpha=0.05, color='black', point_color='cyan', point_size=9)
        #     #plot_chunks_bounds(bins[s], include_next_chunk_start=True, avoid_show=True, avoid_create_frame=True, idxs=chk_idxs[s], alpha=1.0, color='orange', point_color='red', point_size=14)
        #     plot_geometry(g, fill_alpha=0.1, hatch=('\\' if s == 1 else '/'), color=('blue' if s == 1 else 'green'))
        #     plot_geometry_bbox(g, solid=False, color=('blue' if s == 1 else 'green'), linewidth=1)
        #     #For excp cases: plot_coordinates(g, size=40, color=('blue' if s == 1 else 'green'), zorder=100)
        # plot_common_bbox([geom1, geom2], color='red', linewidth=1, zorder=90)
        # from shapely import intersection as s_inter
        # inter_shape = s_inter(geom1, geom2)
        # intersecting_points = [[intersection_points_x[i], intersection_points_y[i]] for i in range(len(intersection_points_x))]
        # plot_geometry(inter_shape, fill_alpha=0.3, hatch='X', color='purple')

        # inter_coords = shapely.get_coordinates(inter_shape)
        # inter_shape_added = list(filter(lambda x: x not in intersecting_points, list(map(lambda y: list(y), inter_coords))))
        # points_in_bounds_idxs1 = list(filter(lambda x: x not in inter_coords, points_in_bounds_idxs1))
        # points_in_bounds_idxs2 = list(filter(lambda x: x not in inter_coords, points_in_bounds_idxs2))

        # plot_intersecting_points(inter_shape_added[:-3], color='red', zorder=100, size=50)
        
        # plot_intersecting_points(inter_shape_added, color='lime', zorder=100, size=30)
        # plot_intersecting_points(points_in_bounds_idxs1, color='red', zorder=100, size=22)
        # plot_intersecting_points(points_in_bounds_idxs2, color='red', zorder=100, size=22)
        # plot_intersecting_points(points_not_in_bounds_idxs1, color='cyan', zorder=100, size=22)
        # plot_intersecting_points(points_not_in_bounds_idxs2, color='cyan', zorder=100, size=22)
        # plot_intersecting_points(intersecting_points, color='black', zorder=100, size=30)
        #plt.show()

        SHOW_COORDINATES = False
        SHOW_GEOMETRIES = False
        SHOW_INTERSECTIONS = False
        SHOW_BOUNDING_INTERSECTION = False
        SHOW_BOUNDING_BOXES = False
            
        legends = []
        if SHOW_BOUNDING_BOXES:
            plt.plot([geom1_bounds[0], geom1_bounds[0],geom1_bounds[2],geom1_bounds[2], geom1_bounds[0]], [geom1_bounds[1], geom1_bounds[3], geom1_bounds[3], geom1_bounds[1],geom1_bounds[1]])
            legends.append("Geometry 1 bounding box")
            plt.plot([geom2_bounds[0], geom2_bounds[0],geom2_bounds[2],geom2_bounds[2],geom2_bounds[0]], [geom2_bounds[1], geom2_bounds[3], geom2_bounds[3], geom2_bounds[1],geom2_bounds[1]])
            legends.append("Geometry 2 bounding box")

        if SHOW_BOUNDING_INTERSECTION:
            plt.plot([bound_xmin, bound_xmin,bound_xmax,bound_xmax,bound_xmin], [bound_ymin, bound_ymax,bound_ymax, bound_ymin,bound_ymin])
            legends.append("Bounding box intersection")

        if SHOW_COORDINATES:
            points_in_bounds_idxs1 = [coords1[i] for i in in_bounds_idxs1]
            points_in_bounds_idxs2 = [coords2[i] for i in in_bounds_idxs2]
            x_coords_geom1 = [point[0] for point in points_in_bounds_idxs1]
            y_coords_geom1 = [point[1] for point in points_in_bounds_idxs1]

            x_coords_geom2 = [point[0] for point in points_in_bounds_idxs2]
            y_coords_geom2 = [point[1] for point in points_in_bounds_idxs2]
        
            plt.scatter(x_coords_geom1, y_coords_geom1, zorder=10)
            legends.append("Geometry 1 coordinates")
            plt.scatter(x_coords_geom2, y_coords_geom2, zorder=10)
            legends.append("Geometry 2 coordinates")
        
        if SHOW_INTERSECTIONS:
            plt.scatter(intersection_points_x, intersection_points_y, zorder=10)
            legends.append("Intersection points")

        if SHOW_GEOMETRIES:
            x1, y1 = geom1.exterior.xy
            x2, y2 = geom2.exterior.xy
            plt.fill(x1, y1, alpha=0.1)
            plt.plot(x1, y1)
            legends.append("Geometry 1")
            plt.fill(x2, y2, alpha=0.1)
            plt.plot(x2, y2)
            legends.append("Geometry 2")

        if SHOW_COORDINATES or SHOW_GEOMETRIES or SHOW_INTERSECTIONS or SHOW_BOUNDING_INTERSECTION or SHOW_BOUNDING_BOXES:
            plt.legend(legends, bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)
            plt.title("Intersection Plot: " + ('False' if len(intersection_points_x) == 0 else 'True'))
            plt.tight_layout
            plt.show()
            
        return len(intersection_points_x) != 0, [(intersection_points_x[i], intersection_points_y) for i in range(len(intersection_points_x))]