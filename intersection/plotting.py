import shapely
import numpy as np
from matplotlib import pyplot as plt
import random
random.seed(13)

xs = lambda ps: [p[0] for p in ps]
ys = lambda ps: [p[1] for p in ps]

inv_color = lambda shp: (1 - color(shp)[0], 1 - color(shp)[1], 1 - color(shp)[2])

colors = {}
def color(shp):
    wkt = shapely.to_wkt(shp)
    if wkt not in colors.keys():
        colors[wkt] = [random.random() for _ in range(3)]
    return colors[wkt]

def bbox_coords(geom):
    x_l, y_b, x_r, y_t = shapely.bounds(geom)
    return [(x_l, y_b), (x_l, y_t), (x_r, y_t), (x_r, y_b), (x_l, y_b)]

def plot_geometry(geom, SHOW_GEOMETRIES=True):
    ps = shapely.get_coordinates(geom)
    if SHOW_GEOMETRIES:
        plt.fill(xs(ps), ys(ps), color=color(geom), alpha=0.1)
        plt.plot(xs(ps), ys(ps), color=color(geom))

def plot_geometry_bbox(geom, SHOW_BOUNDING_BOXES=True, solid=False):
    bbox = bbox_coords(geom)
    if SHOW_BOUNDING_BOXES:
        plt.plot(xs(bbox), ys(bbox), '-' if solid else '--', color=color(geom), zorder=-10)

def plot_common_bbox(geometries, SHOW_COMMON_BOUNDING_BOX=True):
    bbox_shape = shapely.intersection(shapely.Polygon(bbox_coords(geometries[0])), shapely.Polygon(bbox_coords(geometries[1])))
    if SHOW_COMMON_BOUNDING_BOX:
        plot_geometry_bbox(bbox_shape, solid=True)

def plot_coordinates(geom, SHOW_COORDINATES=True):
    ps = shapely.get_coordinates(geom)
    if SHOW_COORDINATES:
        plt.scatter(xs(ps), ys(ps), 13, zorder=10, color=inv_color(geom))
    
def plot_intersecting_points(ps, SHOW_INTERSECTING_POINTS=True):
    if SHOW_INTERSECTING_POINTS:
        plt.scatter(xs(ps), ys(ps), zorder=10)

# FPDE related stuff
from algos.alg_fpd_extended import FpdExtended
from bitarray import bitarray
from shapely import GeometryType as GT
fpd = FpdExtended()


def plot_chunk_bounds(bin_in, include_next_chunk_start=False):
    fig = plt.figure()
    zoom = 2.5
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * zoom, h * zoom)

    chunks = fpd.get_chunks(bin_in)
    for idx, chunk in enumerate(chunks):
        xs, ys = [coord[0] for coord in chunk], [coord[1] for coord in chunk]
        if include_next_chunk_start:
            xs.append(chunks[(idx + 1) % len(chunks)][0][0])
            ys.append(chunks[(idx + 1) % len(chunks)][0][1])
        x_bounds = [min(xs), max(xs), max(xs), min(xs), min(xs)]
        y_bounds = [max(ys), max(ys), min(ys), min(ys), max(ys)]
        chunk_color = (np.random.random(), np.random.random(), np.random.random())
        inverse_chunk_color = (1 - chunk_color[0], 1 - chunk_color[1], 1 - chunk_color[2])
        plt.plot(x_bounds, y_bounds, color=chunk_color)
        plt.fill(x_bounds, y_bounds, color=chunk_color, alpha=0.05)
        plt.scatter(xs, ys, s=10, color=inverse_chunk_color)
    
    _, vertices = fpd.vertices(bin_in) # Avoid problem with ring not connecting in end 
    all_x, all_y = [coord[0] for coord in vertices], [coord[1] for coord in vertices]
    plt.plot(all_x, all_y, zorder=-1)
    plt.title("Chunk Bounds" if not include_next_chunk_start else "Chunk Bounds - With Connecting Borders")
    plt.show()