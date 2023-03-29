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