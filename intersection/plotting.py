import shapely
import numpy as np
from matplotlib import pyplot as plt
import random
random.seed(13)


def xs(pts): return [p[0] for p in pts]
def ys(pts): return [p[1] for p in pts]


def inv_color(shp): return (1 - color(shp)[0], 1 - color(shp)[1], 1 - color(shp)[2])


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


def plot_intersecting_points(pts, SHOW_INTERSECTING_POINTS=True):
    if SHOW_INTERSECTING_POINTS:
        plt.scatter(xs(pts), ys(pts), 22, zorder=10)


# FPDE related stuff
from algos.alg_fpd_extended import FpdExtended
from bitarray import bitarray
from shapely import GeometryType as GT
fpd = FpdExtended()


def calculate_chunks_bounds(bin):
    chunks = fpd.get_chunks(bin)
    chunks_bounds = []
    chunks_vertices = []
    for idx, chunk in enumerate(chunks):
        xs, ys = [coord[0] for coord in chunk], [coord[1] for coord in chunk]
        xs.append(chunks[(idx + 1) % len(chunks)][0][0])
        ys.append(chunks[(idx + 1) % len(chunks)][0][1])
        chunks_bounds.append([min(xs), min(ys), max(xs), max(ys)])
        chunks_vertices.append((xs, ys))
    return chunks_bounds, chunks_vertices


def plot_chunks_bounds(bin_in, include_next_chunk_start=False, avoid_create_frame=False, avoid_show=False, idxs=None):
    if not avoid_create_frame:
        fig = plt.figure()
        zoom = 2.5
        w, h = fig.get_size_inches()
        fig.set_size_inches(w * zoom, h * zoom)

    chunks_bounds, chunks_vertices = calculate_chunks_bounds(bin_in)
    if idxs != None:
        # Filter out based on idxs
        chunks_bounds = [chunks_bounds[i] for i in idxs]
        chunks_vertices = [chunks_vertices[i] for i in idxs]

    for chunk_data in [(chunks_bounds[i], chunks_vertices[i]) for i in range(len(chunks_bounds))]:
        chunk_bounds, chunk_vertices = chunk_data
        xs, ys = chunk_vertices
        x_l, y_b, x_r, y_t = chunk_bounds
        x_bounds = [x_l, x_r, x_r, x_l, x_l]
        y_bounds = [y_t, y_t, y_b, y_b, y_t]
        chunk_color = (np.random.random(), np.random.random(), np.random.random())
        inverse_chunk_color = (1 - chunk_color[0], 1 - chunk_color[1], 1 - chunk_color[2])
        plt.plot(x_bounds, y_bounds, color=chunk_color)
        plt.fill(x_bounds, y_bounds, color=chunk_color, alpha=0.05)
        plt.scatter(xs, ys, s=10, color=inverse_chunk_color)

    _, vertices = fpd.vertices(bin_in)  # Avoid problem with ring not connecting in end
    all_x, all_y = [coord[0] for coord in vertices], [coord[1] for coord in vertices]
    plt.plot(all_x, all_y, zorder=-1)
    plt.title("Chunk Bounds" if not include_next_chunk_start else "Chunk Bounds - With Connecting Borders")
    if not avoid_show:
        plt.show()
