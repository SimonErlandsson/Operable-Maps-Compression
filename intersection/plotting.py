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

def bounds_to_coords(bounds):
    x_l, y_b, x_r, y_t = bounds
    x_bounds = [x_l, x_r, x_r, x_l, x_l]
    y_bounds = [y_t, y_t, y_b, y_b, y_t]
    return x_bounds, y_bounds

def append_polygon(shp, rings):
    exterior_coords = list(shp.exterior.coords)
    interior_coords = [list(interior.coords) for interior in shp.interiors]
    rings += [(True, ring) for ring in ([exterior_coords] + interior_coords)]


def plot_geometry(geom, SHOW_GEOMETRIES=True, solid=True, alpha=1.0, fill_alpha=0.01):
    if SHOW_GEOMETRIES:
        geom_type = shapely.get_type_id(geom)
        rings = []
        if geom_type == GT.LINESTRING or geom_type == GT.POINT or geom_type == GT.MULTILINESTRING or geom_type == GT.MULTIPOINT:
            pts = shapely.get_coordinates(geom)
            rings = [(False, pts)]
        else:
            shps = [geom] if shapely.get_type_id(geom) == GT.POLYGON else geom.geoms
            for shp in shps:
                shp_type = shapely.get_type_id(shp)
                if shp_type == GT.LINESTRING or shp_type == GT.POINT or shp_type == GT.MULTILINESTRING or shp_type == GT.MULTIPOINT:
                    pts = shapely.get_coordinates(geom)
                    rings += [(False, pts)]
                elif shp_type == GT.POLYGON:
                    append_polygon(shp, rings)
                elif shp_type == GT.MULTIPOLYGON:
                    for s in shp.geoms: # Go through polygons in multipoly
                        append_polygon(s, rings)
        for ring_solid, ring_points in rings:
            if ring_solid and solid:
                plt.fill(xs(ring_points), ys(ring_points), color=color(geom), alpha=fill_alpha)
            plt.plot(xs(ring_points), ys(ring_points), '-' if solid else '--', color=color(geom), alpha=alpha)


def plot_geometry_bbox(geom, SHOW_BOUNDING_BOXES=True, solid=False):
    bbox = bbox_coords(geom)
    if SHOW_BOUNDING_BOXES:
        plt.plot(xs(bbox), ys(bbox), '-' if solid else '--', color=color(geom), zorder=-10)

def plot_bounds(bounds, solid=False, color=None, zorder=20, alpha=1.0):
    if color == None:
        color = (np.random.random(), np.random.random(), np.random.random())
    x_coords, y_coords = bounds_to_coords(bounds)
    plt.plot(x_coords, y_coords, '-' if solid else '--', color=color, zorder=zorder, alpha=alpha)


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

def plot_raw_points(pts, color='blue', size=22, alpha=1.0):
    plt.scatter(xs(pts), ys(pts), size, color=color, zorder=20, alpha=alpha)

def plot_raw_point(pt, color='blue', size=22, alpha=1.0):
    plt.scatter(pt[0], pt[1], size, color=color, zorder=20, alpha=alpha)

def plot_line(pt_1, pt_2, color='blue', solid=True, alpha=1.0):
    coords = [pt_1, pt_2]
    plt.plot(xs(coords), ys(coords), '-' if solid else '--', color=color, alpha=alpha)


# FPDE related stuff
from algos.alg_fpd_extended import FpdExtended
from bitarray import bitarray
from shapely import GeometryType as GT
fpd = FpdExtended()


def calculate_chunks_bounds(bin, include_next_chunk_start=True):
    chunks, is_last_chunk_ring = fpd.get_chunks(bin, include_next_chunk_start)
    _, type = fpd.type(bin)
    chunks_bounds = []
    chunks_vertices = []
    chk_idx_ring_start = 0
    for idx, chunk in enumerate(chunks):
        xs, ys = [coord[0] for coord in chunk], [coord[1] for coord in chunk]
        if not is_last_chunk_ring[idx] and type != 'LineString' and include_next_chunk_start:
            xs.append(chunks[idx + 1][0][0])
            ys.append(chunks[idx + 1][0][1])
        chunks_bounds.append([min(xs), min(ys), max(xs), max(ys)])
        chunks_vertices.append((xs, ys))
    return chunks_bounds, chunks_vertices

def create_canvas(zoom=2.5, no_frame=False):
    fig = plt.figure()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * zoom, h * zoom)
    if no_frame:
        plt.axis('off')
        plt.margins(0.01)

def plot_chunks_bounds(bin_in, include_next_chunk_start=False, avoid_create_frame=False, avoid_show=False, idxs=None, txt='', solid=True, alpha=1.0, fill_alpha=0.02):
    if not avoid_create_frame:
        create_canvas()

    chunks_bounds, chunks_vertices = calculate_chunks_bounds(bin_in, include_next_chunk_start)
    if idxs != None:
        # Filter out based on idxs
        chunks_bounds = [chunks_bounds[i] for i in idxs]
        chunks_vertices = [chunks_vertices[i] for i in idxs]

    for chunk_data in [(chunks_bounds[i], chunks_vertices[i]) for i in range(len(chunks_bounds))]:
        chunk_bounds, chunk_vertices = chunk_data
        xs, ys = chunk_vertices
        x_bounds, y_bounds = bounds_to_coords(chunk_bounds)
        chunk_color = (np.random.random(), np.random.random(), np.random.random())
        inverse_chunk_color = (1 - chunk_color[0], 1 - chunk_color[1], 1 - chunk_color[2])
        plot_bounds(chunk_bounds, color=chunk_color, solid=solid, alpha=alpha)
        plt.fill(x_bounds, y_bounds, color=chunk_color, alpha=fill_alpha)
        plt.scatter(xs, ys, s=20, zorder=30, color=inverse_chunk_color)

    _, geom = fpd.decompress(bin_in)
    #plot_geometry(geom, alpha=0.2)

    #plt.title(("Chunk Bounds" if not include_next_chunk_start else "Chunk Bounds - With Connecting Borders") + str(txt))
    if not avoid_show:
        plt.show()
