from algos.fpd_extended_lib.functions import Funcs
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import math
import numpy as np
from shapely import GeometryType as GT
from bitarray import bitarray, util, bits2bytes


class Intersection:
    ALG = None
    func = None
    def __init__(self, ALG) -> None:
        self.ALG = ALG
        self.func = Funcs(ALG)

#----------HELPER METHODS------------

    def get_bounds_intersect(self, bound_min1, bound_min2, bound_max1, bound_max2):
        if bound_min1 <= bound_min2 and bound_max1 <= bound_max2 and bound_max1 >= bound_min2:
            return (bound_min2, bound_max1)
        elif bound_min1 >= bound_min2 and bound_max1 <= bound_max2:
            return (bound_min1, bound_max1)
        if bound_min2 <= bound_min1 and bound_max2 <= bound_max1 and bound_max2 >= bound_min1:
            return (bound_min1, bound_max2)
        elif bound_min2 >= bound_min1 and bound_max2 <= bound_max1:
            return (bound_min2, bound_max2)
        else:
            return (None, None)  # In case they are not overlapping
    # @profile

    def binary_search_min_bound(self, bin, sorted_indicies, bound_right, cache, for_x):
        hi, lo = len(sorted_indicies) - 1, 0
        x_or_y = 0 if for_x else 1
        coord, cache, _ = self.func.access_vertex(bin, sorted_indicies[hi], cache)
        if coord[x_or_y] < bound_right:
            return None
        while lo < hi:
            mid = (lo + hi) // 2
            coord, cache, _ = self.func.access_vertex(bin, sorted_indicies[mid], cache)
            if bound_right <= coord[x_or_y]:
                hi = mid
            else:
                lo = mid + 1
        return lo, cache

    def binary_search_max_bound(self, bin, sorted_indicies, bound_left, cache, for_x):
        hi, lo = len(sorted_indicies) - 1, 0
        x_or_y = 0 if for_x else 1
        coord, cache, _ = self.func.access_vertex(bin, sorted_indicies[lo], cache)
        if coord[x_or_y] > bound_left:
            return None
        while lo < hi:
            mid = math.ceil((lo + hi) / 2)
            coord, cache, _ = self.func.access_vertex(bin, sorted_indicies[mid], cache)
            if coord[x_or_y] > bound_left:
                hi = mid - 1
            else:
                lo = mid
        return lo, cache

    def get_min_max_bounds(self, bin, argsorted_xs, argsorted_ys, bounds, cache):
        xs_min_idx, cache = self.binary_search_min_bound(bin, argsorted_xs, bounds[0], cache=cache, for_x=True)
        xs_max_idx, cache = self.binary_search_max_bound(bin, argsorted_xs, bounds[1], cache=cache, for_x=True)
        ys_min_idx, cache = self.binary_search_min_bound(bin, argsorted_ys, bounds[2], cache=cache, for_x=False)
        ys_max_idx, cache = self.binary_search_max_bound(bin, argsorted_ys, bounds[3], cache=cache, for_x=False)
        return cache, (xs_min_idx, xs_max_idx, ys_min_idx, ys_max_idx)

    def create_linesegments(self, bin, in_bounds_idxs, type, coord_count, cache):
        coords_linestrings = set()
        for idx in in_bounds_idxs:

            coord1, cache, (ring_beg_idx, ring_end_idx, _, _) = self.func.access_vertex(bin, idx, cache, getBoundsData=True)
            if type == GT.LINESTRING:
                if idx != coord_count - 1:
                    coord, cache, _ = self.func.access_vertex(bin, idx + 1, cache)
                    coords_linestrings.add(shapely.LineString([coord1, coord]))
                if idx != 0:
                    coord, cache, _ = self.func.access_vertex(bin, idx - 1, cache)
                    coords_linestrings.add(shapely.LineString([coord, coord1]))

            elif type != GT.LINESTRING:
                if idx == ring_beg_idx:
                    coord, cache, _ = self.func.access_vertex(bin, ring_end_idx, cache)
                    coords_linestrings.add(shapely.LineString([coord, coord1]))
                elif idx == ring_end_idx:
                    coord, cache, _ = self.func.access_vertex(bin, ring_beg_idx, cache)
                    coords_linestrings.add(shapely.LineString([coord1, coord]))
                else:
                    coord, cache, _ = self.func.access_vertex(bin, idx + 1, cache)
                    coords_linestrings.add(shapely.LineString([coord1, coord]))
                    coord, cache, _ = self.func.access_vertex(bin, idx - 1, cache)
                    coords_linestrings.add(shapely.LineString([coord1, coord]))

        return coords_linestrings, cache
    


    #----------------------------------

    def is_intersecting(self, args):
        l_bin_in, r_bin_in = args
        s = time.perf_counter()

        # Convert the binaries to bitarrays
        l_bin = bitarray(endian='big')
        r_bin = bitarray(endian='big')

        l_bin.frombytes(l_bin_in)
        r_bin.frombytes(r_bin_in)

        # Create caches for random accesing for the two binaries
        l_cache, r_cache = {}, {}

        # Retrieve the bounding boxes from the binaries
        l_bounds = list(struct.unpack_from('!dddd', l_bin_in, offset=2))  # Skip first part of header
        r_bounds = list(struct.unpack_from('!dddd', r_bin_in, offset=2))  # Skip first part of header

        # Calculate the intersecting bounding box
        bound_xmin, bound_xmax = self.get_bounds_intersect(l_bounds[0], r_bounds[0], l_bounds[2], r_bounds[2])
        bound_ymin, bound_ymax = self.get_bounds_intersect(l_bounds[1], r_bounds[1], l_bounds[3], r_bounds[3])
        bounds = (bound_xmin, bound_xmax, bound_ymin, bound_ymax)

        # If bounding boxes do not overleap, intersection can not occur
        if bound_xmin == None or bound_ymin == None:
            t = time.perf_counter()
            return t - s, False

        # Decode headers of binaries and save offset to be able to skip decode headers again
        self.ALG.offset = 0
        l_delta_size, l_type, (l_argsorted_x, l_argsorted_y), l_coord_count = self.ALG.decode_header(l_bin, True)
        l_cache['offset'] = self.ALG.offset
        l_cache['header'] = (l_delta_size, l_type)

        self.ALG.offset = 0
        r_delta_size, r_type, (r_argsorted_x, r_argsorted_y), r_coord_count = self.ALG.decode_header(r_bin, True)
        r_cache['offset'] = self.ALG.offset
        r_cache['header'] = (r_delta_size, r_type)

        # Find indicies in sorted list which are inside intersecting bounding box
        l_cache, (l_x_min_idx, l_x_max_idx, l_y_min_idx, l_y_max_idx) = self.get_min_max_bounds(l_bin, l_argsorted_x, l_argsorted_y, bounds, l_cache)
        r_cache, (r_x_min_idx, r_x_max_idx, r_y_min_idx, r_y_max_idx) = self.get_min_max_bounds(r_bin, r_argsorted_x, r_argsorted_y, bounds, r_cache)

        in_bounds_idxs1 = set(l_argsorted_x[l_x_min_idx:l_x_max_idx + 1]).intersection(set(l_argsorted_y[l_y_min_idx:l_y_max_idx + 1]))
        in_bounds_idxs2 = set(r_argsorted_x[r_x_min_idx:r_x_max_idx + 1]).intersection(set(r_argsorted_y[r_y_min_idx:r_y_max_idx + 1]))

        # Create linesegments to check pairwise intersection on
        coords1_linestrings, l_cache = self.create_linesegments(l_bin, in_bounds_idxs1, l_type, l_coord_count, l_cache)
        coords2_linestrings, r_cache = self.create_linesegments(r_bin, in_bounds_idxs2, r_type, r_coord_count, r_cache)

        for i in coords1_linestrings:
            for j in coords2_linestrings:
                if i.intersects(j):
                    t = time.perf_counter()
                    return t - s, True

        t = time.perf_counter()
        return t - s, False
    
    def intersection(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.ALG.decompress(l_bin)
        _, r_geo = self.ALG.decompress(r_bin)
        res = shapely.intersection(l_geo, r_geo)
        t = time.perf_counter()

        return t - s, res
    
    def get_non_looping_coords(self, geometry):
        coords = shapely.get_coordinates(geometry)

        if shapely.get_type_id(geometry) == GT.LINESTRING:
            return shapely.get_coordinates(geometry)

        coords = []
        if shapely.get_type_id(geometry) == GT.POLYGON:
            interiors = geometry.interiors
            coords.extend(geometry.exterior.coords[:-1])
            for ring in interiors:
                coords.extend(ring.coords[:-1])
            return coords

        elif shapely.get_type_id(geometry) == GT.MULTIPOLYGON:
            for polygon in list(geometry.geoms):
                interiors = polygon.interiors
                coords.extend(polygon.exterior.coords[:-1])
                for ring in interiors:
                    coords.extend(ring.coords[:-1])
            return coords
        
    def append_header(self, bits, geometry, d_size):
        # Meta data
        bits.frombytes(self.uchar_to_bytes(d_size))
        bits.frombytes(self.uchar_to_bytes(int(shapely.get_type_id(geometry))))  # 1 byte is enough for storing type
        # Bounding Box
        bounds = shapely.bounds(geometry)
        bits.frombytes(self.double_to_bytes(bounds[0]))
        bits.frombytes(self.double_to_bytes(bounds[1]))
        bits.frombytes(self.double_to_bytes(bounds[2]))
        bits.frombytes(self.double_to_bytes(bounds[3]))

        coords = self.get_non_looping_coords(geometry)
        bits.extend(self.uint_to_ba(int(len(coords)), 4 * 8))  # size of integer
        sorted_idxs = [np.argsort([coord[0] for coord in coords]), np.argsort([coord[1] for coord in coords])]
        idx_bits = math.ceil(math.log2(len(coords)))
        for i in range(2):
            for idx in range(len(coords)):
                bits.extend(self.uint_to_ba(int(sorted_idxs[i][idx]), idx_bits))

    def decode_header(self, bin, get_idxs=False):
        delta_size, type = struct.unpack_from('!BB', bin)
        type = GT(type)
        self.offset += 2 * 8 + 4 * 64  # Offset is 2 bytes for BB + 64 * 4 for bounding box

        # Code segment needed for extracting sorted indexes
        coord_count = self.bytes_to_uint(bin, 4 * 8)
        idx_sizes = math.ceil(math.log2(coord_count))
        sorted_idxs = [[], []]
        if get_idxs:
            for i in range(2):
                for _ in range(coord_count):
                    sorted_idxs[i].append(self.bytes_to_uint(bin, idx_sizes))
            return delta_size, type, sorted_idxs, coord_count
        else:
            self.offset += idx_sizes * 2 * coord_count
        return delta_size, type
