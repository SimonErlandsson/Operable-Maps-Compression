
from base import CompressionAlgorithm
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import sys
import math
from shapely import GeometryType as GT


class Fpd(CompressionAlgorithm):
    MAX_CHUNK_SIZE = 10
    OPTIMAL_DELTA_SIZE = 4  # DEFAULT VALUEß

# ---- HELPER METHODS
    def binary(self, num):
        return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

    def int_repr(self, num):
        return struct.unpack('!i', struct.pack('!f', num))[0]

    def int64_to_int16(self, num):
        return int.from_bytes(struct.pack('!h', num), "big")

    def float_to_bytes(self, x):
        return struct.pack("f", x)

    def int_to_bytes(self, x):
        return struct.pack('!I', x)

    def zz_encode(self, num):
        return -2 * num - 1 if num < 0 else 2 * num

    def get_zz_encoded_delta(self, prev_coord, curr_coord):
        return self.zz_encode(self.int_repr(curr_coord) - self.int_repr(prev_coord))

    def deltas_fit_in_bytes(self, max_bytes, delta_x, delta_y):
        return (delta_x == 0 or math.log2(delta_x) <= max_bytes * 8) and (delta_y == 0 or math.log2(delta_y) <= max_bytes * 8)

    def get_poly_ring_count(self, geometry):
        subpoly_coord_count = deque([])
        if geometry.geom_type == "Polygon":
            subpoly_coord_count.append(len(geometry.exterior.coords))
            for i in range(len(geometry.interiors)):
                subpoly_coord_count.append(len(geometry.interiors[i].coords))

        elif geometry.geom_type == "MultiPolygon":
            for polygon in list(geometry.geoms):
                subpoly_coord_count.append(len(list(polygon.exterior.coords)))
                for i in range(len(polygon.interiors)):
                    subpoly_coord_count.append(len(polygon.interiors[i].coords))
        return subpoly_coord_count
    

    def get_multipoly_count(self, geometry):
        poly_coord_count = deque([])
        if geometry.geom_type == "MultiPolygon":
            poly_coord_count = deque([len(shapely.get_coordinates(polygon)) for polygon in list(geometry.geoms)])
        #Flattens array
        return poly_coord_count
    
    def init_chunk(self, res, chk_size, coords, pnt_idx, chk_ys, chk_xs):
        chk_size_idx = len(res)
        res.append(self.int_to_bytes(chk_size))

        # Add full coordinates
        chk_xs.append(self.float_to_bytes(coords[pnt_idx][0]))
        chk_ys.append(self.float_to_bytes(coords[pnt_idx][1]))
        return chk_size_idx

    def save_clear_state(self, res, chk_xs, chk_ys):
        res += chk_xs + chk_ys
        # Reset current chk state
        chk_xs.clear()
        chk_ys.clear()

    def append_delta_pair(self, chk_xs, chk_ys, deltax, deltay, delta_size):
        delta_bytes_x = deltax.to_bytes(delta_size, 'big')
        delta_bytes_y = deltay.to_bytes(delta_size, 'big')
        chk_xs.append(delta_bytes_x)
        chk_ys.append(delta_bytes_y)


    def init_bytearray(self, geometry, delta_size):
        res = []
        res.append(self.int_to_bytes(delta_size))
        res.append(int(shapely.get_type_id(geometry)).to_bytes(1, 'big')) #1 byte is enought for storing type
        return res
    
    def decode_header(self, bin):
        delta_size, type = struct.unpack_from('!IB', bin)
        type = GT(type)
        return delta_size, type

    def polygon_reset(self, res, chk_xs, chk_ys, chk_size_idx, geom_coord_count, chk_point_nbr):
        #Reset and store previous chk data if interrupted in the middle of a chk
        if(len(chk_xs) != 0):
            res[chk_size_idx] = self.int_to_bytes(chk_point_nbr)
            chk_point_nbr = 0
            self.save_clear_state(res, chk_xs, chk_ys)

        subpoly_point_cnt = geom_coord_count.popleft()
        res.append(self.int_to_bytes(subpoly_point_cnt))
        return chk_point_nbr, subpoly_point_cnt

        
    def fp_delta_encoding(self, geometry, delta_size, max_chk_size):

        coords = shapely.get_coordinates(geometry)

        #Saves interior and exterior lengths if polygon type
        subpoly_coord_count:deque = self.get_poly_ring_count(geometry)
        poly_coord_count:deque = self.get_multipoly_count(geometry)

    
        #State of the current chk state
        chk_xs, chk_ys = [], []
        chk_size_idx = 0  #Used for changing the chk_size of reset occurs

        #State of the current coordiate pair
        chk_point_nbr = 0
        
        #Array for storing list resulting bytes. Initialized with Total coordinates, bounding box, delta_size, chk_size, geom_type
        res = self.init_bytearray(geometry,  delta_size)

        #Polygon specific variables        
        is_multipolygon = len(poly_coord_count) != 0
        is_polygon = len(subpoly_coord_count) != 0
        
        poly_point_cnt = 0
        subpoly_point_cnt = 0

        #Loop all coordinates
        for i in range(len(coords)):
            if is_multipolygon and poly_point_cnt == 0:
               chk_point_nbr, poly_point_cnt = self.polygon_reset(res, chk_xs, chk_ys, chk_size_idx, poly_coord_count, chk_point_nbr)


            if is_polygon and subpoly_point_cnt == 0:
               chk_point_nbr, subpoly_point_cnt = self.polygon_reset(res, chk_xs, chk_ys, chk_size_idx, subpoly_coord_count, chk_point_nbr)

            #----CHUNK HEADER INFORMATION (CHUNK SIZE, FULL FIRST COORDINATES)
            if chk_point_nbr == 0:
                #save chk size and later change if reset occurs
                chk_size_idx = self.init_chunk(res, self.MAX_CHUNK_SIZE, coords, i, chk_ys, chk_xs)

            else: #Loop for delta
                zig_delta_x = self.get_zz_encoded_delta(coords[i - 1][0],coords[i][0])
                zig_delta_y = self.get_zz_encoded_delta(coords[i - 1][1],coords[i][1])
                if self.deltas_fit_in_bytes(delta_size,zig_delta_x,zig_delta_y):
                    self.append_delta_pair(chk_xs, chk_ys, zig_delta_x, zig_delta_y, delta_size)
                    chk_point_nbr += 1
                    
                else:
                    res[chk_size_idx] = self.int_to_bytes(chk_point_nbr)
                    self.save_clear_state(res, chk_xs, chk_ys)
                    chk_size_idx = self.init_chunk(res, self.MAX_CHUNK_SIZE, coords, i, chk_ys, chk_xs)

                    chk_point_nbr = 1

                if chk_point_nbr == max_chk_size:
                    chk_point_nbr = 0
                    self.save_clear_state(res, chk_xs, chk_ys)

            subpoly_point_cnt -= 1
            poly_point_cnt -= 1


        if chk_point_nbr > 0:
            self.save_clear_state(res, chk_xs, chk_ys)
            res[chk_size_idx] = self.int_to_bytes(chk_point_nbr)

        res = b''.join(res)
        return res, len(res)

    def fp_delta_decoding(self, bin):
        geometry = shapely.LineString()
        delta_size, type = self.decode_header(bin)
        
        if type == GT.LINESTRING:
            pass
        elif type == GT.POLYGON:
            pass
        elif type == GT.MULTIPOLYGON:
            pass


        # Remove header information about total nbr of nodes and bounding box
        bin = bin[4 + 4 * 2:]
        delta_size = int.from_bytes(bin[0:4], byteorder='big')
        reset_next = True
        chk_deltas = 0


        # Parse to shapely object
        # first retrieve the coordinates from your geojson dictionary
        #rawcoords = candidate["coordinates"]

        # # define the conversion functionß
        # def PrepCoordsForShapely(rawcoords):
        #     preppedcoords = []
        #     #according to the geojson specs, a multipolygon is a list of linear rings, so we loop each
        #     for eachpolygon in rawcoords:
        #         #the first linear ring is the coordinates of the polygon, and shapely needs it to be a tuple
        #         tupleofcoords = tuple(eachpolygon[0])
        #         #the remaining linear rings, if any, are the coordinates of inner holes, and shapely needs these to be nested in a list
        #         if len(eachpolygon) > 1:
        #             listofholes = list(eachpolygon[1:])
        #         else:
        #             listofholes = []
        #         #shapely defines each polygon in a multipolygon with the polygoon coordinates and the list of holes nested inside a tuple
        #         eachpreppedpolygon = (tupleofcoords, listofholes)
        #         #so append each prepped polygon to the final multipolygon list
        #         preppedcoords.append(eachpreppedpolygon)
        #     #finally, the prepped coordinates need to be nested inside a list in order to be used as a star-argument for the MultiPolygon constructor.
        #     return [preppedcoords]

        # # use the function to prepare coordinates for MultiPolygon
        # preppedcoords = PrepCoordsForShapely(rawcoords)
        # # use the prepped coordinates as a star-argument for the MultiPolygon constructor
        # shapelymultipolygon = shapely.MultiPolygon(*preppedcoords)

        return geometry



    def calculate_delta_size(self, geometry):
         coords = shapely.get_coordinates(geometry)
         prev = [0, 0]
         bit_cnts = {}
         for coord in coords:
             for i in range(2):
                 delta = self.get_zz_encoded_delta(prev[i], coord[i])
                 bit_cnt = 1 if delta == 0 else math.ceil(
                     math.log2(delta))
                 if bit_cnt not in bit_cnts:
                     bit_cnts[bit_cnt] = 1
                 else:
                     bit_cnts[bit_cnt] += 1
             prev = coord
         bit_cnts = dict(sorted(bit_cnts.items(), reverse=True))

         tot_size = {}
         upper_cnt = 0
         lower_cnt = 2 * len(coords)
         for n in bit_cnts.keys():
             tot_size[n] = n * lower_cnt + 32 * upper_cnt
             lower_cnt -= bit_cnts[n]
             upper_cnt += bit_cnts[n]

         self.OPTIMAL_DELTA_SIZE = min(tot_size, key=tot_size.get)
         return self.OPTIMAL_DELTA_SIZE

    def compress(self, geometry):

        # Create pre computed values to store as metadata
        s = time.perf_counter()
        optimal_size = self.calculate_delta_size(geometry)
        if optimal_size % 8 != 0:
            optimal_size += 8
        optimal_size /= 8
        
        res, _ = self.fp_delta_encoding(
            geometry, 2, self.MAX_CHUNK_SIZE)
        t = time.perf_counter()
        return t - s, res

    def decompress(self, bin):
        s = time.perf_counter()
        geometry = self.fp_delta_decoding(bin)
        t = time.perf_counter()
        return t - s, geometry


# ---- UNARY ---- #

    def vertices(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        coords = shapely.get_coordinates(geometry)
        t = time.perf_counter()
        return t - s, coords

    def type(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        type = geometry.geom_type
        t = time.perf_counter()

    def bounding_box(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds

    def add_vertex(self, args):
        res = None
        bin, idx, pos = args
        s = time.perf_counter()

        return t - s, res

    def is_intersecting(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        _, r_geo = self.decompress(r_bin)
        res = shapely.intersects(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, res

    def intersection(self, args):
        l_bin, r_bin = args
        s = time.perf_counter()
        _, l_geo = self.decompress(l_bin)
        _, r_geo = self.decompress(r_bin)
        res = shapely.intersection(l_geo, r_geo)
        t = time.perf_counter()
        return t - s, res


def main():
    x = Fpd()
    geom1 = shapely.wkt.loads("MULTIPOLYGON (((13.1910341 55.7064716, 13.1909909 55.7064777, 13.1903007 55.7065758, 13.1902852 55.7065781, 13.1903071 55.7066218, 13.190342 55.7067242, 13.190363 55.7067227, 13.1903945 55.7067205, 13.1903987 55.7067157, 13.1904352 55.7067094, 13.1910138 55.7066293, 13.1909872 55.7065641, 13.1910247 55.7065596, 13.1910719 55.7065528, 13.1910341 55.7064716)), ((13.1913015 55.7064068, 13.1912723 55.706411, 13.1912616 55.7064183, 13.191259 55.7064293, 13.191248 55.7064396, 13.1911833 55.7064522, 13.1911141 55.7064603, 13.1911564 55.7065469, 13.1911868 55.7065431, 13.1912061 55.7065406, 13.1912298 55.7066002, 13.1912449 55.7065981, 13.1916827 55.7065379, 13.1916619 55.7064795, 13.1916589 55.7064711, 13.1916541 55.7064574, 13.1916482 55.7064411, 13.1913299 55.7064803, 13.1913015 55.7064068)))")
    t, bin3 = x.compress(geom1)
    decomp = x.decompress(bin3)

    #print(x.calculate_delta_size(geom1))
    #t, bin = x.decompress(bin)


if __name__ == "__main__":
    main()
