
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
    offset = 0 # Used when parsing

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
        return delta_size, type, 5 # Offset is 5 bytes for IB

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


    def bytes_to_float32(self, bin):
        val = struct.unpack_from('!f', bin, self.offset)[0] # Nbr of deltas
        self.offset += 4
        return val
    
    def bytes_to_float(self, bin, size=32):
        # For now only bytes:
        size //= 8
        bin = bin[self.offset : self.offset + size]
        print(bin)
        self.offset += size
        return bin

    def fp_delta_decoding(self, bin):
        geometry = shapely.LineString()
        self.offset = 0
        delta_size, type, offset = self.decode_header(bin)
        
        if type == GT.LINESTRING:
            while (self.offset < len(bin)): # While != EOF
                chk_size = struct.unpack_from('!I', bin, offset)[0] 
                self.offset += 4
                coords = []

                # Extract reset point
                x = self.bytes_to_float32(bin)
                y = self.bytes_to_float32(bin)
                coords.append((x, y))

                # Loop through deltas in chunk
                for i in range(chk_size):
                    d_x = self.bytes_to_float(bin, size=delta_size * 8)
                    d_y = self.bytes_to_float(bin, size=delta_size * 8)
                    x += d_x # Add delta to old x value
                    y += d_y
                    coords.append((x, y))
                #print(coords)

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
    geom1 = shapely.wkt.loads("MULTIPOLYGON (((13.193709 55.7021381, 13.1937743 55.7021279, 13.1938355 55.7021184, 13.1938461 55.702109, 13.1938566 55.7020984, 13.1938611 55.7020902, 13.1938655 55.7020774, 13.1938655 55.7020633, 13.1938583 55.7020408, 13.1938402 55.7020014, 13.1937184 55.7017259, 13.1937008 55.7016876, 13.1936836 55.7016654, 13.1936537 55.7016428, 13.1936223 55.7016242, 13.1935741 55.7016036, 13.1935354 55.7015911, 13.1935006 55.701584, 13.1934829 55.701598, 13.1934673 55.7016115, 13.1934736 55.7016164, 13.1934776 55.7016216, 13.1934875 55.7016633, 13.1934985 55.7016898, 13.1935196 55.7017337, 13.1935659 55.7018353, 13.1936162 55.7018282, 13.1936551 55.7019155, 13.1936651 55.7019377, 13.1936955 55.7020047, 13.1936497 55.7020119, 13.193709 55.7021381)), ((13.1938175 55.7017126, 13.1938602 55.7017068, 13.1939048 55.7017007, 13.1938998 55.7016861, 13.193892 55.7016685, 13.1938831 55.7016589, 13.193871 55.701651, 13.1938602 55.701646, 13.1938405 55.7016438, 13.193822 55.7016456, 13.1938062 55.7016517, 13.1937985 55.7016571, 13.1937953 55.7016646, 13.1937979 55.7016746, 13.1938017 55.7016836, 13.1938052 55.7016908, 13.1938175 55.7017126)), ((13.1940245 55.7019788, 13.19398 55.7019848, 13.1939372 55.7019907, 13.1939585 55.7020383, 13.1939692 55.7020479, 13.1939841 55.7020512, 13.1939975 55.7020519, 13.1940079 55.702051, 13.1940198 55.7020497, 13.1940317 55.7020463, 13.1940395 55.7020422, 13.1940435 55.7020369, 13.1940452 55.7020314, 13.1940457 55.7020218, 13.1940245 55.7019788)), ((13.1939779 55.7015541, 13.1939529 55.701555, 13.1939622 55.7015658, 13.1939755 55.7015942, 13.194075 55.7018201, 13.1941382 55.7019637, 13.1941483 55.7019866, 13.194164 55.7020087, 13.1941899 55.7020304, 13.1942142 55.7020424, 13.1942291 55.7020486, 13.1942638 55.702042, 13.195019 55.7018988, 13.1948681 55.7018923, 13.1944181 55.7018687, 13.1944172 55.7018717, 13.194395 55.7018706, 13.1942164 55.7018622, 13.194172 55.7017564, 13.1941218 55.701761, 13.1941279 55.7017262, 13.1941357 55.7016818, 13.1940872 55.7015737, 13.1940769 55.7015503, 13.1939779 55.7015541), (13.1942341 55.7020059, 13.1942075 55.7020095, 13.1941895 55.7019673, 13.1941696 55.701921, 13.1941936 55.7019177, 13.1941884 55.7019055, 13.19426 55.7018958, 13.1942645 55.7019063, 13.1943172 55.7018991, 13.1943567 55.7019912, 13.1943041 55.7019984, 13.1943086 55.7020089, 13.1942394 55.7020183, 13.1942341 55.7020059)))")
    geom2 = shapely.wkt.loads("LINESTRING (13.199378 55.7034667, 13.1999441 55.7033986, 13.200125 55.7033882, 13.2002723 55.7033936, 13.2004383 55.7034097, 13.2005935 55.7034211, 13.2007699 55.703423, 13.2011275 55.7034136, 13.2012413 55.7034103, 13.2012947 55.7034088)")
    #t, bin3 = x.compress(geom2)
    #print(bin3)
    #decomp = x.decompress(bin3)
    print(bytes(b'\xcf\x84\xcf\x81\xce\xbd\xcf\x82').binary())
    #print(x.calculate_delta_size(geom1))
    #t, bin = x.decompress(bin)


if __name__ == "__main__":
    main()
