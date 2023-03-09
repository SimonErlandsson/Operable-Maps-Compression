
from base import CompressionAlgorithm
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import sys
import math
from shapely import GeometryType as GT

## char: 8 bits
## float: 32 bits
## double: 64 bits
## int: 32 bits
## long: 64 bits


class Fpd(CompressionAlgorithm):
    MAX_CHUNK_SIZE = 2
    offset = 0 # Used when parsing

# ---- HELPER METHODS
    def double_as_long(self, num):
        return struct.unpack('!q', struct.pack('!d', num))[0]
    
    def long_as_double(self, num):
        return struct.unpack('!d', struct.pack('!q', num))[0]

    def double_to_bytes(self, x):
        return struct.pack("!d", x)

    def uint_to_bytes(self, x):
        return struct.pack('!I', x)

    def zz_encode(self, num):
        return -2 * num - 1 if num < 0 else 2 * num
    
    def zz_decode(self, num):
        return - num // 2 if num & 1 == 1 else num // 2

    def get_zz_encoded_delta(self, prev_coord, curr_coord):
        return self.zz_encode(self.double_as_long(curr_coord) - self.double_as_long(prev_coord))

    def deltas_fit_in_bits(self, max_bits, delta_x, delta_y):
        return (delta_x == 0 or math.log2(delta_x) <= max_bits) and (delta_y == 0 or math.log2(delta_y) <= max_bits)

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
            poly_coord_count = deque([len(list(polygon.interiors)) + 1 for polygon in list(geometry.geoms)])
        return poly_coord_count
    
    def init_chunk(self, res, chk_size, coords, pnt_idx):
        chk_size_idx = len(res)
        res.append(self.uint_to_bytes(chk_size))

        # Add full coordinates
        res.append(self.double_to_bytes(coords[pnt_idx][0]))
        res.append(self.double_to_bytes(coords[pnt_idx][1]))
        return chk_size_idx

    def append_delta_pair(self, res, d_x, d_y, delta_size):
        delta_bytes_x = d_x.to_bytes(delta_size // 8, 'big')
        delta_bytes_y = d_y.to_bytes(delta_size // 8, 'big')
        res.extend([delta_bytes_x, delta_bytes_y])


    def init_bytearray(self, geometry, delta_size):
        res = []
        # Meta data
        res.append(self.uint_to_bytes(delta_size))
        res.append(int(shapely.get_type_id(geometry)).to_bytes(1, 'big')) #1 byte is enough for storing type

        return res
    
    def decode_header(self, bin):
        delta_size, type = struct.unpack_from('!IB', bin)
        type = GT(type)
        self.offset += 5
        return delta_size, type # Offset is 5 bytes for IB

    def polygon_reset(self, res, geom_coord_count, chk_point_cnt, chk_size_idx):
        #Reset and store previous chk data if interrupted in the middle of a chk
        if chk_size_idx != 0:
            res[chk_size_idx] = self.uint_to_bytes(chk_point_cnt - 1)
        subpoly_point_cnt = geom_coord_count.popleft()
        res.append(self.uint_to_bytes(subpoly_point_cnt))
        chk_point_cnt = 0
        return chk_point_cnt, subpoly_point_cnt

        
    def fp_delta_encoding(self, geometry, delta_size):
        coords = shapely.get_coordinates(geometry)

        #Saves interior and exterior lengths if polygon type
        subpoly_coord_count:deque = self.get_poly_ring_count(geometry)
        poly_coord_count:deque = self.get_multipoly_count(geometry)

        chk_size_idx = 0  #Used for changing the chk_size of reset occurs

        #State of the current coordiate pair
        chk_point_cnt = 0
        
        #Array for storing list resulting bytes. Initialized with Total coordinates, bounding box, delta_size, chk_size, geom_type
        res = self.init_bytearray(geometry,  delta_size)

        #Polygon specific variables        
        is_multipolygon = len(poly_coord_count) != 0
        is_polygon = len(subpoly_coord_count) != 0
        poly_ring_cnt, subpoly_point_cnt = 0, 0

        #Loop all coordinates
        for i in range(len(coords)):
            if is_multipolygon and poly_ring_cnt == 0:
               chk_point_cnt, poly_ring_cnt = self.polygon_reset(res, poly_coord_count, chk_point_cnt, chk_size_idx)

            if is_polygon and subpoly_point_cnt == 0:
               chk_point_cnt, subpoly_point_cnt = self.polygon_reset(res, subpoly_coord_count, chk_point_cnt, chk_size_idx)

            #----CHUNK HEADER INFORMATION (CHUNK SIZE, FULL FIRST COORDINATES)
            if chk_point_cnt == 0:
                #save chk size and later change if reset occurs
                chk_size_idx = self.init_chunk(res, self.MAX_CHUNK_SIZE, coords, i)
                chk_point_cnt += 1

            else: #Loop for delta
                zig_delta_x = self.get_zz_encoded_delta(coords[i - 1][0],coords[i][0])
                zig_delta_y = self.get_zz_encoded_delta(coords[i - 1][1],coords[i][1])

                if self.deltas_fit_in_bits(delta_size,zig_delta_x,zig_delta_y):
                    self.append_delta_pair(res, zig_delta_x, zig_delta_y, delta_size)
                    chk_point_cnt += 1
                    
                else:
                    res[chk_size_idx] = self.uint_to_bytes(chk_point_cnt - 1)
                    chk_size_idx = self.init_chunk(res, self.MAX_CHUNK_SIZE, coords, i)
                    chk_point_cnt = 1

                if chk_point_cnt == self.MAX_CHUNK_SIZE:
                    chk_point_cnt = 0

            subpoly_point_cnt -= 1
       
            if subpoly_point_cnt == 0:
                poly_ring_cnt -= 1

        if chk_point_cnt > 0:
            res[chk_size_idx] = self.uint_to_bytes(chk_point_cnt - 1)
        print([int.from_bytes(i, 'big') for i in res])
        res = b''.join(res)
        return res, len(res)

    
    def bytes_to_decoded_coord(self, bin, prev_coord, input_size=64):
        # For now only bytes:
        size = input_size // 8
        bin = bin[self.offset : self.offset + size]
        val = int.from_bytes(bin, 'big', signed=False)
        val = self.zz_decode(val) + self.double_as_long(prev_coord)
        val = self.long_as_double(val)
        self.offset += size
        return val
    
    def bytes_to_double(self, bin):
        val = struct.unpack_from('!d', bin, self.offset)[0]
        self.offset += 8
        return val
    
    def bytes_to_uint(self, bin):
        val = struct.unpack_from('!I', bin, self.offset)[0]
        self.offset += 4
        return val

    def fp_delta_decoding(self, bin):
        self.offset = 0
        delta_size, type = self.decode_header(bin)
        if type == GT.LINESTRING:
            coords = [] 
            while (self.offset < len(bin)): # While != EOF                
                chk_size = self.bytes_to_uint(bin)
                # Extract reset point
                x = self.bytes_to_double(bin)
                y = self.bytes_to_double(bin)
                coords.append((x, y))

                # Loop through deltas in chunk
                for i in range(chk_size):
                    x = self.bytes_to_decoded_coord(bin, x, delta_size)
                    y = self.bytes_to_decoded_coord(bin, y, delta_size)
                    coords.append((x, y))
            
            # All coords added
            geometry = shapely.LineString(coords)
        elif type == GT.POLYGON:
            coords = []
            print(delta_size)
            while (self.offset < len(bin)): # While != EOF
                chks_in_ring = self.bytes_to_uint(bin) 
                print(chks_in_ring)   

                for i in range(chks_in_ring):   
                    chk_size = self.bytes_to_uint(bin)
                    # Extract reset point
                    x = self.bytes_to_double(bin)
                    y = self.bytes_to_double(bin)
                    coords.append((x, y))
                    print(chk_size)
                    # Loop through deltas in chunk
                    for i in range(chk_size):
                        #print(i)
                        x = self.bytes_to_decoded_coord(bin, x, delta_size)
                        y = self.bytes_to_decoded_coord(bin, y, delta_size)
                        coords.append((x, y))
            
            # All coords added
            geometry = shapely.LineString(coords)
        elif type == GT.MULTIPOLYGON:
            pass



        return geometry



    def calculate_delta_size(self, geometry):
         RESET_POINT_SIZE = 64 + 32
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
             tot_size[n] = n * lower_cnt + RESET_POINT_SIZE * upper_cnt
             lower_cnt -= bit_cnts[n]
             upper_cnt += bit_cnts[n]

         #    res_sizes.append(variable_size)
         # for i in range(len(delta_sizes)):
         #    print("delta size of", delta_sizes[i], "bytes -> total size of geometry:", res_sizes[i], "bytes")
         return min(tot_size, key=tot_size.get)

    def compress(self, geometry):

        # Create pre computed values to store as metadata
        s = time.perf_counter()
        optimal_size = self.calculate_delta_size(geometry)
        optimal_size = math.ceil(optimal_size / 8) * 8
        
        res, _ = self.fp_delta_encoding(
            geometry, optimal_size)
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
        return t - s, type

    def bounding_box(self, bin):
        s = time.perf_counter()
        _, geometry = self.decompress(bin)
        bounds = shapely.bounds(geometry)
        t = time.perf_counter()
        return t - s, bounds

    def add_vertex(self, args):
        bin, idx, pos = args
        s = time.perf_counter()

        _, geometry = self.decompress(bin)
        ragged = shapely.to_ragged_array([geometry])
        points = np.concatenate((ragged[1][:idx, :], np.array([pos]), ragged[1][idx:, :]), axis=0) #Appearentlyq

        offsets = list(ragged[2])
        offsets[0] = np.array(list(map(lambda x: x + 1 if x > idx else x, offsets[0])))

        res = shapely.from_ragged_array(geometry_type=shapely.get_type_id(geometry), coords=points, offsets=offsets)[0]
        t = time.perf_counter()
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
    geom3 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848603 55.7056619, 13.1846238 55.7056422, 13.1846085 55.7057159, 13.1846356 55.7057179, 13.1848537 55.7057363), (13.1846694 55.705714, 13.1846543 55.7057128, 13.1846563 55.705705, 13.1846714 55.7057062, 13.1846694 55.705714), (13.1847425 55.7057123, 13.1847405 55.7057201, 13.1847254 55.7057188, 13.1847274 55.705711, 13.1847425 55.7057123), (13.1848001 55.7057179, 13.1848152 55.7057192, 13.1848131 55.705727, 13.1847981 55.7057258, 13.1848001 55.7057179), (13.1848068 55.7056929, 13.1848088 55.7056851, 13.1848239 55.7056863, 13.1848218 55.7056941, 13.1848068 55.7056929), (13.1847507 55.7056878, 13.1847356 55.7056865, 13.1847377 55.7056787, 13.1847528 55.70568, 13.1847507 55.7056878), (13.1846811 55.7056732, 13.184679 55.705681, 13.184664 55.7056798, 13.184666 55.705672, 13.1846811 55.7056732))')
    geom4 = shapely.wkt.loads('POLYGON ((13.1848537 55.7057363, 13.1848861 55.705646, 13.1848537 55.7057363), (13.1847425 55.7057123, 13.1847405 55.7057201, 13.1847254 55.7057188, 13.1847274 55.705711, 13.1847425 55.7057123))')
    t, bin3 = x.compress(geom4)
    #print(x.bytes_to_float32(b'0xa70x300x530x41'))
    print(bin3.hex(sep=' '))
    #print("-DELTASIZE-_TY_POINTSINCHK_XFIRST-----------------_YFIRST-----------------")
    print("-DELTASIZE-_TY_RINGS------_POINTSINCHK_XFIRST-----------------_YFIRST-----------------_--XD1---------_--YD1---------_--XD2---------_--YD2---------_--XD3---------_--YD3---------_--XD4---------_--YD4---------_--XD5---------_--YD5---------_--XD6---------_--YD6---------_--XD7---------_--YD7---------_--XD8---------_--YD8---------_--XD9---------_--YD9---------")
    print(shapely.to_ragged_array([geom4]))
    t, decomp = x.decompress(bin3)
    #print(x.bytes_to_float(bin3, size=16))
    #print(x.calculate_delta_size(geom1))
    #t, bin = x.decompress(bin)

    #print(bin(0b01100111 | 0b0000000000000000))
    #t = '{0:0{1}b}'.format(0b101, 10)
    #print(t)
    #print(int(t, 2))

    print(len(bytes(shapely.to_wkt(geom2), 'utf-8')), shapely.to_wkt(geom2))
    print(len(bin3), shapely.to_wkt(decomp))

    


if __name__ == "__main__":
    main()
