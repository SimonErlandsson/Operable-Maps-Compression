
from base import CompressionAlgorithm
from collections import deque
import shapely
import shapely.wkt
import time
import struct
import sys
import math
class Fpd(CompressionAlgorithm):
    MAX_CHUNK_SIZE = 10
    OPTIMAL_DELTA_SIZE = 4 #DEFAULT VALUE
    type_int_convertion = {"Point": 0, "LineString": 1, "Polygon": 2, "MultiPolygon": 3, "MultiLineString": 4, "MultiPoint": 5, "LinearRing": 6, 
                  0: "Point", 1: "LineString", 2: "Polygon", 3: "MultiPolygon", 4: "MultiLineString", 5: "MultiPoint", 6: "LinearRing"} 

#---- HELPER METHODS
    def binary(self, num):
        return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))
    
    def int_repr(self, num):
        return struct.unpack('!i',struct.pack('!f',num))[0]
    
    def int64_to_int16(self, num):
        return int.from_bytes(struct.pack('!h', num), "big")
    
    def float_to_bytes(self, x):
        return struct.pack("f", x)
    
    def int_to_bytes(self, x):
        return struct.pack('>I', x)

    def zz_encode(self,num):
        return -2 * num - 1 if num < 0 else 2 * num
        
    def get_zz_encoded_delta(self, prev_coord, curr_coord):
        return self.zz_encode(self.int_repr(curr_coord) - self.int_repr(prev_coord))
    
    def deltas_fit_in_bytes(self, max_bytes, delta_x, delta_y):
        return math.log2(delta_x) <= max_bytes * 8 and math.log2(delta_y) <= max_bytes * 8

    def get_poly_ring_count(self, geometry):
        subpoly_coord_count = deque([])
        if geometry.geom_type == "Polygon" or geometry.geom_type == "MultiPolygon":
            subpoly_coord_count.append(len(geometry.exterior.coords))
            for i in range(len(geometry.interiors)):
                subpoly_coord_count.append(len(geometry.interiors[i].coords))
        return subpoly_coord_count
    
    def init_chunk(self, res, chk_size, coords, pnt_idx, chk_ys, chk_xs):
        chk_size_idx = len(res)
        res.append(self.int_to_bytes(chk_size))

        #Add full coordinates
        chk_xs.append(self.float_to_bytes(coords[pnt_idx][0]))
        chk_ys.append(self.float_to_bytes(coords[pnt_idx][1]))
        return chk_size_idx


    def save_clear_state(self,res, chk_xs, chk_ys):
        res += chk_xs + chk_ys
        #Reset current chk state
        chk_xs.clear()
        chk_ys.clear()

    def append_delta_pair(self, chk_xs, chk_ys, deltax, deltay, delta_size):
        delta_bytes_x = deltax.to_bytes(delta_size, 'big')
        delta_bytes_y = deltay.to_bytes(delta_size, 'big')
        chk_xs.append(delta_bytes_x)
        chk_ys.append(delta_bytes_y)


    def init_bytearray(self, geometry, delta_size, chk_size):
        res = []
        total_coords = len(shapely.get_coordinates(geometry))  
        # Operation metadata
        res.append(self.int_to_bytes(total_coords))
        #Add bounding box
        for bound in geometry.bounds:
            res.append(self.float_to_bytes(bound))

        #Meta data
        res.append(self.int_to_bytes(delta_size))
        res.append(self.int_to_bytes(chk_size))
        res.append(self.type_int_convertion[geometry.geom_type].to_bytes(1, 'big')) #1 byte is enought for storing type
        return res


    def fp_delta_encoding(self, geometry, max_chk_size, delta_size):
        coords = shapely.get_coordinates(geometry)

        #Saves interior and exterior lengths if polygon type
        subpoly_coord_count:deque = self.get_poly_ring_count(geometry)
    
        #State of the current chk state
        chk_xs, chk_ys = [], []
        chk_size_idx = 0  #Used for changing the chk_size of reset occurs

        #State of the current coordiate pair
        chk_point_nbr = 0
        
        #Array for storing list resulting bytes. Initialized with Total coordinates, bounding box, delta_size, chk_size, geom_type
        res = self.init_bytearray(geometry,  delta_size, chk_size)

        #Polygon specific variables        
        is_polygon = len(subpoly_coord_count) != 0
        subpoly_point_cnt = 0

        #Loop all coordinates
        for i in range(len(coords)):
            if is_polygon and subpoly_point_cnt == 0:
                #Reset and store previous chk data if interrupted in the middle of a chk
                if(len(chk_xs) != 0):
                    res[chk_size_idx] = self.int_to_bytes(chk_point_nbr)
                    chk_point_nbr = 0
                    self.save_clear_state(res, chk_xs, chk_ys)


                subpoly_point_cnt = subpoly_coord_count.popleft()
                res.append(self.int_to_bytes(subpoly_point_cnt))


            #----CHUNK HEADER INFORMATION (CHUNK SIZE, FULL FIRST COORDINATES)
            if chk_point_nbr == 0:
                #save chk size and later change if reset occurs
                chk_size_idx = self.init_chunk(res, chk_size, coords, i, chk_ys, chk_xs)
                chk_point_nbr += 1

            else: #Loop for delta
                zig_delta_x = self.get_zz_encoded_delta(coords[i - 1][0],coords[i][0])
                zig_delta_y = self.get_zz_encoded_delta(coords[i - 1][1],coords[i][1])

                if self.deltas_fit_in_bytes(delta_size,zig_delta_x,zig_delta_y):
                    self.append_delta_pair(chk_xs, chk_ys, zig_delta_x, zig_delta_y, delta_size)
                    chk_point_nbr += 1
                    
                else:
                    res[chk_size_idx] = self.int_to_bytes(chk_point_nbr)
                    self.save_clear_state(res, chk_xs, chk_ys)
                    chk_size_idx = self.init_chunk(res, chk_size, coords, i, chk_ys, chk_xs)

                    chk_point_nbr = 1

                if chk_point_nbr == max_chk_size:
                    chk_point_nbr = 0
                    self.save_clear_state(res, chk_xs, chk_ys)

            subpoly_point_cnt -= 1

        if chk_point_nbr > 0:
            self.save_clear_state(res, chk_xs, chk_ys)
            res[chk_size_idx] = self.int_to_bytes(chk_point_nbr)

        res = b''.join(res)
        return res, len(res)
    

    def fp_delta_decoding(self, bin):
        bin = bin[4 + 4 * 2:] #Remove header information about total nbr of nodes and bounding box
        delta_size = int.from_bytes(bin[0:4],byteorder='big')
        print(delta_size)
        reset_next = True
        chk_deltas = 0

        return bin, None



    def calculate_delta_sizes(self, geometry):
        _, variable_size = self.fp_delta_encoding(geometry, self.MAX_CHUNK_SIZE, delta_size=None)
        
        #    res_sizes.append(variable_size)
        #for i in range(len(delta_sizes)):
        #    print("delta size of", delta_sizes[i], "bytes -> total size of geometry:", res_sizes[i], "bytes")
        #self.OPTIMAL_DELTA_SIZE = delta_sizes[res_sizes.index(min(res_sizes))]
        return self.OPTIMAL_DELTA_SIZE



    def compress(self, geometry):

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        optimal_size = self.getOptimalDeltaSize(geometry)
        res, _ = self.fp_delta_encoding(geometry, optimal_size, self.MAX_CHUNK_SIZE)
        t = time.perf_counter()
        return t - s, res


    def decompress(self, bin):
        res = None
        s = time.perf_counter()
        res, _ = self.fp_delta_decoding(bin, self.OPTIMAL_DELTA_SIZE, self.MAX_CHUNK_SIZE)
        t = time.perf_counter()
        return t - s, res     


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
    geom1 = shapely.wkt.loads("POLYGON ((13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362))")
   
    t, bin3 = x.compress(geom1)
    print("___")
    print(bin3)

    #t, bin = x.decompress(bin)   




if __name__ == "__main__":
    main()
