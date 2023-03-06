
from base import CompressionAlgorithm
from pympler import asizeof
import shapely
import shapely.wkt
import time
import struct
import sys
import math
class Fpd(CompressionAlgorithm):

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

    def fp_delta_encoding(self, geometry, delta_size, chunk_size):
        coords = shapely.get_coordinates(geometry)
        total_coords = len(coords)

        #State of the current coordiate pair
        delta_nbr = 0 
        
        #State of the current chunk state
        chunk_xs = []
        chunk_ys = []
        chunk_size_idx = 0  #Used for changing the chunk_size of reset occurs


        #Array for storing list of full coords, deltas, and metadata.
        res = []


        #----FIRST HEADER INFORMATION(TOTAL_SIZE + BOUNDING BOX)----
        res.append(self.int_to_bytes(total_coords))

        #Add bounding box
        for bound in geometry.bounds:
            res.append(self.float_to_bytes(bound))

        #Loop all coordinates
        for i in range(total_coords):
            
            #----CHUNK HEADER INFORMATION (CHUNK SIZE, FULL FIRST COORDINATES)
            if delta_nbr == 0:
                #save chunk size and later change if reset occurs
                chunk_size_idx = len(res)
                res.append(self.int_to_bytes(chunk_size))

                #Add full coordinates
                chunk_xs.append(self.float_to_bytes(coords[i][0]))
                chunk_ys.append(self.float_to_bytes(coords[i][1]))

                delta_nbr += 1

            else: #Loop for delta
                zig_delta_x = self.get_zz_encoded_delta(coords[i - 1][0],coords[i][0])
                zig_delta_y = self.get_zz_encoded_delta(coords[i - 1][1],coords[i][1])

                if self.deltas_fit_in_bytes(delta_size,zig_delta_x,zig_delta_y):
                    delta_bytes_x = zig_delta_x.to_bytes(delta_size, 'big')
                    delta_bytes_y = zig_delta_y.to_bytes(delta_size, 'big')
                    chunk_xs.append(delta_bytes_x)
                    chunk_ys.append(delta_bytes_y)
                    delta_nbr += 1
                     
                else:

                    res[chunk_size_idx] = self.int_to_bytes(delta_nbr)
                    
                    res += chunk_xs
                    res += chunk_ys
                    #Reset current chunk state
                    chunk_xs.clear()
                    chunk_ys.clear()
                    
                    chunk_size_idx = len(res)
                    res.append(self.int_to_bytes(chunk_size))
                    #Add full coordinates
                    chunk_xs.append(self.float_to_bytes(coords[i][0]))
                    chunk_ys.append(self.float_to_bytes(coords[i][1]))

                    delta_nbr = 1

                if delta_nbr == chunk_size:
                    delta_nbr = 0
                    res += chunk_xs
                    res += chunk_ys

                    #Reset current chunk state
                    chunk_xs.clear()
                    chunk_ys.clear()

        if delta_nbr > 0:
            res += chunk_xs
            res += chunk_ys
            res[chunk_size_idx] = self.int_to_bytes(delta_nbr)
        res = b''.join(res)
        return res, len(res)




    def getMinimumDeltaSize(self, geometry):
        delta_sizes  = [1,2,3,4,5,6,7,8]
        res_sizes = []
        for delta_size in delta_sizes:
            _, variable_size = self.fp_delta_encoding(geometry, delta_size, 10)
            res_sizes.append(variable_size)
        return delta_sizes[res_sizes.index(min(res_sizes))]

        



    CHUNK_SIZE = 10
    def compress(self, geometry):

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        optimal_size = self.getMinimumDeltaSize(geometry)
        res, _ = self.fp_delta_encoding(geometry, optimal_size, 10)
        t = time.perf_counter()
        return t - s, res


    def decompress(self, bin):
        res = None
        s = time.perf_counter()
        
        t = time.perf_counter()
        return t - s, res     


# ---- UNARY ---- #

    def vertices(self, bin):
        res = None
        s = time.perf_counter()

        t = time.perf_counter()
        return t - s, res


    def type(self, bin): 
        res = None
        s = time.perf_counter()

        return t - s, res 
    
    def bounding_box(self, bin):
        res = None
        s = time.perf_counter()

        return t - s, res
    
    def add_vertex(self, args):
        res = None        
        bin, idx, pos = args
        s = time.perf_counter()
  
        return t - s, res
    
    def is_intersecting(self, args):
        res = None        
        l_bin, r_bin = args
        s = time.perf_counter()

        t = time.perf_counter()
        return t - s, res

    def intersection(self, args):
        res = None        
        l_bin, r_bin = args
        s = time.perf_counter()

        t = time.perf_counter()
        return t - s, res
    
def main():
    x = Fpd()
    geom = shapely.wkt.loads("POLYGON ((13.1635138 55.7053599, 13.1637569 55.7053536, 13.1635571 55.705336, 13.1635158 55.7053284, 13.1635184 55.7053437, 13.1635138 55.7053599), (13.1667021 55.7046362, 13.1667117 55.7046498, 13.1667021 55.7046362))")
    t, bin = x.compress(geom)
    print(bin)
   




if __name__ == "__main__":
    main()
