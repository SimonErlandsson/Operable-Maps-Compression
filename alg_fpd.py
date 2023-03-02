
from base import CompressionAlgorithm
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
    

    def zz_encode(self,num):
        return -2 * num - 1 if num < 0 else 2 * num
        

    def fp_delta_encoding(self, geometry, delta_size, chunk_size):
        coords = shapely.get_coordinates(geometry)
        variable_size = 0

        #State of the current coordiate pair
        delta_nbr = 0 

        total_coords = len(coords)

        #Array for storing list of full coords, deltas, and metadata.
        res = []

        #Arrays for storing current chunk state
        current_chunk_x = []
        current_chunk_y = []
        chunk_size_idx = 0

        res.append(total_coords)
        #Add bounding box
        for bound in geometry.bounds: #Always two coordinates to O(1)
            res.append(bound)

        #Loop all coordinates
        for i in range(total_coords):

            #For all chunks beginnings
            if delta_nbr == 0:
                #save chunk size and later change if reset occurs
                chunk_size_idx = len(res)
                res.append(chunk_size)
                variable_size += 32 #Assuming int32 chunksize

                #Reset current chunk state
                current_chunk_x.clear()
                current_chunk_y.clear()

                #Add full coordinates
                current_chunk_x.append(coords[i][0])
                current_chunk_y.append(coords[i][1])
                variable_size += 32 * 2 #Assuming float32 bits
                delta_nbr += 1

            else: #Loop for delta
                delta_x = self.int_repr(coords[i][0]) - self.int_repr(coords[i - 1][0])
                delta_y = self.int_repr(coords[i][1]) - self.int_repr(coords[i - 1][1])
              
                zig_delta_x:int = self.zz_encode(delta_x)
                zig_delta_y:int = self.zz_encode(delta_y)
                
                if math.log2(zig_delta_x) <= delta_size and math.log2(zig_delta_y) <= delta_size:
                    delta_bytes_x = zig_delta_x.to_bytes(delta_size, 'big')
                    delta_bytes_y = zig_delta_y.to_bytes(delta_size, 'big')
                    current_chunk_x.append(delta_bytes_x)
                    current_chunk_y.append(delta_bytes_y)
                    delta_nbr += 1
                    variable_size += delta_size * 2 * 8
                     
                else:
                    res[chunk_size_idx] = delta_nbr
                    chunk_size_idx = len(res)
                    res.append(chunk_size)
                    variable_size += 32 #Assuming int32 chunksize
                    res += current_chunk_x
                    res += current_chunk_y

                    #Reset current chunk state
                    current_chunk_x.clear()
                    current_chunk_y.clear()

                    #Add full coordinates
                    current_chunk_x.append(coords[i][0])
                    current_chunk_y.append(coords[i][1])
                    variable_size += 32 * 2 #Assuming float32 bits
                    delta_nbr = 1

                if delta_nbr == self.CHUNK_SIZE:
                    delta_nbr = 0
                    res += current_chunk_x
                    res += current_chunk_y

        return res, variable_size




    def getMinimumDeltaSize(self, geometry):
        delta_sizes  = [2,3,4,5,6]
        res_sizes = []
        for delta_size in delta_sizes:
            _, variable_size = self.fp_delta_encoding(geometry, delta_size, 10)
            res_sizes.append(variable_size)
        print(res_sizes)
        return delta_sizes[res_sizes.index(min(res_sizes))]

        



    CHUNK_SIZE = 10
    def compress(self, geometry):

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        optimal_size = self.getMinimumDeltaSize(geometry)
        print(optimal_size)
        res, _ = self.fp_delta_encoding(geometry, 4, 10)
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
