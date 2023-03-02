
from base import CompressionAlgorithm
import shapely
import shapely.wkt
import time
import struct

class Fpd(CompressionAlgorithm):


    def binary(self, num):
        return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))
    
    
    def int_repr(self, num):
        return struct.unpack('!q', struct.pack('!f', num))[0]
    
    def int_to_float16(self, num):
        return int.from_bytes(struct.pack('!q', num), "big")
        




    CHUNK_SIZE = 10
    def compress(self, geometry):

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        
        #State of the current coordiate pair
        chunk_nbr= 0
        delta_nbr = 0 
        
        #Information about the geometry in general
        coords = shapely.get_coordinates(geometry)
        total_coords = len(coords)

        allighned_chunks = total_coords // self.CHUNK_SIZE
        deltas_on_last_chunk = total_coords % self.CHUNK_SIZE

        #Calcualte number of chunks depending on unaligned number of deltas
        total_chunks = allighned_chunks        
        if deltas_on_last_chunk > 0:
            total_chunks += 1

        #Array for storing list of full coords, deltas, and metadata.
        res = []

        #Arrays for storing current chunk state
        current_chunk_x = []
        current_chunk_y = []

        #Loop all coordinates
        for i in range(total_coords):

            #Iteration for first chunk header for first block
            if chunk_nbr == 0 and delta_nbr == 0:
                #Append number of chunks in geometry
                res.append(total_coords)
               
                #Add bounding box
                for bound in geometry.bounds: #Always two coordinates to O(1)
                    res.append(bound)

            #For all chunks
            if delta_nbr == 0:
                if chunk_nbr <= allighned_chunks:
                    res.append(self.CHUNK_SIZE)
                else:
                    res.append(deltas_on_last_chunk)

                #Reset current chunk state
                current_chunk_x.clear()
                current_chunk_y.clear()

                #Add full coordinates
                current_chunk_x.append(coords[i][0])
                current_chunk_y.append(coords[i][1])
                delta_nbr += 1

            else: #Loop for delta
                delta_nbr += 1
                delta_x = self.int_repr(coords[i][0]) - self.int_repr(coords[i - 1][0])
                delta_y = self.int_repr(coords[i][1]) - self.int_repr(coords[i - 1][1])
              
                zig_delta_x = self.zz_encode(delta_x)
                zig_delta_y = self.zz_encode(delta_y)

                print(zig_delta_x, zig_delta_y)
                current_chunk_x.append(self.int_to_float16(zig_delta_x))
                current_chunk_y.append(self.int_to_float16(zig_delta_y))
                print(zig_delta_x, self.int_to_float16(zig_delta_x))
                if delta_nbr >= self.CHUNK_SIZE:
                    delta_nbr = 0
                    res += current_chunk_x
                    res += current_chunk_y
            




        t = time.perf_counter()
        return t - s, res
    
    def zz_encode(self,num):
        return -2 * num - 1 if num < 0 else 2 * num
        
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
   




if __name__ == "__main__":
    main()
