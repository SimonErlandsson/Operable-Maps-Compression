
from base import CompressionAlgorithm
import shapely
import shapely.wkt
import time

class Fpd(CompressionAlgorithm):
    CHUNK_SIZE = 10
    def compress(self, geometry):

        #Create pre computed values to store as metadata
        s = time.perf_counter()
        
        chunk_nbr, delta_nbr = 0 #State of the current coordiate
        
        #Information about the geometry in general
        coords = shapely.get_coordinates(geometry)
        total_coords = len(coords)
        allighned_chunks = total_coords // self.CHUNK_SIZE
        deltas_on_last_chunk = total_coords % self.CHUNK_SIZE

        total_chunks = allighned_chunks        
        if deltas_on_last_chunk != 0:
            total_chunks += 1

        #Encoding loop
        res = []
        for coords in shapely.get_coordinates(geometry):

            #Iteration for first chunk header for first block
            if chunk_nbr == 0 and delta_nbr == 0:
                #Append number of chunks in geometry
                res.append(len(coords))
               
                #Add bounding box
                for bound in geometry.bounds: #Always two coordinates to O(1)
                    res.append(bound)

            #For all chunks
            if delta_nbr == 0:
                if chunk_nbr <= allighned_chunks:
                    res.append(self.CHUNK_SIZE)
                else:
                    res.append(deltas_on_last_chunk)


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

if __name__ == "__main__":
    main()
