from abc import ABC, abstractmethod


class CompressionAlgorithm(ABC):

    @abstractmethod
    def compress(self, geometry):
        pass

    @abstractmethod
    def decompress(self, bin):
        pass

# ---- UNARY ---- #
    @abstractmethod
    def vertices(self, bin):
        pass

    @abstractmethod  # Geometry type
    def type(self, bin): 
        pass

    @abstractmethod
    def bounding_box(self, bin):
        pass

    @abstractmethod  # For Polygon
    def add_vertex(self, args):
        pass

    @abstractmethod
    def add_vertex(self, args):
        pass

# ---- BINARY ---- #
    @abstractmethod  # For Polygon
    def is_intersecting(self, args):
        pass

    @abstractmethod
    def intersection(self, args):
        pass