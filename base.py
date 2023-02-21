from abc import ABC, abstractmethod


class CompressionAlgorithm(ABC):

    @abstractmethod
    def compress(self, file_uncomp, file_comp):
        pass

    @abstractmethod
    def decompress(self, file_comp, file_decomp):
        pass

# ---- UNARY ---- #
    @abstractmethod
    def vertices(self, idx):
        pass

    @abstractmethod  # Geometry type
    def type(self, idx):
        pass

    @abstractmethod
    def vertex_count(self, idx):
        pass

    @abstractmethod  # For Polygon
    def area(self, idx):
        pass

    @abstractmethod
    def length(self, idx):
        pass

    # @abstractmethod
    # def end_points(self, idx):  # For LineString
    #     pass

    # @abstractmethod
    # def add_vertex(self, idx, vertex):  # vertex: (INT_idx, (FLOAT_x, FLOAT_y))
    #     pass

    # @abstractmethod
    # def remove_vertex(self, idx, vertex_idx):
    #     pass

    # @abstractmethod
    # def move_vertex(self, idx, vertex):  # vertex: (INT_idx, (FLOAT_x, FLOAT_y))
    #     pass

    # @abstractmethod
    # def bounding_box(self, idx):
    #     pass

    # @abstractmethod
    # def scale(self, idx, factor):
    #     pass

    # @abstractmethod
    # def center_of_mass(self, idx):
    #     pass

    # @abstractmethod
    # def convex_hull(self, idx):
    #     pass

# ---- BINARY ---- #
# https://en.wikipedia.org/wiki/DE-9IM
# TODO: Add more as we advance

    # @abstractmethod
    # def equals(self, idx_l, idx_r):
    #     pass

    # @abstractmethod
    # def contains(self, idx_l, idx_r):
    #     pass

    # @abstractmethod
    # def is_intersecting(self, idx_l, idx_r):  # Predicate: returns boolean
    #     pass

    # @abstractmethod
    # def intersection(self, idx_l, idx_r):  # Returns the intersecting polygon
    #     pass

    # @abstractmethod
    # def union(self, idx_l, idx_r):  # Returns the union polygon
    #     pass
