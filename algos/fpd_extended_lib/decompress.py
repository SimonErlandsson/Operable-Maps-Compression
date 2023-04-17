from bitarray import bitarray, util
import time
import shapely
import struct
import algos.fpd_extended_lib.intersection_chunk_bbox_wrapper 
from algos.fpd_extended_lib.low_level import *
import algos.fpd_extended_lib.cfg as cfg
from algos.fpd_extended_lib.cfg import *
from shapely import GeometryType as GT


# Structural things (per type):
def sequence_decoder(bin, seq_list, delta_size):
    cfg.offset += D_CNT_SIZE
    chk_size = bytes_to_uint(bin, D_CNT_SIZE)
    # Extract reset point
    x = bytes_to_double(bin)
    y = bytes_to_double(bin)
    seq_list.append((x, y))
    # Loop through deltas in chunk
    for _ in range(chk_size):
        x = bytes_to_decoded_coord(bin, x, delta_size)
        y = bytes_to_decoded_coord(bin, y, delta_size)
        seq_list.append((x, y))

def ring_decoder(bin, polygon_list, delta_size):
    # Extract number of chunks for a ring
    chks_in_ring = bytes_to_uint(bin, RING_CHK_CNT_SIZE)
    ring_coords = []
    # Loop through chunks in ring
    for i in range(chks_in_ring):
        sequence_decoder(bin, ring_coords, delta_size)
    polygon_list.append(ring_coords)

def polygon_decoder(bin, multipolygon_coords, delta_size):
    # Extract number of rings for a polygon
    rings_in_poly = bytes_to_uint(bin, POLY_RING_CNT_SIZE)
    polygon_coords = []
    # Loop through rings in polygon
    for _ in range(rings_in_poly):
        ring_decoder(bin, polygon_coords, delta_size)
    multipolygon_coords.append(shapely.Polygon(shell=polygon_coords[0], holes=polygon_coords[1:]))

def decode_header(bin):
    delta_size, type = struct.unpack_from('!BB', bin)
    type = GT(type)
    cfg.offset += 2 * 8 + 4 * FLOAT_SIZE  # Offset is 2 bytes for BB + 64 * 4 for bounding box
    algos.fpd_extended_lib.intersection_chunk_bbox_wrapper.intersection_skip_header(bin) # Circular import

    return delta_size, type

def fp_delta_decoding(bin_in):
    cfg.offset = 0
    bin = bitarray(endian='big')
    bin.frombytes(bin_in)

    delta_size, type = decode_header(bin)

    binary_length = len(bin)
    coords = []
    if type == GT.LINESTRING:
        while (cfg.offset + EOF_THRESHOLD <= binary_length):  # While != EOF
            sequence_decoder(bin, coords, delta_size)
        geometry = shapely.LineString(coords)

    elif type == GT.POLYGON:
        while (cfg.offset + EOF_THRESHOLD <= binary_length):  # While != EOF, i.e. at least one byte left
            ring_decoder(bin, coords, delta_size)
        geometry = shapely.Polygon(shell=coords[0], holes=coords[1:])

    elif type == GT.MULTIPOLYGON:
        while (cfg.offset + EOF_THRESHOLD <= binary_length):  # While != EOF
            polygon_decoder(bin, coords, delta_size)
        geometry = shapely.MultiPolygon(coords)
    return geometry

def decompress(self, bin):
    s = time.perf_counter()
    geometry = fp_delta_decoding(bin)
    t = time.perf_counter()
    return t - s, geometry