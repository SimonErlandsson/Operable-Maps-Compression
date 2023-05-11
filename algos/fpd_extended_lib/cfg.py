

USE_DEFAULT_DOUBLE = True
USE_DOUBLE_INTEGER = False

USE_ENTROPY = False
AUTO_SELECT_METHOD, ENTROPY_STATE = False, (None, None, None)
ENTROPY_METHOD = "Golomb" #Golomb or Huffman
ENTROPY_PARAM = None
CODES, DECODE_TREES = None, None
GOLOMB_MIN_BRUTEFORCE = False

FLOAT_SIZE = 64 if USE_DEFAULT_DOUBLE else 32
EXPONENT = 6
D_CNT_SIZE = 32 if USE_DEFAULT_DOUBLE else 16
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16
MAX_NUM_DELTAS = 10  # Max number of deltas in a chunk before split
EOF_THRESHOLD = 2 * D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing


COMPRESS_CHUNK = False
COMPRESSION_METHOD = "arith" #"zlib" or "gzip"
CONTEXT_SIZE = 4
ARITHMETIC_ENCODER = None


OPTIMIZED = True
POINT_ERR_TOL = 1e-12 if USE_DEFAULT_DOUBLE else 1e-10

binary_length = 0 # Used when parsing
offset = 0  # Used when parsing
