required_bits = lambda x: int(np.ceil(np.log2(x + 1)))
import numpy as np

# If true: use "normal" FPD, else use "32-bit integer reprs.".
USE_DEFAULT_DOUBLE = False

USE_ENTROPY = True
AUTO_SELECT_METHOD = False # TODO: Seems broken
ENTROPY_METHOD = "Huffman" #Golomb or Huffman
GOLOMB_MIN_BRUTEFORCE = False

# Size of data-structures
FLOAT_SIZE = 64 if USE_DEFAULT_DOUBLE else 32
EXPONENT = 6
D_BITSIZE_SIZE = 20 # Size used to store the size of the chunks OPTIMIZE
POLY_RING_CNT_SIZE = required_bits(127)
RING_CHK_CNT_SIZE = 10
MAX_NUM_DELTAS = 31  # Max number of deltas in a chunk before split
D_CNT_SIZE = required_bits(MAX_NUM_DELTAS + 1)

EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing

# Whole chunk compression (all deltas in one, applied after entropy coding)
COMPRESS_CHUNK = False
COMPRESSION_METHOD = "zlib" #"zlib" or "gzip"

# Disable optimized FPDE operations
DISABLE_OPTIMIZED_INTERSECTION = False
DISABLE_OPTIMIZED_BOUNDING_BOX = False

# Do not tweak manually
binary_length = 0 # Used when parsing
offset = 0  # Used when parsing
CODES, DECODE_TREES = None, None # Used as global data
ENTROPY_STATE = (None, None, None)
ENTROPY_PARAM = None
ARITHMETIC_ENCODER = None
USE_FPINT = not USE_DEFAULT_DOUBLE