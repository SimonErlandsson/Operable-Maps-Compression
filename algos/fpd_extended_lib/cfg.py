# If true: use "normal" FPD, else use "32-bit integer reprs.".
USE_DEFAULT_DOUBLE = True

USE_ENTROPY = False
AUTO_SELECT_METHOD = False
ENTROPY_METHOD = "Huffman" #Golomb or Huffman
GOLOMB_MIN_BRUTEFORCE = False

# Size of data-structures
FLOAT_SIZE = 64 if USE_DEFAULT_DOUBLE else 32
EXPONENT = 6
D_CNT_SIZE = 16
D_BITSIZE_SIZE = 32 # Size used to store the size of the chunks OPTIMIZE
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16
MAX_NUM_DELTAS = 15  # Max number of deltas in a chunk before split
EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing

# Whole chunk compression (all deltas in one)
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
USE_DOUBLE_INTEGER = not USE_DEFAULT_DOUBLE