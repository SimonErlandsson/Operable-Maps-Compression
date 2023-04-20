USE_DEFAULT_DOUBLE = False
USE_DOUBLE_INTEGER = True


#Only seem to work for USE_DOUBLE_INTEGER = True with 32 bit FLOAT_SIZE (offset error)
COMPRESS_CHUNK = False
COMPRESSION_METHOD = "zlib" #"zlib" or "gzip"

USE_ENTROPY = False
CODES = None
DECODE_TREES = None

FLOAT_SIZE = 64
EXPONENT = 6
D_CNT_SIZE = 10 
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16
MAX_NUM_DELTAS = 15  # Max number of deltas in a chunk before split
EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing

binary_length = 0 # Used when parsing
offset = 0  # Used when parsing