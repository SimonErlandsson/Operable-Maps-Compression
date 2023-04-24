USE_DEFAULT_DOUBLE = False
USE_DOUBLE_INTEGER = True

USE_ENTROPY = False
ENTROPY_METHOD = "Golomb" #Golomb or Huffman
ENTROPY_PARAM = None
CODES, DECODE_TREES = None, None
GOLOMB_MIN_BRUTEFORCE = False

FLOAT_SIZE = 32
EXPONENT = 6
D_CNT_SIZE = 16
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16
MAX_NUM_DELTAS = 10000  # Max number of deltas in a chunk before split
EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing

#Only seem to work for USE_DOUBLE_INTEGER = True with 32 bit FLOAT_SIZE (offset error)
COMPRESS_CHUNK = False
COMPRESSION_METHOD = "zlib" #"zlib" or "gzip"


binary_length = 0 # Used when parsing
offset = 0  # Used when parsing
