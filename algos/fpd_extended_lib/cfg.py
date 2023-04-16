USE_DEFAULT_DOUBLE = False
USE_DOUBLE_INTEGER = True

USE_ENTROPY = True
CODES = None
DECODE_TREES = None

FLOAT_SIZE = 32
EXPONENT = 6
D_CNT_SIZE = 16
POLY_RING_CNT_SIZE = 16
RING_CHK_CNT_SIZE = 16
MAX_NUM_DELTAS = 15  # Max number of deltas in a chunk before split
EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing

offset = 0  # Used when parsing