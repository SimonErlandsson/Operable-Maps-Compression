from enum import Enum
import numpy as np
EM = Enum('EM', ['HUFFMAN', 'GOLOMB', 'AUTO', 'NONE'])
CM = Enum('CM', ['ZLIB', 'GZIP', 'ARITHMETIC', 'NONE'])
required_bits = lambda x: int(np.ceil(np.log2(x + 1)))

# If true: use "normal" FPD, else use "32-bit integer reprs.".
USE_DEFAULT_DOUBLE = False

# Enable per delta entropy compression
ENTROPY_METHOD = EM.NONE # 'HUFFMAN', 'GOLOMB', 'AUTO', 'NONE'.
# Whole chunk compression (all deltas in one, applied after entropy coding if active)
CHUNK_COMP_METHOD = CM.NONE # 'ZLIB', 'GZIP', 'NONE'

# Size of data-structures
D_BITSIZE_SIZE = 12 # Size used to store the size of the chunks OPTIMIZE
POLY_RING_CNT_SIZE = required_bits(2048)
RING_CHK_CNT_SIZE = 30
MAX_NUM_DELTAS = 31  # Max number of deltas in a chunk before split
D_CNT_SIZE = required_bits(MAX_NUM_DELTAS + 1)

# Disable optimized FPDE operations
DISABLE_OPTIMIZED_INTERSECTION = False
DISABLE_OPTIMIZED_BOUNDING_BOX = False

# Do not tweak manually
binary_length = 0 # Used when parsing
offset = 0  # Used when parsing
FLOAT_SIZE = 64 if USE_DEFAULT_DOUBLE else 32
ENTROPY_PARAM = None # Per geometry Golomb parameter
USE_FPINT = not USE_DEFAULT_DOUBLE
USE_ENTROPY = ENTROPY_METHOD != EM.NONE
COMPRESS_CHUNK = CHUNK_COMP_METHOD != CM.NONE


# Use "custom float size": TODO: Seems broken now
EXPONENT = 6 # Always defined, but not used if below are commented out
# FLOAT_SIZE = 35
# USE_FPINT = False
# USE_DEFAULT_DOUBLE = False

# Can be tweaked, though should not be needed to
EOF_THRESHOLD = D_CNT_SIZE + FLOAT_SIZE * 2  # Number of bits required to continue parsing