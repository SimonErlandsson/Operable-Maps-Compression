from bitarray import bitarray, bits2bytes, util
from var_float import VarFloat
from algos.fpd_extended_lib.cfg import *
import algos.fpd_extended_lib.cfg as cfg
import struct

var_float = VarFloat(EXPONENT, FLOAT_SIZE)

def zz_encode(num):
    return -2 * num - 1 if num < 0 else 2 * num

def zz_decode(num):
    return - (num + 1) // 2 if num & 1 == 1 else num // 2

# Split float to int
def double_as_fpint(num):
    integral_part, decimal_part = format(num + 180, '.7f').split(".")
    res = int(integral_part + decimal_part)
    return res

def fpint_as_double(num):
    str_nbr = str(num)
    integral_part, decimal_part = str_nbr[:-7], str_nbr[-7:]
    res = round(float(integral_part + "." + decimal_part) - 180, 7)
    return res

# Decode
def bin_to_double(bin):
    if USE_DEFAULT_DOUBLE:
        return struct.unpack('!d', bin)[0]
    elif USE_FPINT:
        return long_as_double(util.ba2int(bin, signed=False))
    else:
        return var_float.bin_to_float(bin)
    
def long_as_double(long):
    if USE_DEFAULT_DOUBLE:
        return struct.unpack('!d', struct.pack('!q', long))[0]
    elif USE_FPINT:
        return fpint_as_double(long)
    else:
        return var_float.bin_to_float(var_float.long_to_bits(long))

def bytes_to_decoded_coord(bin, prev_coord, input_size=64):
    from algos.fpd_extended_lib.entropy_coder import decode
    if cfg.USE_ENTROPY:
        bin, input_size = decode(bin[cfg.offset:], input_size)
    else: 
        bin = bin[cfg.offset: cfg.offset + input_size]

    val = util.ba2int(bin, signed=False)
    val = zz_decode(val) + double_as_long(prev_coord)
    val = long_as_double(val)
    cfg.offset += input_size
    return val

def bytes_to_double(bin):
    bin = bin[cfg.offset:cfg.offset + FLOAT_SIZE]
    val = bin_to_double(bin)
    cfg.offset += FLOAT_SIZE
    return val

def bytes_to_uint(bin, len):
    val = util.ba2int(bin[cfg.offset:cfg.offset + len], signed=False)
    cfg.offset += len
    return val

# Encode
def double_to_bytes(x):
    if USE_DEFAULT_DOUBLE:
        return struct.pack("!d", x)
    elif USE_FPINT:
        return uint_to_ba(double_as_long(x), 32)
    else:
        return var_float.float_to_bin(x)
    
def double_as_long(num):
    if USE_DEFAULT_DOUBLE:
        return struct.unpack('!q', struct.pack('!d', num))[0]
    elif USE_FPINT:
        return double_as_fpint(num)
    else:
        return var_float.bits_to_long(var_float.float_to_bin(num))

# Inline refactorized from https://github.com/ilanschnell/bitarray/blob/master/bitarray/util.py
def uint_to_ba(x, length):
    if x == 0:
        return util.zeros(length or 1, "big")

    a = bitarray(0, "big")

    a.frombytes(x.to_bytes(bits2bytes(x.bit_length()), byteorder="big"))
    la = len(a)
    if la == length:
        return a

    return a[-length:] if la > length else util.zeros(length - la, "big") + a

def uchar_to_bytes(x):
    return x.to_bytes(1, 'big')