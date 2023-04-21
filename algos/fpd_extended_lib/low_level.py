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
def double_as_float_int(num):
    neg = False
    num = str(num)
    #num = format(num, '.80f')
    if "." not in num:
        num = num + ".0"
    [integral_part, decimal_part] = str(num).split(".")
    if integral_part[0] == "-":
        neg = True
        integral_part = integral_part[1:]

    decimal_part = decimal_part[0:min(7, len(decimal_part))]
    decimal_part = decimal_part+ str(0) * (7 - len(decimal_part))
    res = int(integral_part + decimal_part)
    return  zz_encode(res * (-1) if neg else res)

def float_int_as_double(num):
    neg = False
    num = str(zz_decode(num))
    if num[0] == "-":
        neg = True
        num = num[1:]
    integral_part, decimal_part= num[:-7], num[-7:]
    integral_part = str(integral_part)
    decimal_part = str(decimal_part).rstrip("0")
    num = float(integral_part + "." + decimal_part)
    #print(integral_part, decimal_part)
    return num * (-1) if neg else num

# Decode
def bin_to_double(bin):
    if USE_DEFAULT_DOUBLE:
        return struct.unpack('!d', bin)[0]
    elif USE_DOUBLE_INTEGER:
        return long_as_double(util.ba2int(bin, signed=False))
    else:
        return var_float.bin_to_float(bin)

def bytes_to_decoded_coord(bin, prev_coord, input_size=64):
    from algos.fpd_extended_lib.entropy_coder import decode
    if USE_ENTROPY:
        bin , input_size  =  decode(bin[cfg.offset:], input_size)
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
def double_as_long(num):
    if USE_DEFAULT_DOUBLE:
        return struct.unpack('!q', struct.pack('!d', num))[0]
    elif USE_DOUBLE_INTEGER:
        return double_as_float_int(num)
    else:
        return var_float.bits_to_long(var_float.float_to_bin(num))

def long_as_double(num):
    if USE_DEFAULT_DOUBLE:
        return struct.unpack('!d', struct.pack('!q', num))[0]
    elif USE_DOUBLE_INTEGER:
        return float_int_as_double(num)
    else:
        return var_float.bin_to_float(var_float.long_to_bits(num))
    
def double_to_bytes(x):
    if USE_DEFAULT_DOUBLE:
        return struct.pack("!d", x)
    elif USE_DOUBLE_INTEGER:
        return uint_to_ba(double_as_long(x),32)
    else:
        return var_float.float_to_bin(x)

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