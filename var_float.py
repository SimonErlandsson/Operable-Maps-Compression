import bitarray
from bitarray import util, bits2bytes, bitarray
import math

class VarFloat:
    def __init__(self, exponent, float_size):
        self.EXPONENT = exponent
        self.FLOAT_SIZE = float_size

    def uint_to_ba(self, x, length):
            if x == 0:
                return util.zeros(length or 1, "big")

            a = bitarray(0, "big")

            a.frombytes(x.to_bytes(bits2bytes(x.bit_length()), byteorder="big"))
            la = len(a)
            if la == length:
                return a

            return a[-length:] if la > length else util.zeros(length - la, "big") + a

    def float_to_bin(self, float):
        decimal_part, integral_part = math.modf(float)
        sign = 0 if float > 0 else 1
        decimal_part = decimal_part * (1 - sign * 2)
        integral_part = int(integral_part) * (1 - sign * 2)

        integral_bin = bitarray()
        while(integral_part >= 1):
            integral_bin.append(integral_part % 2)
            integral_part = integral_part >> 1
        integral_bin.reverse()

        decimal_bin = bitarray()
        beg_bits_zero, had_1s = 0, False
        while(decimal_part != 0):
            decimal_part *= 2
            rem = 0 if decimal_part < 1 else 1
            if not had_1s and len(integral_bin) == 0:
                if rem == 0:
                    beg_bits_zero += 1
                    continue
                else:
                    had_1s = True
                    
            decimal_bin.append(rem)
            if decimal_part >= 1:
                decimal_part -= 1
        
        if len(integral_bin) == 0:
            exponent = (pow(2,(self.EXPONENT - 1)) - 1) - beg_bits_zero - 1
        else:
            exponent = (pow(2,(self.EXPONENT - 1)) - 1) + len(integral_bin) - 1

        res = bitarray()
        res.append(sign)
        res.extend(self.uint_to_ba(exponent, self.EXPONENT ))
        integral_bin.extend(decimal_bin)
        while len(integral_bin) < self.FLOAT_SIZE - (self.EXPONENT + 1) + 1:
            integral_bin.append(0)

        res.extend(integral_bin[1:self.FLOAT_SIZE - (self.EXPONENT + 1)  + 1])
        return res


    def bin_to_float(self, precision_float):
        sign = (1 + (-2 * precision_float[0]))
        exponent = util.ba2int(precision_float[1:(self.EXPONENT + 1) ], signed=False) - (pow(2,(self.EXPONENT - 1)) - 1)
        decimal_part = precision_float[(self.EXPONENT + 1) :]
        decimal = 1
        for i in range(1, self.FLOAT_SIZE - (self.EXPONENT + 1)  + 1):
            decimal += decimal_part[i - 1] * math.pow(2, -i)
        return round(decimal * sign * math.pow(2,exponent),7)


    def bits_to_long(self, bits):
        res = 0
        for ele in bits:
            res = (res << 1) | ele
        return res

    def long_to_bits(self, num):      
        binary = bin(num)[2:]
        res =  bitarray()
        res.extend('0' * (self.FLOAT_SIZE - len(binary)))
        res.extend(binary)
        return res