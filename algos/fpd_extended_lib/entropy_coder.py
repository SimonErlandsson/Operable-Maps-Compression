from collections import defaultdict
import bench_utils
from bitarray import bitarray, util, decodetree
from shapely.geometry import shape
from collections import defaultdict
import math
from algos.fpd_extended_lib.low_level import *
from algos.fpd_extended_lib.cfg import *
import algos.fpd_extended_lib.cfg as cfg


def set_entropy_code(path):
    from algos.fpd_extended_lib.compress import calculate_delta_size
    from algos.fpd_extended_lib.helpers import k_est
    if ENTROPY_METHOD == "Huffman":
        total_deltas = []
        df, unary_idxs = bench_utils.read_dataset(path)
        deltas_by_bits = defaultdict(list)
        for idx in unary_idxs: # List of single idxs
            opt_size, _, delta_list = calculate_delta_size(shape(df.iloc[idx]), True)
            for delta in delta_list[1]:
                #If coordinate not delta_encoded
                if delta != 0 and math.log2(delta) > opt_size:
                    continue
                total_deltas.append(delta)
                deltas_by_bits[opt_size].append(uint_to_ba(delta, opt_size).to01())
        
    if ENTROPY_METHOD == "Huffman":
        delta_freqs = __get_bit_seq_freqs(deltas_by_bits)
        codes = {}
        decode_trees = {}
        for opt_size in delta_freqs:
            codes[opt_size] = __get_entropy_codes(delta_freqs[opt_size])
            decode_trees[opt_size] = decodetree(codes[opt_size])
        cfg.CODES = codes
        cfg.DECODE_TREES = decode_trees

    # if ENTROPY_METHOD == "Golomb":
    #     cfg.ENTROPY_PARAM = k_est(total_deltas)
    #     print(cfg.ENTROPY_PARAM)


def __get_entropy_codes(frequency_list):
        return util.huffman_code(frequency_list)
    
def encode(msg, delta_len):
    if ENTROPY_METHOD == "Huffman":
        return huffman_encode(msg, delta_len)
    elif ENTROPY_METHOD == "Golomb":
        return golomb_encode(msg)


#decode also defined in low level (circular import)
def decode(msg, delta_len):
    if ENTROPY_METHOD == "Huffman":
        return huffman_decode(msg, delta_len)
    elif ENTROPY_METHOD == "Golomb":
        return golomb_decode(msg, delta_len)            
    

def golomb_encode(msg):
    value = util.ba2int(msg)
    M = int(math.pow(2,cfg.ENTROPY_PARAM))
    r = uint_to_ba(value & int(M - 1),cfg.ENTROPY_PARAM)
    q = __unary_encode(value >> cfg.ENTROPY_PARAM)
    return q + r, len(r + q)

def golomb_decode(msg, delta_len):
    M = math.pow(2,cfg.ENTROPY_PARAM) 
    q = 0
    while(msg[q] == 1):
        q += 1
    r = msg[q + 1:q + 1 + cfg.ENTROPY_PARAM]
    r = util.ba2int(r) # + 1 for the ending 0
    res =  q * M + r
    return uint_to_ba(int(res), delta_len), q + 1 + cfg.ENTROPY_PARAM


def huffman_decode(msg, delta_len):
    for i in range(1, len(msg) + 100):
        try:
            value = msg[:i].decode(cfg.DECODE_TREES[delta_len])[0]
            return bitarray(value), i
        except:
            pass    

def huffman_encode(msg, delta_len):
    value = cfg.CODES[delta_len][msg.to01()]
    return value, len(value)

def __unary_encode(value):
    return bitarray("1" * value + "0")

def __get_bit_seq_freqs(deltas):
    freq_by_bits = defaultdict(dict)
    for opt_size in deltas: 
        for delta in deltas[opt_size]:  
            delta_len = len(delta)
            if delta in freq_by_bits[delta_len]:
                freq_by_bits[delta_len][delta] += 1
            else:
                freq_by_bits[delta_len][delta] = 1
    return freq_by_bits
