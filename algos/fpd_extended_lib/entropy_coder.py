from collections import defaultdict
import bench_utils
from bitarray import bitarray, util, decodetree
from shapely.geometry import shape
from collections import defaultdict
import math
from algos.fpd_extended_lib.low_level import *
from algos.fpd_extended_lib.helpers import calculate_delta_size
from algos.fpd_extended_lib.cfg import *
import algos.fpd_extended_lib.cfg as cfg

import pickle

def set_entropy_code(path):
    df, unary_idxs = bench_utils.read_dataset(path)
    deltas_by_bits = defaultdict(list)
    for idx in unary_idxs: # List of single idxs
        opt_size, _, delta_list = calculate_delta_size(shape(df.iloc[idx]), True)
        for delta in delta_list[1]:
            #If coordinate not delta_encoded
            if delta != 0 and math.log2(delta) > opt_size:
                continue
            deltas_by_bits[opt_size].append(uint_to_ba(delta, opt_size).to01())
          
    delta_freqs = __get_bit_seq_freqs(deltas_by_bits)
    codes = {}
    decode_trees = {}
    for opt_size in delta_freqs:
        codes[opt_size] = __get_entropy_codes(delta_freqs[opt_size])
        decode_trees[opt_size] = decodetree(codes[opt_size])
    
    cfg.CODES = codes
    cfg.DECODE_TREES = decode_trees


def __get_entropy_codes(frequency_list):
        return util.huffman_code(frequency_list)
    
def encode(msg, delta_len):
    value = cfg.CODES[delta_len][msg.to01()]
    return value, len(value)

#decode also defined in low level (circular import)
def decode(msg, delta_len):
    return decode_msg(msg, delta_len)              
    
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

