import bench_utils
import math
import zlib
import gzip
import glob
import dill as pickle
from collections import defaultdict
from bitarray import bitarray, util, decodetree
from shapely.geometry import shape
from collections import defaultdict
from algos.fpd_extended_lib.low_level import *
import algos.fpd_extended_lib.cfg as cfg

ENTROPY_BASE_PATH = "data/entropy_models"
huffman_codecs = {}
huffman_decodecs = {}
## Load Huffman trees from disk
for f in glob.glob(ENTROPY_BASE_PATH + '/*'):
    with open(f, 'rb') as f:
            bit_size, codec = pickle.load(f)
            huffman_codecs[bit_size] = codec
            huffman_decodecs[bit_size] = decodetree(codec)

def uint_to_ba(x, length):
    from algos.fpd_extended_lib.low_level import uint_to_ba as uint_to_ba_optimized
    return uint_to_ba_optimized(x, length)
    
def encode(msg, delta_len):
    if cfg.ENTROPY_METHOD == cfg.EM.HUFFMAN:
        if delta_len in huffman_codecs:
            return huffman_encode(msg, delta_len)
        else:
            return msg, len(msg)
    elif cfg.ENTROPY_METHOD == cfg.EM.GOLOMB:
        return golomb_encode(msg)
        
def decode(msg, delta_len):
    if cfg.ENTROPY_METHOD == cfg.EM.HUFFMAN:
        if delta_len in huffman_decodecs:
            return huffman_decode(msg, delta_len)
        else:
            return msg[:delta_len], delta_len
    elif cfg.ENTROPY_METHOD == cfg.EM.GOLOMB:
        return golomb_decode(msg, delta_len)            

#Huffman Encoding
def huffman_decode(msg, delta_len):
    for i in range(1, len(msg) + 100):
        try:
            dt = msg[:i].decode(huffman_decodecs[delta_len])[0]
            if dt == -1:
                dt = msg[i:i + delta_len]
                return dt, i + delta_len
            else:
                return util.int2ba(dt, delta_len, endian='big', signed=False), i
        except:
            pass    

def huffman_encode(msg, delta_len): 
    dt_long = util.ba2int(msg, signed=False)
    if dt_long in huffman_codecs[delta_len]:
        bits = huffman_codecs[delta_len][dt_long]
    else:
        # Prepend 'missing value' symbol and add full delta
        bits = huffman_codecs[delta_len][-1] + msg.to01()
    return bits, len(bits)

#Golomb-Rise Encoding
def golomb_encode(msg):
    value = util.ba2int(msg)
    M = int(math.pow(2, cfg.ENTROPY_PARAM))
    r = uint_to_ba(value & int(M - 1), cfg.ENTROPY_PARAM)
    q = __unary_encode(value >> cfg.ENTROPY_PARAM)
    res = q + r
    return res, len(res)

def golomb_decode(msg, delta_len):
    M = math.pow(2,cfg.ENTROPY_PARAM) 
    q = 0
    while msg[q] == 1:
        q += 1
    r = msg[q + 1:q + 1 + cfg.ENTROPY_PARAM]
    r = util.ba2int(r) # + 1 for the ending 0
    res =  q * M + r
    return uint_to_ba(int(res), delta_len), q + 1 + cfg.ENTROPY_PARAM

def __unary_encode(value):
    return bitarray("1" * value + "0")

def k_est(deltas):
    delta_mean = sum(deltas) / len(deltas)
    golden_ratio = (math.sqrt(5) + 1) / 2
    try:
        return max(0, 1 + math.floor(math.log2(math.log(golden_ratio - 1) / math.log(delta_mean /(delta_mean + 1)))))
    except:
        return 0
    
#Overall
def get_entropy_metadata(deltas, delta_size):
    '''
        Returns data appended to header and settings used for compression/decomp.
    '''
    deltas = [d for d in deltas if (d == 0 or math.log2(d) <= delta_size)]
    if cfg.ENTROPY_METHOD == cfg.EM.NONE:
         return (cfg.EM.NONE, False, 0)
    if cfg.ENTROPY_METHOD == cfg.EM.AUTO:
        return get_best_strategy(deltas, delta_size)
    if cfg.ENTROPY_METHOD == cfg.EM.GOLOMB:
        return (cfg.EM.GOLOMB, True, k_est(deltas))
    if cfg.ENTROPY_METHOD == cfg.EM.HUFFMAN:
        return (cfg.EM.HUFFMAN, True, 255)
    
def get_best_strategy(deltas, delta_size):
    '''
        Calculates the compression method which gives the lowest total size.
    '''
    old_entropy_param = cfg.ENTROPY_PARAM
    golomb_estimated = k_est(deltas)
    cfg.ENTROPY_PARAM = golomb_estimated 
    golomb_tot_size, huffman_tot_size, normal_tot_size = 0, 0, 0
    
    for d in deltas:
        d_bits = uint_to_ba(d, delta_size)
        golomb_tot_size += len(golomb_encode(d_bits)[0])
        normal_tot_size += len(d_bits)
        if delta_size in huffman_codecs:
            huffman_tot_size += len(huffman_encode(d_bits, delta_size)[0])
    cfg.ENTROPY_PARAM = old_entropy_param
    
    if delta_size not in huffman_codecs:
        huffman_tot_size = 99999999999
    if golomb_tot_size < normal_tot_size and golomb_tot_size < huffman_tot_size:
        return (cfg.EM.GOLOMB, True, golomb_estimated)
    elif huffman_tot_size < normal_tot_size:
        return (cfg.EM.HUFFMAN, True, 255)
    else:
        return (cfg.EM.NONE, False, 0)

def decode_entropy_param(value, delta_size):
    if value == 0:
        cfg.USE_ENTROPY = False
        cfg.ENTROPY_METHOD = cfg.EM.NONE
        cfg.ENTROPY_PARAM = 0
    elif value == 255:
        cfg.USE_ENTROPY = True
        cfg.ENTROPY_METHOD = cfg.EM.HUFFMAN
        cfg.ENTROPY_PARAM = delta_size
    else:
        cfg.USE_ENTROPY = True
        cfg.ENTROPY_METHOD = cfg.EM.GOLOMB
        cfg.ENTROPY_PARAM = value

def train_arith_model(model, dataset, iter = None):
    '''
        NOTE: Arithmetic encoding was discountinued during development
        due to time constraints. But may be used instead of Huffman for
        encoding a whole chunk of deltas.
    '''
    import bench_utils, tqdm
    from algos.fpd_extended_lib.compress import calculate_delta_size
    
    df, unary_idxs = bench_utils.read_dataset(dataset)
    unary_idxs = list(set(unary_idxs))
    
    # Compress files, benchmark unaries
    cnt = 0
    iter == len(unary_idxs) if iter == None else iter
    for idx in tqdm.tqdm(range(iter)): # List of single idxs
        idx = unary_idxs[idx]
        opt_size, _, deltas = calculate_delta_size(shape(df.iloc[idx]), True)
        for d in deltas[1]:
            for bit in uint_to_ba(d, opt_size):
                model.update(bit)
        cnt += 1
        if iter != None and cnt > iter:
            break

# Chunk based compression
def compress_chunk(bits, chk_hdr_offset, delta_bytes_size):
    chk_dt_offset = chk_hdr_offset + cfg.D_CNT_SIZE + cfg.D_BITSIZE_SIZE + 2 * cfg.FLOAT_SIZE
    before_bits = bits[:chk_dt_offset]
    coords_bits = bits[chk_dt_offset:chk_dt_offset + delta_bytes_size]
    comp_dt_bits = bitarray()

    if cfg.CHUNK_COMP_METHOD == cfg.CM.ZLIB:
        comp_dt_bytes = zlib.compress(coords_bits.tobytes())
        comp_dt_bits.frombytes(comp_dt_bytes)

    elif cfg.CHUNK_COMP_METHOD == cfg.CM.GZIP:
        comp_dt_bytes = gzip.compress(coords_bits.tobytes())
        comp_dt_bits.frombytes(comp_dt_bytes)

    elif cfg.CHUNK_COMP_METHOD == cfg.CM.ARITHMETIC:
        comp_dt_bits = bitarray()
        comp_dt_bits.extend(cfg.ARITHMETIC_ENCODER.compress([int(x) for x in [*coords_bits.to01()]]))
        
    else:
        comp_dt_bytes = coords_bits.tobytes()
        comp_dt_bits.frombytes(comp_dt_bytes)

    before_bits += comp_dt_bits
    return before_bits, len(comp_dt_bits)


def decompress_chunk(bits, chk_dt_offset, chk_dt_bitsize):
    before_bits = bits[:chk_dt_offset]
    comp_dt_bits = bits[chk_dt_offset:chk_dt_offset + chk_dt_bitsize]
    after_bits = bits[chk_dt_offset + chk_dt_bitsize:]
    decompressed_bits = bitarray()

    if cfg.CHUNK_COMP_METHOD == cfg.CM.ZLIB:
        coords_bits = zlib.decompress(comp_dt_bits.tobytes())
        decompressed_bits.frombytes(coords_bits)

    elif cfg.CHUNK_COMP_METHOD == cfg.CM.GZIP:
        coords_bits = gzip.decompress(comp_dt_bits.tobytes())
        decompressed_bits.frombytes(coords_bits)

    elif cfg.CHUNK_COMP_METHOD == cfg.CM.ARITHMETIC:
        indata = [int(x) for x in [*comp_dt_bits.to01()]]
        decompressed_bits.extend(cfg.ARITHMETIC_ENCODER.decompress(indata, len(indata)))
    else:
        decompressed_bits = comp_dt_bits
    
    chk_header_offset = chk_dt_offset - 2 * cfg.FLOAT_SIZE - cfg.D_BITSIZE_SIZE - cfg.D_CNT_SIZE
    before_bits[chk_header_offset + cfg.D_CNT_SIZE:chk_header_offset + cfg.D_CNT_SIZE + cfg.D_BITSIZE_SIZE] = int_to_ba(0, cfg.D_BITSIZE_SIZE) 

    before_bits += decompressed_bits
    before_bits += after_bits
    cfg.binary_length = len(before_bits)
    return before_bits, len(decompressed_bits)