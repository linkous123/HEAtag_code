import pysam
from argparse import ArgumentParser
import os
from collections import Counter
import pandas as pd
import time
from memory_profiler import profile
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor

parser = ArgumentParser(description='Read demultiplexer')
parser.add_argument('-r1', help='input sampletag r1 file')
parser.add_argument('-r2', help='input sampletag r2 file')
parser.add_argument('-c', help='input config file')
parser.add_argument('-b', help='input barcode whitelidt file')
parser.add_argument('-p', help='the prefix of output file')
parser.add_argument('-t', help='the number of threads used')
args = parser.parse_args()
inr1 = args.r1
inr2 = args.r2
cfig = args.c
outf = args.p + ".csv"
barcode=args.b
threads_used= args.t

def hamming_distance(str1, str2, threshold=1):
    """计算两个字符串的汉明距离"""
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))
 
def char_to_hamming_binary_char(char):
    """
    Convert a nucleotide to its corresponding binary array.

    Parameters:
    - nucleotide (str): A single nucleotide ('A', 'T', 'C', 'G', or 'N').

    Returns:
    - list: Binary array representing the input nucleotide.
    """
    dict_mapping = {'A': [0, 0, 0, 0, 1], 'T': [0, 0, 0, 1, 0], 'C': [0, 0, 1, 0, 0], 'G': [0, 1, 0, 0, 0], 'N': [1,0, 0, 0, 0]}

    return dict_mapping[char]



def char_to_hamming_binary_array(strings):
    """
    Convert a list of DNA strings to their corresponding Hamming binary arrays.

    Parameters:
    - strings (list): A list of DNA strings.

    Returns:
    - numpy.ndarray: Hamming binary arrays representing the input DNA strings.
    """
    result = np.array([[char_to_hamming_binary_char(char) for char in string] for string in strings], dtype=np.uint8)
    return result.reshape(result.shape[0], -1)


def hamming_group(observed_sample, sample_group):
    sample = observed_sample[np.newaxis, :]
    matrix = (np.abs(sample_group - sample) != 0).sum(axis=1, dtype=np.uint8)
    return matrix

def hamming_matrix_parallel(check_arr, white_arr, num_processes=50):
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        matrices = list(executor.map(hamming_group, check_arr, [white_arr] * len(check_arr)))

    return np.vstack(matrices)


def hamming_distance_dict(set_check, set_white, threshold=0, threads=50):
    binary_check=char_to_hamming_binary_array(set_check)
    binary_white=char_to_hamming_binary_array(set_white)

    start_time_4 = time.time()
    matrix=hamming_matrix_parallel(binary_check, binary_white, threads)
    end_time_4 = time.time()
    print(f"Execution time for matrix processing: {end_time_4 - start_time_4} seconds")

    result = {}
    dict_C = {i: string for i, string in enumerate(set_white)}

    for b, m in zip(set_check, matrix):
        indices = np.where(m <= threshold * 2)[0]
        if indices.size > 0:
            min_index = np.argmin(m[indices])
            min_index_global = indices[min_index]
            values = dict_C[min_index_global]
            result[b] = values
        else:
            result[b] = "nomatch"
    return result

def parse_config_file(cfig):
    config_values = {}
    sbs = list()
    with open(cfig, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # 跳过空行和注释
            key, value = map(str.strip, line.split(None, 1))
            config_values[key] = value
    cbl = config_values['-cb_location'].split(':')
    sbl = config_values['-sb_location'].split(':')
    umil = config_values['-umi_location'].split(':')
    sbs = set(config_values['-sb_seq'].split(','))
    cbm = int(config_values['-cb_mismatch'])
    sbm = int(config_values['-sb_mismatch'])
    return cbl, sbl, umil, sbs, cbm, sbm



def get_barcode(barcode):
    print("start get barcode")
    barcodes = set()
    with open(barcode, 'r') as inf:
        for bc in inf:
            bc = bc.strip()
            barcodes.add(bc)
    return barcodes

def parse_reads(inr1, inr2, cbl, sbl, umil):
    print("parse_reads")
    mat = set()
    cbl=slice(int(cbl[0]), int(cbl[1]))
    sbl=slice(int(sbl[0]), int(sbl[1]))
    umil=slice(int(umil[0]), int(umil[1]))
    with pysam.FastxFile(inr1) as fin_a, pysam.FastxFile(inr2) as fin_b:
        for reada in fin_a:
            readb = next(fin_b)

            if reada.name != readb.name:
                print("read names not equal.")
                os.exit()
            else:
                cbc = reada.sequence[cbl]
                sbc = readb.sequence[sbl]
                umi = reada.sequence[umil]
                mat.add((umi, sbc, cbc))
    return mat

def modify(mat, included_sbc, included_cbc, outf, sbm, cbm, threads):
    print("modify")
    print(f"The number of multiple threads used{threads}")
    mat_dict = Counter([(sbc, cbc) for _, sbc, cbc in mat])
    values=list(mat_dict.values())
    keys_list = list(mat_dict.keys())
    sbc, cbc = zip(*keys_list)

    sbc_set=set(sbc)
    cbc_set=set(cbc)

    #with open("sbc_set.txt", 'w') as file:
    #    for element in sbc_set:
    #        file.write(element + '\n')


    if cbm > 0:
        #print(len(cbc_set - included_cbc))
        cb_dict = hamming_distance_dict(cbc_set - included_cbc, included_cbc, threshold=cbm, threads=threads)
        cbc=[cb_dict.get(item, item) for item in cbc]
    else:
        cbc=[item if item in included_cbc else "nomatch" for item in cbc]
    
    if sbm >0:
        sb_dict = hamming_distance_dict(sbc_set - included_sbc, included_sbc, threshold=sbm)
        sbc=[sb_dict.get(item, item) for item in sbc]
    else:
        sbc= [item if item in included_sbc else "nomatch" for item in sbc]
    
    data = {'cellbarcode': cbc, 'Column': sbc, 'Value': values}
    df = pd.DataFrame(data)
    df = df.pivot_table(index='cellbarcode', columns='Column', values='Value', aggfunc='sum',fill_value=0)
    df.sort_index(inplace=True)
    df.to_csv(outf)


if __name__=="__main__":

    cbl, sbl, umil, sbs, cbm, sbm = parse_config_file(cfig)

    start_time = time.time()   
    mat=parse_reads(inr1, inr2, cbl, sbl, umil)
    end_time = time.time()
    print(f"Execution time for barcodes df processing: {end_time - start_time} seconds")


    start_time_1 = time.time()
    barcodes=get_barcode(args.b)
    end_time_1 = time.time()
    print(f"Execution time for barcodes df processing: {end_time_1 - start_time_1} seconds")

    start_time_2 = time.time()
    modify(mat, sbs, barcodes, outf, sbm, cbm, threads_used)
    end_time_2 = time.time()
    print(f"Execution time for modify df processing: {end_time_2 - start_time_2} seconds")