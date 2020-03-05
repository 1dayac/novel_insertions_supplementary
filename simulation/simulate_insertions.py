#Usage
#simulate_insertions.py num_of_insertions min_len max_len out_file

import numpy as np
import sys
from random import choice
from random import random, randrange

chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22", "chrX", "chrY"]

lengths = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
           135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
           46709983, 50818468, 156040895]

cum_lengths = np.cumsum(lengths)

def get_chrom_pos(pos):
    for i in range(len(chrs)):
        if cum_lengths[i] > pos:
            if i != 0:
                return chrs[i], pos - cum_lengths[i-1]
            else:
                return chrs[i], pos


num_of_insertions = int(sys.argv[1])
min_length = int(sys.argv[2])
max_length = int(sys.argv[3])
out_file = sys.argv[4]


with open(out_file, 'w') as out:
    for i in range(num_of_insertions):
        chrom, pos = get_chrom_pos(randrange(sum(lengths)))
        ins_len = randrange(min_length, max_length)
        ins = ''.join(choice(["A", "C", "G", "T"]) for i in range(ins_len))
        out.write(str(chrom) + "/" + str(pos) + "\t" + ins + "\n")
