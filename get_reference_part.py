import sys
from Bio import SeqIO, SeqRecord
#usage:
#python get_reference_part.py <reference.fasta> <chromosome name> <start_position> <end_position> <out_dir>
#or
#python get_reference_part.py <reference.fasta> <chrname_start_end> <out_dir>
if len(sys.argv) == 6:
    chr_name = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    out_dir = sys.argv[5]
elif len(sys.argv) == 4:
    splitted = sys.argv[2].split("_")
    chr_name = splitted[0]
    start = int(splitted[1])
    end = int(splitted[2])
    out_dir = sys.argv[3]

subseq = None
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if record.id == chr_name:
        subseq = record.seq[start:end]

out_filename = out_dir + "/" + chr_name + "_" + str(start) + "_" + str(end) + ".fasta"
SeqIO.write([SeqRecord("subseq", subseq)], out_filename, "fasta")


