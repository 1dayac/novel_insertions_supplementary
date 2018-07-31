import sys
from Bio import SeqIO, SeqRecord
#usage:
#python get_reference_part.py <reference.fasta> <chromosome name> <start_position> <end_position> <out_dir>
chr_name = sys.argv[2]
subseq = None
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if record.id == chr_name:
        subseq = record.seq[int(sys.argv[3]):int(sys.argv[4])]

out_filename = sys.argv[5] + "/" + chr_name + "_" + sys.argv[3] + "_" + sys.argv[4] + ".fasta"
SeqIO.write([SeqRecord("subseq", subseq)], out_filename, "fasta")


