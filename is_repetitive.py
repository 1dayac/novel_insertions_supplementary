from Bio import SeqIO
import sys

genome = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

def get_len(line):
    start_pos = line.find("SVLEN")
    ans = ""
    start_pos += 6
    while line[start_pos].isdigit():
        ans += line[start_pos]
        start_pos += 1
    return int(ans)
def get_popins_len(line):
    start_pos = line.find("length_")
    ans = ""
    start_pos += 7
    while line[start_pos].isdigit():
        ans += line[start_pos]
        start_pos += 1
    return int(ans) - 55
def get_my_len(s):
    return len(s.split("\t")[4])

is_repeat = 0
non_repeat = 0
with open("insertions_setcover.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        if get_len(r) < 300:
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        s = genome[chrom][pos - 300 : pos + 300].seq
        if all(s[i].islower() for i in range(len(s))):
            is_repeat += 1
        else:
            non_repeat += 1

print(is_repeat)
print(non_repeat)