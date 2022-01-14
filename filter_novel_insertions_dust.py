from Bio import pairwise2

class SV:
    def __init__(self, chrom, pos, length, seq = ""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.seq = seq.upper()


def align_sequences(s1, s2):
    alignments = pairwise2.align.globalxx(s1, s2)
    print(float(alignments[0][2])/max(len(alignments[0][0]), len(alignments[0][1])))
    print(pairwise2.format_alignment(*alignments[0]))

sv_dict = {}


def get_seq(line):
    start_pos = line.find("SEQ")
    ans = ""
    start_pos += 4
    while line[start_pos].isalpha():
        ans += line[start_pos]
        start_pos += 1
    return ans


def get_len(line):
    start_pos = line.find("SVLEN")
    if start_pos == -1:
        return 0
    ans = ""
    start_pos += 6
    while line[start_pos].isdigit() or line[start_pos] == "+" or line[start_pos] == "-":
        ans += line[start_pos]
        start_pos += 1
    if ans == "":
        return 0
    return abs(int(ans))

def get_popins_len(line):
    start_pos = line.find("length_")
    ans = ""
    start_pos += 7
    while line[start_pos].isdigit():
        ans += line[start_pos]
        start_pos += 1
    return int(ans) - 55

positions = []
with open("results/NA/records_to_keep.txt", "r") as dust:
    for r in dust.readlines():
        positions.append(r.strip())

print(len(positions))

with open("results/NA/pacbio.bed", "r") as pacbio:
    with open("results/NA/pacbio_filtered.bed", "w") as pacbio_filtered:

        for r in pacbio.readlines():
             if r.find("DEL") != -1:
                 continue
             if r.startswith("#"):
                 pacbio_filtered.write(r)
                 continue
             splitted = r.split("\t")
             if splitted[0] + "_" + splitted[1] not in positions:
                 pacbio_filtered.write(r)
