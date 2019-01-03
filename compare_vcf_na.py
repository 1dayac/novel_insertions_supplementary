import sys
from Bio import pairwise2

import numpy as np
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

def get_len_nui(line):
    return int(line[4])

def get_len(line):
    start_pos = line.find("SVLEN")
    if start_pos == -1:
        return 0
    ans = ""
    start_pos += 6
    while line[start_pos].isdigit() or line[start_pos] == '-' or line[start_pos] == '+':
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

def Near(sv1, sv2):
    return abs(sv1.pos - sv2.pos) <= 100
len_300 = 0

print(len_300)
print(sv_dict)
# print(len(sv_dict))
len_50_300 = 0
len_300_500 = 0
len_500 = 0
with open("./../results/NA/remapped_pacbio.vcf", "r") as pacbio:
     for r in pacbio.readlines():
         if r.startswith("#") or r.find("DEL") != -1:
             continue
         splitted = r.split("\t")
         length = get_len(r)
         if length == 0:
             continue
         sv = SV(splitted[0], int(splitted[1]), get_len(r))
         if get_len(r) < 300:
             continue
         elif get_len(r) < 300:
             len_50_300 += 1
         elif get_len(r) < 500:
             len_300_500 += 1
         else:
             len_500 += 1


         if splitted[0] not in sv_dict:
             sv_dict[splitted[0]] = []
         sv_dict[splitted[0]].append(sv)

print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )

len_50_300 = 0
len_300_500 = 0
len_500 = 0
near = 0
with open("./../results/NA/nui.txt", "r") as nui:
    for r in nui.readlines():
        splitted = r.split("\t")
        length = get_len_nui(splitted)
        if length == 0:
            continue
        new_sv = SV("chr" + splitted[0], int(splitted[1]), get_len_nui(splitted))

        chrom = "chr" + splitted[0]
        pos = int(splitted[1])
        if get_len_nui(splitted) < 300:
            len_50_300 += 1
        elif get_len_nui(splitted) < 500:
            len_300_500 += 1
        else:
            len_500 += 1

        for sv in sv_dict[chrom]:
            if sv.checked:
                pass
                continue
            if Near(sv, new_sv):


                found = True
                near += 1
                #align_sequences(sv.seq, new_sv.seq)
                # seq = r.split("\t")[7].split(";")[5][4:]
                # print(seq)
                sv.checked = True
                print(chrom)
                print(pos)

print("NUI")
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )



total = 0
near = 0
not_near = 0
num = 0
ins = 0

lengths = []
anchors = []
with open("./../results/NA/novelx.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, len(r.split("\t")[4]), r.split("\t")[4])


        total += len(r.split("\t")[4])
        num += 1
        found = False
        if chrom not in sv_dict:
            sv_dict[chrom] = []
        if len(r.split("\t")[4]) >= 300:
            ins += 1
            lengths.append(len(r.split("\t")[4]))
            anchors.append(int(r.split("\t")[9]) + int(r.split("\t")[10]))
        if len(r.split("\t")[4]) < 300:
            len_50_300 += 1
        elif len(r.split("\t")[4]) < 500:
            len_300_500 += 1
        else:
            len_500 += 1
        for sv in sv_dict[chrom]:
            if sv.checked:
                continue
            if Near(sv, new_sv):
                pass


                near += 1
                #align_sequences(sv.seq, new_sv.seq)
                # seq = r.split("\t")[7].split(";")[5][4:]
                # print(seq)
                sv.checked = True
                print(chrom)
                print(pos)


print(ins)

print("Novel-X")
print("Total - " + str(ins))
print("Shared - " + str(near))
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )

print("Mean length - " + str(np.mean(lengths)))
print("Mean anchor sum - " + str(np.mean(anchors)))
print("Anchor sd - " + str(np.std(anchors)))
print("Max anchor sum - " + str(max(anchors)))
print("Max ins length - " + str(max(lengths)))



