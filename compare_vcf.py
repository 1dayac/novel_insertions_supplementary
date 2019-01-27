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


def get_len(line):
    start_pos = line.find("SVLEN")
    ans = ""
    start_pos += 6
    while line[start_pos].isdigit():
        ans += line[start_pos]
        start_pos += 1
    if ans == "":
        return 0
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
with open("./../pacbio.txt", "r") as pacbio:
     for r in pacbio.readlines():
         break
         #if r[0] == "#":
         #    continue
         if r.split("\t")[1].split("/")[0] not in sv_dict:
             sv_dict[r.split("\t")[1].split("/")[0]] = []
         sv = SV(r.split("\t")[1].split("/")[0], int(r.split("\t")[1].split("/")[1]), len(r.split("\t")[2]), r.split("\t")[2])
         if len(r.split("\t")[2]) > 300:
            sv_dict[r.split("\t")[1].split("/")[0]].append(sv)
            len_300 += 1

print(len_300)
print(sv_dict)
# print(len(sv_dict))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000 = 0

with open("./../results/19240/pacbio.vcf", "r") as pacbio:
     for r in pacbio.readlines():
         #break
         if r.startswith("#") or r.find("DEL") != -1:
             continue
         splitted = r.split("\t")
         sv = SV(splitted[0], int(splitted[1]), get_len(r))
         if get_len(r) < 300:
             pass
         elif get_len(r) < 300:
             len_50_300 += 1
         elif get_len(r) < 500:
             len_300_500 += 1
         elif get_len(r) < 1000:
             len_500_1000 += 1
         elif get_len(r) < 2000:
             len_1000_2000 += 1
         else:
             len_2000 += 1


         if splitted[0] not in sv_dict:
             sv_dict[splitted[0]] = []
         sv_dict[splitted[0]].append(sv)


print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print(">=2000 " + str(len_2000))

pamir_dict = {}
near = 0
not_near = 0

# print(len(sv_dict))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000 = 0
try:
    with open("./../results/19240/pamir.vcf", "r") as pamir_vcf:
        for r in pamir_vcf.readlines():
            #break
            if r.startswith("#"):
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, get_len(r), get_seq(r))


            found = False
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if chrom not in pamir_dict:
                pamir_dict[chrom] = []
            pamir_dict[chrom].append(new_sv)
            if get_len(r) < 50:
                pass
            elif get_len(r) < 300:
                len_50_300 += 1
            elif get_len(r) < 500:
                len_300_500 += 1
            elif get_len(r) < 1000:
                len_500_1000 += 1
            elif get_len(r) < 2000:
                len_1000_2000 += 1
            else:
                len_2000 += 1

            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if Near(sv, new_sv):
                    if get_len(r) >= 300:
                        near += 1
                        sv.checked = True

                            #align_sequences(sv.seq, new_sv.seq)
                    #seq = r.split("\t")[7].split(";")[5][4:]
                    #print(seq)
                        print(chrom)
                        print(pos)
except:
    pass
for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False


print("pamir")
print(near)


print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print(">=2000 " + str(len_2000))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000 = 0

popins_dict = {}
near = 0
not_near = 0
ins = 0
try:

    with open("./../results/19240/popins.vcf", "r") as popins_vcf:
        for r in popins_vcf.readlines():
            if r.startswith("#"):
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, get_popins_len(r))
            found = False

            if get_popins_len(r) >= 50:
                ins += 1
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if chrom not in popins_dict:
                popins_dict[chrom] = []
            popins_dict[chrom].append(new_sv)

            if get_popins_len(r) < 50:
                pass
            elif get_popins_len(r) < 300:
                len_50_300 += 1
            elif get_popins_len(r) < 500:
                len_300_500 += 1
            elif get_popins_len(r) < 1000:
                len_500_1000 += 1
            elif get_popins_len(r) < 2000:
                len_1000_2000 += 1
            else:
                len_2000 += 1
            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if get_popins_len(r) >= 300:
                    if Near(sv, new_sv):
                        near += 1



                        #seq = r.split("\t")[7].split(";")[5][4:]
                        #print(seq)
                        sv.checked = True
                        print(chrom)
                        print(pos)
                        break
except:
    pass

for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False
print("popins")
print(near)
print(ins)

print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print(">=2000 " + str(len_2000))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000 = 0


def get_len_nui(line):
    return int(line[4])


near = 0

try:
    with open("./../results/19240/nui.txt", "r") as nui:
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
            elif get_len_nui(splitted) < 1000:
                len_500_1000 += 1
            elif get_len_nui(splitted) < 2000:
                len_1000_2000 += 1
            else:
                len_2000 += 1


            for sv in sv_dict[chrom]:
                if Near(sv, new_sv):
                    if sv.checked:
                        pass
                        continue

                    sv.checked = True
                    break
except:
    pass

print("NUI")

print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print(">=2000 " + str(len_2000))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000 = 0

for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False


total = 0
near = 0
not_near = 0
num = 0
ins = 0

lengths = []
anchors = []
with open("./../results/19240/novelx.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, len(r.split("\t")[4]), r.split("\t")[4])

        if chrom.startswith("chrUn_"):
            continue
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
        elif len(r.split("\t")[4]) < 1000:
            len_500_1000 += 1
        elif len(r.split("\t")[4]) < 2000:
            len_1000_2000 += 1
        else:
            len_2000 += 1

        for sv in sv_dict[chrom]:
            if sv.checked:
                continue
            if Near(sv, new_sv):


                if len(r.split("\t")[4]) >=  300:
                    found = True
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
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print(">=2000 " + str(len_2000))

print("Mean length - " + str(np.mean(lengths)))
print("Sum length - " + str(sum(lengths)))

print("Mean anchor sum - " + str(np.mean(anchors)))
print("Anchor sd - " + str(np.std(anchors)))
print("Max anchor sum - " + str(max(anchors)))
print("Max ins length - " + str(max(lengths)))



