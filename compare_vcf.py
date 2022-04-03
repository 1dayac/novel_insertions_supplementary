import sys
from Bio import pairwise2
try:
    dataset = sys.argv[1]
except:
    dataset = "HG002_stlfr"

try:
    dataset_simple = dataset[:dataset.index("_")]
except:
    dataset_simple = dataset
postfix = ""

import numpy as np
class SV:
    def __init__(self, chrom, pos, length, seq = ""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.seq = seq.upper()

    def __str__(self):
        return ">" + self.chrom + "_" + str(self.pos) + "\n" + str(self.seq.strip())


def align_sequences(s1, s2):
    alignments = pairwise2.align.globalxx(s1, s2)
    print(float(alignments[0][2])/max(len(alignments[0][0]), len(alignments[0][1])))
    return float(alignments[0][2])/max(len(alignments[0][0]), len(alignments[0][1]))
    #print(pairwise2.format_alignment(*alignments[0]))

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

def get_len_bed(line):
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
    return int(ans) - 47

def Near(sv1, sv2):
    return abs(sv1.pos - sv2.pos) <= 100


len_300 = 0
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0

total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0



try:
    with open("./results/simulated/insertions.txt", "r") as pacbio:
         for r in pacbio.readlines():
             if not dataset.startswith("simulated"):
                 break
             #if r[0] == "#":
             #    continue
             if r.split("\t")[0].split("/")[0] not in sv_dict:
                 sv_dict[r.split("\t")[0].split("/")[0]] = []
             sv = SV(r.split("\t")[0].split("/")[0], int(r.split("\t")[0].split("/")[1]), len(r.split("\t")[1]), r.split("\t")[1])
             if len(r.split("\t")[1]) >= 50:
                sv_dict[r.split("\t")[0].split("/")[0]].append(sv)
                len_300 += 1
                if len(r.split("\t")[1]) < 300:
                    len_50_300 += 1
                elif len(r.split("\t")[1]) < 500:
                    len_300_500 += 1
                elif len(r.split("\t")[1]) < 1000:
                    len_500_1000 += 1
                elif len(r.split("\t")[1]) < 2000:
                    len_1000_2000 += 1
                elif len(r.split("\t")[1]) < 5000:
                    len_2000_5000 += 1
                else:
                    len_5000 += 1
except:
    pass

print(len_300)
print(sv_dict)
# print(len(sv_dict))
max_len = 0

pacbio_name = "pacbio_filtered.vcf"

try:
    with open("results/" + dataset_simple + "/" + pacbio_name, "r") as pacbio:
         for r in pacbio.readlines():
             if dataset.startswith("simulated"):
                 break
             #if r.startswith("#") or r.find("DEL") != -1 or r.find("Tandem") != -1 or r.find("Alu") != -1 or r.find("ANN=L1") != -1:
             #    continue
             if r.startswith("#"):
                 continue

             splitted = r.split("\t")

             length = len(splitted[-1])
             chrom = splitted[0]
             if dataset.startswith("HG002"):
                chrom = "chr" + chrom
                length = len(splitted[4])
             if dataset == "NA19240" or dataset.startswith("chm") :
                length = get_len(r)

             if pacbio_name == "t2t.vcf":
                 length = len(splitted[4])

             sv = SV(chrom, int(splitted[1]), length) #get_len(r)
             if length > max_len:
                 max_len = length
             if length < 300:
                 pass
             elif length < 300:
                 len_50_300 += 1
             elif length < 500:
                 len_300_500 += 1
             elif length < 1000:
                 len_500_1000 += 1
             elif length < 2000:
                 len_1000_2000 += 1
             else:
                 len_5000 += 1


             if chrom not in sv_dict:
                 sv_dict[chrom] = []
             sv_dict[chrom].append(sv)
except:
    pass

print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000))
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))
print("Max len " + str(max_len))

pamir_dict = {}
near = 0
not_near = 0

# print(len(sv_dict))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0

total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0
max_len = 0
ins = 0
identity_vector = []
try:
    with open("results/" + dataset + "/pamir.vcf", "r") as pamir_vcf:
        for r in pamir_vcf.readlines():
            #break
            if r.startswith("#"):
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, get_len(r), get_seq(r))
            #new_sv.length = len(r.split("\t")[4])
#            new_sv = SV(chrom, pos, get_len(r), r.split("\t")[4])
            if new_sv.length > max_len:
                max_len = new_sv.length
            if (new_sv.length < 300):
                continue
            ins += 1
            found = False
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if chrom not in pamir_dict:
                pamir_dict[chrom] = []
            pamir_dict[chrom].append(new_sv)

            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if Near(sv, new_sv):

                    if new_sv.length >= 50  and abs(sv.length - new_sv.length) <= 0.05 * sv.length:
                        near += 1
                        sv.checked = True
                        if new_sv.length < 50:
                            pass
                        elif new_sv.length < 300:
                            len_50_300 += 1
                        elif new_sv.length < 500:
                            len_300_500 += 1
                        elif new_sv.length < 1000:
                            len_500_1000 += 1
                        elif new_sv.length < 2000:
                            len_1000_2000 += 1
                        elif new_sv.length < 5000:
                            len_2000_5000 += 1
                        else:
                            len_5000 += 1

                        #identity_vector.append(align_sequences(sv.seq, new_sv.seq))
                        seq = r.split("\t")[7].split(";")[5][4:]
                        #print(seq)
                        #print(chrom)
                        #print(pos)
                        break
except:
    pass
for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False

identity_vector = []
print("Incorrect sequence - " + str(sum(i < 0.8 for i in identity_vector)))

print("pamir")
print(near)
print(ins)


print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000))
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))
print("Max insertion length " + str(max_len))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0

total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0



manta_dict = {}
near = 0
not_near = 0
ins = 0
# print(len(sv_dict))
max_len = 0

try:
    with open("results/" + dataset + "/manta.vcf", "r") as manta_vcf:
        for r in manta_vcf.readlines():
            #break
            if r.startswith("#"):
                continue
            if r.find("MantaBND") == -1:
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, 500, get_seq(r))

            #if get_len(r) > max_len:
            #    max_len = get_len(r)
            #if get_len(r) < 300:
            #    continue

            found = False
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if chrom not in manta_dict:
                manta_dict[chrom] = []
            manta_dict[chrom].append(new_sv)
            ins += 1
            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if Near(sv, new_sv):

                    #if get_len(r) >= 50:#  and abs(sv.length - new_sv.length) <= 0.05 * sv.length:
                        near += 1
                        sv.checked = True
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

                        #align_sequences(sv.seq, new_sv.seq)
                        seq = r.split("\t")[7].split(";")[5][4:]
                        #print(seq)
                        #print(chrom)
                        #print(pos)
                        break
except:
    pass
for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False


print("manta")
print(ins)
print(near)


print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))
print("Max insertion length " + str(max_len))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0

total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0
max_len = 0
ins = 0

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

popins_dict = {}
near = 0
not_near = 0
ins = 0
max_len = 0

identity_vector = []
threshold_size = 300

try:

    with open("results/" + dataset + "/popins2_extended.vcf", "r") as popins_vcf:
        for r in popins_vcf.readlines():
            #break
            if r.startswith("#"):
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, len(r.strip().split("\t")[-1]), r.strip().split("\t")[-1])
            found = False
            if new_sv.length >= threshold_size:
                ins += 1
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if chrom not in popins_dict:
                popins_dict[chrom] = []
            popins_dict[chrom].append(new_sv)
            if new_sv.length < 50:
                pass
            elif new_sv.length < 300:
                total_len_50_300 += 1
            elif new_sv.length < 500:
                total_len_300_500 += 1
            elif new_sv.length < 1000:
                total_len_500_1000 += 1
            elif new_sv.length < 2000:
                total_len_1000_2000 += 1
            elif new_sv.length < 5000:
                total_len_2000_5000 += 1
            else:
                total_len_5000 += 1
            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if Near(sv, new_sv):

                    if new_sv.length >= threshold_size and abs(sv.length - new_sv.length) <= 0.05 * sv.length:
                        near += 1
                        if new_sv.length > max_len:
                            max_len = new_sv.length

                        if new_sv.length < 50:
                            pass
                        elif new_sv.length < 300:
                            len_50_300 += 1
                        elif new_sv.length < 500:
                            len_300_500 += 1
                        elif new_sv.length < 1000:
                            len_500_1000 += 1
                        elif new_sv.length < 2000:
                            len_1000_2000 += 1
                        elif new_sv.length < 5000:
                            len_2000_5000 += 1
                        else:
                            len_5000 += 1

                        #seq = r.split("\t")[7].split(";")[5][4:]
                        #print(seq)
                        sv.checked = True
                        #identity_vector.append(align_sequences(sv.seq, new_sv.seq))

                        #print(chrom)
                        #print(pos)
                        break
                    else:
                        print(sv.length - new_sv.length)
except:
    print("Hewe")
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
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))
print("Max insertion length " + str(max_len))

print("Total: 50-300 " + str(total_len_50_300) )
print("Total: 300-500 " + str(total_len_300_500))
print("Total: 500-1000 " + str(total_len_500_1000) )
print("Total: 1000-2000 " + str(total_len_1000_2000) )
print("Total: 2000-5000 " + str(total_len_2000_5000))
print("Total: >=5000 " + str(total_len_5000))

print("Incorrect sequence - " + str(sum(i < 0.8 for i in identity_vector)))

identity_vector = []

len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0

total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0
max_len = 0

def get_len_nui(line):
    return int(line[4])

identity_vector = []
near = 0
ins = 0
try:
    with open("results/" + dataset + "/NUI.txt", "r") as nui:
        for r in nui.readlines():
            if r.startswith("ref_chr"):
                continue
            splitted = r.split("\t")
            length = get_len_nui(splitted)
            if length < 300:
                continue
            chrom = splitted[0]
            if dataset.startswith("NA19240"):
                chrom = "chr"+ chrom
            new_sv = SV(chrom, int(splitted[1]), get_len_nui(splitted), splitted[-1])
            ins+=1
            pos = int(splitted[1])
            if get_len_nui(splitted) < 300:
                total_len_50_300 += 1
            elif get_len_nui(splitted) < 500:
                total_len_300_500 += 1
            elif get_len_nui(splitted) < 1000:
                total_len_500_1000 += 1
            elif get_len_nui(splitted) < 2000:
                total_len_1000_2000 += 1
            elif get_len_nui(splitted) < 5000:
                total_len_2000_5000 += 1
            else:
                total_len_5000 += 1


            for sv in sv_dict[chrom]:
                if Near(sv, new_sv) and abs(sv.length - new_sv.length) <= 0.05 * sv.length: #and abs(sv.length - new_sv.length) < 10:
                    if get_len_nui(splitted) > max_len:
                        max_len = get_len_nui(splitted)

                    if sv.checked:
                        continue
                    near += 1
                    if get_len_nui(splitted) < 300:
                        len_50_300 += 1
                    elif get_len_nui(splitted) < 500:
                        len_300_500 += 1
                    elif get_len_nui(splitted) < 1000:
                        len_500_1000 += 1
                    elif get_len_nui(splitted) < 2000:
                        len_1000_2000 += 1
                    elif get_len_nui(splitted) < 5000:
                        len_2000_5000 += 1
                    else:
                        len_5000 += 1
                        print(chrom)
                        print(pos)
                    #identity_vector.append(align_sequences(sv.seq, new_sv.seq))
                    sv.checked = True
                    break
except:
    pass
print("NUI")
print("Total " + str(ins))
print("Shared " + str(near))
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))

print("Total: 50-300 " + str(total_len_50_300) )
print("Total: 300-500 " + str(total_len_300_500))
print("Total: 500-1000 " + str(total_len_500_1000) )
print("Total: 1000-2000 " + str(total_len_1000_2000) )
print("Total: 2000-5000 " + str(total_len_2000_5000))
print("Total: >=5000 " + str(total_len_5000))

print("Max insertion length " + str(max_len))


print("Incorrect sequence - " + str(sum(i < 0.8 for i in identity_vector)))
identity_vector = []

len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0

total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0

for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False


total = 0
near = 0
not_near = 0
num = 0
ins = 0
max_len = 0
lengths = []
anchors = []
try:
    with open("results/" + dataset + "/paftools_filtered.vcf", "r") as my_vcf:
        for r in my_vcf.readlines():
            #break
            if r.startswith("#"):
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, len(r.split("\t")[4]), r.split("\t")[4])

            if len(chrom) > 10:
                continue
            total += len(r.split("\t")[4])
            if len(r.split("\t")[4]) > max_len:
                max_len = len(r.split("\t")[4])
            num += 1
            found = False
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if len(r.split("\t")[4]) >= 300:
                ins += 1
            else:
                continue
            if len(r.split("\t")[4]) < 300:
                total_len_50_300 += 1
            elif len(r.split("\t")[4]) < 500:
                total_len_300_500 += 1
            elif len(r.split("\t")[4]) < 1000:
                total_len_500_1000 += 1
            elif len(r.split("\t")[4]) < 2000:
                total_len_1000_2000 += 1
            elif len(r.split("\t")[4]) < 5000:
                total_len_2000_5000 += 1

            else:
                total_len_5000 += 1
                print(chrom + " " + str(pos))

            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if Near(sv, new_sv):
                    if len(r.split("\t")[4]) >=  300 and abs(sv.length - new_sv.length) <= 0.05 * sv.length:
                        found = True
                        near += 1

                        if len(r.split("\t")[4]) < 300:
                            len_50_300 += 1
                        elif len(r.split("\t")[4]) < 500:
                            len_300_500 += 1
                        elif len(r.split("\t")[4]) < 1000:
                            len_500_1000 += 1
                        elif len(r.split("\t")[4]) < 2000:
                            len_1000_2000 += 1
                        elif len(r.split("\t")[4]) < 5000:
                            len_2000_5000 += 1

                        else:
                            len_5000 += 1
                            print(chrom + " " + str(pos))
                        #align_sequences(sv.seq, new_sv.seq)
                        #seq = r.split("\t")[7].split(";")[5][4:]
                        #print(seq)
                        sv.checked = True
                        #print(chrom)
                        #print(pos)
                        break
                    else:
                        print(abs(sv.length - new_sv.length))
except:
    pass

print(ins)

print("Paftools")
print("Total - " + str(ins))
print("Shared - " + str(near))
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))
print("Max insertion length " + str(max_len))


print("Total: 50-300 " + str(total_len_50_300) )
print("Total: 300-500 " + str(total_len_300_500))
print("Total: 500-1000 " + str(total_len_500_1000) )
print("Total: 1000-2000 " + str(total_len_1000_2000) )
print("Total: 2000-5000 " + str(total_len_2000_5000))
print("Total: >=5000 " + str(total_len_5000))


len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000_5000 = 0
len_5000 = 0


total_len_50_300 = 0
total_len_300_500 = 0
total_len_500_1000 = 0
total_len_1000_2000 = 0
total_len_2000_5000 = 0
total_len_5000 = 0

for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False


total = 0
near = 0
not_near = 0
num = 0
ins = 0
max_len = 0
lengths = []
anchors = []
identity_vector = []
try:
    with open("results/" + dataset + "/novelx.vcf", "r") as my_vcf:
        for r in my_vcf.readlines():
            #break
            if r.startswith("#"):
                continue
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, len(r.split("\t")[4]), r.split("\t")[4])


            if len(chrom) > 10:
                continue
            total += len(r.split("\t")[4])
            num += 1
            found = False
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            if int(r.split("\t")[9]) == 0 or int(r.split("\t")[10]) == 0:
                continue

            if len(r.split("\t")[4]) < 300:
                total_len_50_300 += 1

            elif len(r.split("\t")[4]) < 500:
                total_len_300_500 += 1
            # print(new_sv)
            #                            print(sv)
            elif len(r.split("\t")[4]) < 1000:
                total_len_500_1000 += 1
            elif len(r.split("\t")[4]) < 2000:
                total_len_1000_2000 += 1
            elif len(r.split("\t")[4]) < 5000:
                total_len_2000_5000 += 1
            else:
                total_len_5000 += 1

            if len(r.split("\t")[4]) >= 300:
                ins += 1
                lengths.append(len(r.split("\t")[4]))
                anchors.append(int(r.split("\t")[9]) + int(r.split("\t")[10]))
            else:
                continue
            for sv in sv_dict[chrom]:
                if sv.checked:
                    continue
                if Near(sv, new_sv):
                    if len(r.split("\t")[4]) >=  50 and abs(sv.length - new_sv.length) <= 0.05 * sv.length:
                        if len(r.split("\t")[4]) > max_len:
                            max_len = len(r.split("\t")[4])
                        found = True
                        near += 1
                        if len(r.split("\t")[4]) < 300:
                            len_50_300 += 1

                        elif len(r.split("\t")[4]) < 500:
                            len_300_500 += 1
                        #    print(new_sv)
                        #                            print(sv)
                        elif len(r.split("\t")[4]) < 1000:
                            len_500_1000 += 1
                            print(new_sv)
                        elif len(r.split("\t")[4]) < 2000:
                            len_1000_2000 += 1
                        elif len(r.split("\t")[4]) < 5000:
                            len_2000_5000 += 1
                        else:
                            len_5000 += 1

                            #print(new_sv)

                        #identity_vector.append(align_sequences(sv.seq, new_sv.seq))
                        #seq = r.split("\t")[7].split(";")[5][4:]
                        #print(seq)
                        sv.checked = True
                        #print(chrom)
                        #print(pos)
                        break
                    else:
                        print(abs(sv.length - new_sv.length))
                        print(0.05 * sv.length)

            #if not found:
            #    print(r)
except:
    pass


print("Incorrect sequence - " + str(sum(i < 0.8 for i in identity_vector)))

print(ins)

print("Novel-X")
print("Total - " + str(ins))
print("Shared - " + str(near))
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500-1000 " + str(len_500_1000) )
print("1000-2000 " + str(len_1000_2000) )
print("2000-5000 " + str(len_2000_5000))
print(">=5000 " + str(len_5000))


print("Total: 50-300 " + str(total_len_50_300) )
print("Total: 300-500 " + str(total_len_300_500))
print("Total: 500-1000 " + str(total_len_500_1000) )
print("Total: 1000-2000 " + str(total_len_1000_2000) )
print("Total: 2000-5000 " + str(total_len_2000_5000))
print("Total: >=5000 " + str(total_len_5000))

print("Max insertion length " + str(max_len))
print("Mean length - " + str(np.mean(lengths)))
print("Sum length - " + str(sum(lengths)))

print("Mean anchor sum - " + str(np.mean(anchors)))
print("Anchor sd - " + str(np.std(anchors)))
print("Max anchor sum - " + str(max(anchors)))
print("Max ins length - " + str(max(lengths)))



#for sv_vect in sv_dict.values():
#    for sv in sv_vect:
#        if not sv.checked:
#            print(sv)
