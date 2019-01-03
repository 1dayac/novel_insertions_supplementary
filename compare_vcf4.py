import sys


class SV:
    def __init__(self, chrom, pos, length):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.in_pacbio = False
        self.in_novel = False
        self.in_pamir = False


sv_dict = {}


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


def Near(sv1, sv2):
    return abs(sv1.pos - sv2.pos) <= 100


# with open("pacbio.txt", "r") as pacbio:
#    for r in pacbio.readlines():
#        if r.split("\t")[1].split("/")[0] not in sv_dict:
#            sv_dict[r.split("\t")[1].split("/")[0]] = []
#        sv = SV(r.split("\t")[1].split("/")[0], int(r.split("\t")[1].split("/")[1]), len(r.split("\t")[2]))
#        if len(r.split("\t")[2]) > 500:
#            sv_dict[r.split("\t")[1].split("/")[0]].append(sv)

with open("test.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, len(r.split("\t")[4]))
        if len(r.split("\t")[4]) < 300:
            continue

        if chrom not in sv_dict:
            sv_dict[chrom] = []
        sv_dict[chrom].append(new_sv)



with open("insertions_setcover.vcf", "r") as pamir_vcf:
    for r in pamir_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, get_len(r))
        if get_len(r) < 300:
            continue
        if chrom not in sv_dict:
            sv_dict[chrom] = []
        for sv in sv_dict[chrom]:
            if Near(sv, new_sv):
                sv.in_pamir = True

print(len(sv_dict))
with open("CHM1_final_genotypes.annotated.vcf", "r") as pacbio:
    for r in pacbio.readlines():
        if r.startswith("#") or r.find("deletion") != -1:
            continue
        splitted = r.split("\t")
        new_sv = SV(splitted[0], int(splitted[1]), get_len(r))
        if get_len(r) < 300:
           continue
        chrom = splitted[0]
        if chrom not in sv_dict:
            sv_dict[chrom] = []
        for sv in sv_dict[chrom]:
            if Near(sv, new_sv):
                sv.in_pacbio = True

pamir_dict = {}
near = 0
not_near = 0

print("pamir")
print(near)

len_50_300 = 0
len_300_500 = 0
len_500 = 0

popins_dict = {}
near = 0
not_near = 0

print("popins")
print(near)

print("50-300 " + str(len_50_300))
print("300-500 " + str(len_300_500))
print("500 " + str(len_500))
len_50_300 = 0
len_300_500 = 0
len_500 = 0
near = 0
not_near = 0
with open("remapped_popins.vcf", "r") as popins_vcf:
    for r in popins_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom =  r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, get_popins_len(r))
        found = False
        if get_popins_len(r) < 300:
            continue
        found = False
        if chrom not in sv_dict:
            sv_dict[chrom] = []
        for sv in sv_dict[chrom]:
            if Near(sv, new_sv):
                sv.in_novel = True

                near += 1
                # seq = r.split("\t")[7].split(";")[5][4:]
                # print(seq)
                print(chrom)
                print(pos)
print("me")
print(near)

only_popins = 0
popins_smrt = 0
popins_novel = 0
popins_pamir = 0
popins_novel_pamir = 0
popins_novel_smrt = 0
popins_smrt_pamir = 0
all = 0

for chr in sv_dict.keys():
    for sv in sv_dict[chr]:
        if sv.in_novel == False and sv.in_pacbio == False and sv.in_pamir == False:
            only_popins += 1
        if sv.in_novel == True and sv.in_pacbio == False and sv.in_pamir == False:
            popins_novel += 1
        if sv.in_novel == False and sv.in_pacbio == True and sv.in_pamir == False:
            popins_smrt += 1
        if sv.in_novel == False and sv.in_pacbio == False and sv.in_pamir == True:
            popins_pamir += 1
        if sv.in_novel == True and sv.in_pacbio == True and sv.in_pamir == False:
            popins_novel_smrt += 1
        if sv.in_novel == True and sv.in_pacbio == False and sv.in_pamir == True:
            popins_novel_pamir += 1
        if sv.in_novel == False and sv.in_pacbio == True and sv.in_pamir == True:
            popins_smrt_pamir += 1
        if sv.in_novel == True and sv.in_pacbio == True and sv.in_pamir == True:
            all += 1

print(only_popins)
print(popins_smrt)
print(popins_novel)
print(popins_pamir)
print(popins_novel_pamir)
print(popins_novel_smrt)
print(popins_smrt_pamir)
print(all)
