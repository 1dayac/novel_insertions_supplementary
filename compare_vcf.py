import sys

class SV:
    def __init__(self, chrom, pos, length):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False


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
len_300 = 0
with open("./../pacbio.txt", "r") as pacbio:
     for r in pacbio.readlines():
         break
         if r.split("\t")[1].split("/")[0] not in sv_dict:
             sv_dict[r.split("\t")[1].split("/")[0]] = []
         sv = SV(r.split("\t")[1].split("/")[0], int(r.split("\t")[1].split("/")[1]), len(r.split("\t")[2]))
         if len(r.split("\t")[2]) > 300:
            sv_dict[r.split("\t")[1].split("/")[0]].append(sv)
            len_300 += 1

print(len_300)
print(sv_dict)
# print(len(sv_dict))
len_50_300 = 0
len_300_500 = 0
len_500 = 0
with open("CHM1_final_genotypes.annotated.vcf", "r") as pacbio:
     for r in pacbio.readlines():
         #break
         if r.startswith("#") or r.find("deletion") != -1:
             continue
         splitted = r.split("\t")
         sv = SV(splitted[0], int(splitted[1]), get_len(r))
         if get_len(r) < 300:
             pass
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
pamir_dict = {}
near = 0
not_near = 0

len_50_300 = 0
len_300_500 = 0
len_500 = 0
with open("insertions_setcover.vcf", "r") as pamir_vcf:
    for r in pamir_vcf.readlines():
        #break
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, get_len(r))


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
                if get_len(r) < 50:
                    pass
                elif get_len(r) < 300:
                    len_50_300 += 1
                elif get_len(r) < 500:
                    len_300_500 += 1
                else:
                    len_500 += 1

                near += 1
                sv.checked = True
                #seq = r.split("\t")[7].split(";")[5][4:]
                #print(seq)
                print(chrom)
                print(pos)
for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False


print("pamir")
print(near)


print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )
len_50_300 = 0
len_300_500 = 0
len_500 = 0

popins_dict = {}
near = 0
not_near = 0
ins = 0
with open("popins.vcf", "r") as popins_vcf:
    for r in popins_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = "chr" + r.split("\t")[0]
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


        for sv in sv_dict[chrom]:
            if sv.checked:
                continue
            if Near(sv, new_sv):
                near += 1
                if get_popins_len(r) < 50:
                    pass
                elif get_popins_len(r) < 300:
                    len_50_300 += 1
                elif get_popins_len(r) < 500:
                    len_300_500 += 1
                else:
                    len_500 += 1
                #seq = r.split("\t")[7].split(";")[5][4:]
                #print(seq)
                sv.checked = True
                print(chrom)
                print(pos)
                break


for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False
print("popins")
print(near)
print(ins)
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )
len_50_300 = 0
len_300_500 = 0
len_500 = 0
total = 0
near = 0
not_near = 0
num = 0
ins = 0
with open("test.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, len(r.split("\t")[4]))


        total += len(r.split("\t")[4])
        num += 1
        found = False
        if chrom not in sv_dict:
            sv_dict[chrom] = []
        if len(r.split("\t")[4]) >= 300:
            ins += 1

        for sv in sv_dict[chrom]:
            if sv.checked:
                continue
            if Near(sv, new_sv):

                if len(r.split("\t")[4]) < 300:
                    len_50_300 += 1
                elif len(r.split("\t")[4]) < 500:
                    len_300_500 += 1
                else:
                    len_500 += 1

                near += 1
                # seq = r.split("\t")[7].split(";")[5][4:]
                # print(seq)
                sv.checked = True
                print(chrom)
                print(pos)


print(ins)

print("Novel-X")
print("Total - " + str(num))
print("Shared - " + str(near))
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )





