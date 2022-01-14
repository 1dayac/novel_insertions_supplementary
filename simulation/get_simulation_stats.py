
class SV:
    def __init__(self, chrom, pos, length, seq = ""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.seq = seq.upper()

sv_dict = {}
with open("chm1_novel.txt") as file:
    for line in file.readlines():
        position = int(line.split("\t")[0].split("/")[1])
        chrom = line.split("\t")[0].split("/")[0]
        nucls = line.split("\t")[1].strip()
        length = len(line.split("\t")[1].strip())
        if length >= 0:

            sv = SV(chrom, position, length, nucls)
            if chrom not in sv_dict:
                sv_dict[chrom] = []
            sv_dict[chrom].append(sv)

def Near(sv1, sv2):
    if abs(sv1.pos - sv2.pos) >= 100 and abs(sv1.pos - sv2.pos) < 1000:
        print("Here")
    return abs(sv1.pos - sv2.pos) <= 100

def get_popins_len(line):
    start_pos = line.find("length_")
    ans = ""
    start_pos += 7
    while line[start_pos].isdigit():
        ans += line[start_pos]
        start_pos += 1
    return int(ans) - 55

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
def get_seq(line):
    start_pos = line.find("SEQ")
    ans = ""
    start_pos += 4
    while line[start_pos].isalpha():
        ans += line[start_pos]
        start_pos += 1
    return ans

pamir_dict = {}
near = 0
pamir_total = 0
# print(len(sv_dict))
len_50_300 = 0
len_300_500 = 0
len_500_1000 = 0
len_1000_2000 = 0
len_2000 = 0
try:
    with open("pamir_simulated.vcf", "r") as pamir_vcf:
        for r in pamir_vcf.readlines():
            #break
            if r.startswith("#") or get_len(r) < 50:
                continue
            pamir_total += 1
            chrom = r.split("\t")[0]
            pos = int(r.split("\t")[1])
            new_sv = SV(chrom, pos, get_len(r), get_seq(r))


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
                    if get_len(r) >= 50:
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
print(pamir_total)

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

    with open("popins_simulated.vcf", "r") as popins_vcf:
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
                    elif get_popins_len(r) >= 500:
                        len_500_1000 += 1

                    #seq = r.split("\t")[7].split(";")[5][4:]
                    #print(seq)
                    sv.checked = True
                    print(chrom)
                    print(pos)
                    break
                print("!" + str(chrom) + "_" + str(pos))
except:
    pass

for sv_vect in sv_dict.values():
    for sv in sv_vect:
        sv.checked = False
print("popins")
print("Total (<300) - " + str(len_50_300))
print("Total (<500) - " + str(len_300_500))
print("Total (>500) - " + str(len_500_1000))


near = 0
len_50_300 = 0
len_300_500 = 0
len_500 = 0
with open("novel-x.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, len(r.split("\t")[4]), r.split("\t")[4])
        if chrom not in sv_dict:
            sv_dict[chrom] = []
        for sv in sv_dict[chrom]:
            if sv.checked:
                pass
                continue
            if Near(sv, new_sv):

                if len(r.split("\t")[4]) < 300:
                    len_50_300 += 1
                elif len(r.split("\t")[4]) < 500:
                    len_300_500 += 1
                else:
                    len_500 += 1

                found = True
                near += 1
                #align_sequences(sv.seq, new_sv.seq)
                # seq = r.split("\t")[7].split(";")[5][4:]
                # print(seq)
                sv.checked = True
                print(chrom)
                print(pos)
        if not found:
            print("Here")


total = sum([len(a) for a in sv_dict.values()])

total_50 = 0
total_50_300 = 0
total_300_500 = 0
total_500 = 0

for sv_arr in sv_dict.values():
    for sv in sv_arr:
        if sv.length < 50:
            total_50 += 1
        elif sv.length < 300:
            total_50_300 += 1
        elif sv.length < 500:
            total_300_500 += 1
        else:
            total_500 += 1



print("Novel-X")
print("Total - " + str(total))
print("Total (<50) - " + str(total_50))
print("Total (<300) - " + str(total_50_300))
print("Total (<500) - " + str(total_300_500))
print("Total (>500) - " + str(total_500))

print("Shared - " + str(near))
print("50-300 " + str(len_50_300) )
print("300-500 " + str(len_300_500))
print("500 " + str(len_500) )
