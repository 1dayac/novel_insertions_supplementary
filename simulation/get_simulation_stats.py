
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
