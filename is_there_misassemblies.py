class SV:
    def __init__(self, chrom, pos, pos_start, pos_end):
        self.chrom = chrom
        self.pos = pos
        self.pos_start = pos_start
        self.pos_end = pos_end

sv_dict = {}

with open("CHM1_180GB_CrG_GRCh38_phased_possorted.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        start_pos = int(r.split("\t")[11])
        end_pos = int(r.split("\t")[12])
        node = r.split("\t")[8]
        new_sv = SV(chrom, pos, start_pos, end_pos)
        if len(r.split("\t")[4]) < 300:
            continue
        if node not in sv_dict.keys():
            sv_dict[node] = []
        sv_dict[node].append(new_sv)

#print(sv_dict)
non_misassembled1 = {}
non_misassembled2 = {}

with open("CHM1.1.tsv", "r") as alignments:
    lines = alignments.readlines()
    for i in range(len(lines)):

        r = lines[i]

        if r[0].isdigit():
            node = r.split("\t")[5]
            if node not in sv_dict.keys():
                continue
            else:
                for sv in sv_dict[node]:
                    start_pos = min(int(r.split("\t")[2]), int(r.split("\t")[3]))
                    end_pos = max(int(r.split("\t")[2]), int(r.split("\t")[3]))
                    if sv.pos_start > start_pos and  sv.pos_end < end_pos:
                        if node not in non_misassembled1.keys():
                            non_misassembled1[node] = []
                        non_misassembled1[node].append(r)


with open("CHM1.2.tsv", "r") as alignments:
    lines = alignments.readlines()
    for i in range(len(lines)):

        r = lines[i]

        if r[0].isdigit():
            node = r.split("\t")[5]
            if node not in sv_dict.keys():
                continue
            else:
                for sv in sv_dict[node]:
                    start_pos = min(int(r.split("\t")[2]), int(r.split("\t")[3]))
                    end_pos = max(int(r.split("\t")[2]), int(r.split("\t")[3]))
                    if sv.pos_start > start_pos and  sv.pos_end < end_pos:
                        if node not in non_misassembled2.keys():
                            non_misassembled2[node] = []
                        non_misassembled2[node].append(r)

#print(non_misassembled1)
#print(non_misassembled2)

answer = []

for node in non_misassembled1.keys():
    for r in non_misassembled1[node]:
        answer.append(node)

for node in non_misassembled2.keys():
    if node not in answer:
        for r in non_misassembled2[node]:
            answer.append(node)
print(answer)