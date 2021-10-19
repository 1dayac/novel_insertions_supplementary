import sys
class SV:
    def __init__(self, chrom, pos, length, seq=""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.seq = seq.upper()

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

with open(sys.argv[2], "w") as pacbio_fasta:
        with open(sys.argv[1], "r") as pacbio:
            for r in pacbio.readlines():
                #break
                # if r.startswith("#") or r.find("DEL") != -1 or r.find("Tandem") != -1 or r.find("Alu") != -1 or r.find("ANN=L1") != -1:
                #    continue
                if r.startswith("#") or r.find("DEL") != -1:
                    continue

                splitted = r.split("\t")
                sv = SV(splitted[0], int(splitted[1]), len(splitted[4]), splitted[4])
                if len(sv.seq) >= 300:
                    pacbio_fasta.write(">" + splitted[0] + "_" + splitted[1] + "\n")
                    pacbio_fasta.write(sv.seq + "\n")