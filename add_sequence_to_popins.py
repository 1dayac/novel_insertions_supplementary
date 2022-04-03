from Bio import SeqIO
import sys

record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

def parse_contig_line(line):
    r = line.split()
    position = -1
    start = r[4].index("c")
    do_rc = False
    try:
        end = r[4].index("r")
        do_rc = True
    except:
        end = r[4].index("f")

    contigname = r[4][start:end]
    try:
        start_pos = r[4].index(":") + 1
        start_pos2 = start_pos
        while r[4][start_pos].isdigit():
            start_pos += 1
        position = int(r[4][start_pos2 : start_pos])
    except:
        pass
    return contigname, do_rc, position


with open(sys.argv[2], 'r') as input_vcf:
    with open(sys.argv[3], 'w') as output_vcf:
        all_lines = input_vcf.readlines()
        current = 0
        while current < len(all_lines):
            if all_lines[current].startswith("#"):
                output_vcf.write(all_lines[current])
                current += 1
                continue
            if current == len(all_lines) - 1:
                current += 1
                continue

            contigname1, do_rc1, position1 = parse_contig_line(all_lines[current])
            contigname2, do_rc2, position2 = parse_contig_line(all_lines[current + 1])
            if contigname1 == contigname2 and position1 >= 0 and position2 >= 0:
                if position1 > position2:
                    position1, position2 = position2, position1
                r = all_lines[current].split()
                try:
                    if not do_rc1:
                        r.append(str(record_dict[contigname1].seq)[position1 : position2])
                    else:
                        r.append(str(record_dict[contigname1].reverse_complement().seq)[position1 : position2])
                except:
                    if not do_rc1:
                        r.append(str(record_dict[contigname1].seq))
                    else:
                        r.append(str(record_dict[contigname1].reverse_complement().seq))
                r = "\t".join(r)
                output_vcf.write(r)
                output_vcf.write("\n")
            current += 1



