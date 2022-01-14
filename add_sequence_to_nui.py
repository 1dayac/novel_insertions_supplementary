from Bio import SeqIO
import sys

record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
record_dict2 = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))

with open(sys.argv[3], 'r') as input_nui:
    with open(sys.argv[4], 'w') as output_nui:
        for r in input_nui.readlines():
            if r.startswith("#"):
                continue
            r = r.strip().split()
            contigname = r[3]
            start = int(r[9].split(":")[1].split("-")[0]) - int(r[3].split(":")[1].split("-")[0])
            end = int(r[9].split(":")[1].split("-")[1]) - int(r[3].split(":")[1].split("-")[0])
            if contigname.endswith("_1") or contigname.endswith("_2"):
                contigname = contigname[:-2]
            if contigname in record_dict.keys():
                seq = record_dict[contigname]
                subseq = str(seq.seq[int(start):int(end)])
                if r[5] == "-":
                    subseq = subseq.reverse_complement()

                r.append(subseq)
                output_nui.write("\t".join(r))
                output_nui.write("\n")
            elif contigname in record_dict2.keys():
                seq = record_dict2[contigname]
                subseq = str(seq.seq[int(start):int(end)])
                if r[5] == "-":
                    subseq = subseq.reverse_complement()
                r.append(subseq)
                output_nui.write("\t".join(r))
                output_nui.write("\n")
            else:
                r.append("-")
                output_nui.write("\t".join(r))
                output_nui.write("\n")
                print("Not found")