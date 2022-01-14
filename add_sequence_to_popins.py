from Bio import SeqIO
import sys

record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

with open(sys.argv[2], 'r') as input_vcf:
    with open(sys.argv[3], 'w') as output_vcf:
        for r in input_vcf.readlines():
            if r.startswith("#"):
                output_vcf.write(r)
                continue
            r = r.split()
            start = r[4].index("p")
            do_rc = False
            try:
                end =  r[4].index("r")
                do_rc = True
            except:
                end =  r[4].index("f")

            contigname = r[4][start:end]
            if not do_rc:
                r.append(str(record_dict[contigname].seq))
            else:
                r.append(str(record_dict[contigname].reverse_complement().seq))

            r = "\t".join(r)
            output_vcf.write(r)
            output_vcf.write("\n")