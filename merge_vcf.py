import sys


filename1 = sys.argv[1]
filename2 = sys.argv[2]
out_filename = sys.argv[3]
vcf_dict = {}

with open(filename1, "r") as file1, open(filename2, "r") as file2:
    for line in file1.readlines():
        if line[0] == "#":
            continue
        chrom = line.split("\t")[0]
        pos =  line.split("\t")[1]
        ins_len = len(line.split("\t")[4])
        if ins_len < 50:
            continue
        vcf_dict[chrom + "_" + pos] = line
    for line in file2.readlines():
        if line[0] == "#":
            continue
        chrom = line.split("\t")[0]
        pos =  line.split("\t")[1]
        ins_len = len(line.split("\t")[4])
        if ins_len < 50:
            continue
        if chrom + "_" + pos not in vcf_dict:
            vcf_dict[chrom + "_" + pos] = line

with open(out_filename, "w") as out_file:
    for key, value in vcf_dict.items():
        out_file.write(value)
