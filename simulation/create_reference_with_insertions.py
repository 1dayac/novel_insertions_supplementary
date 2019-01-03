insertions = []
with open("chm1_novel.txt") as file:
    for line in file.readlines():
        position = int(line.split("\t")[0].split("/")[1])
        chrom = line.split("\t")[0].split("/")[0]
        nucls = line.split("\t")[1].strip()
        insertions.append((chrom, position, nucls))

insertions.sort()

total_nucls = {}

from Bio import SeqIO, Seq

records_dict = SeqIO.to_dict(SeqIO.parse("hg38.fa", "fasta"))

with open("hg.withoutadditional.fa", "w") as output_handle:
    for record in records_dict.values():
        if len(record.id) < 10:
            SeqIO.write(record, output_handle, "fasta")

#records_dict = SeqIO.to_dict(SeqIO.parse("insertions_with_anchors.fasta", "fasta"))
#print(len(records_dict["chr1"]))
#records_dict["chr1"] = records_dict["chr1"][:100] + Seq.Seq("aaaaaaa") + records_dict["chr1"][100:]
#print(len(records_dict["chr1"]))
#exit()
for insertion in insertions:
    if insertion[0] not in total_nucls:
        total_nucls[insertion[0]] = 0
    print("Old " + str(len(records_dict[insertion[0]])))
    new_seq = records_dict[insertion[0]][ : insertion[1] + total_nucls[insertion[0]]] + Seq.Seq(insertion[2]) + records_dict[insertion[0]].seq[insertion[1] + total_nucls[insertion[0]]:]
    print("New " + str(len(new_seq)))

    records_dict[insertion[0]] = new_seq
    total_nucls[insertion[0]] += len(insertion[2])
    print(insertion[2])
    print(total_nucls)


with open("hg.extended.fa", "w") as output_handle:
    for record in records_dict.values():
        if len(record.id) < 10:
            SeqIO.write(record, output_handle, "fasta")


