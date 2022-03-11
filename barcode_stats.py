import pysam
import sys
import os
import portion as P
bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
try:
    os.makedirs(sys.argv[2])
except:
    pass

iter = bamfile.fetch("chr1", 1000000, 3000000)

max_barcodes = 20000
current_barcodes = 0

superdict = {}

for record in iter:
    try:
        bx_tag = record.get_tag("BX")
    except:
        continue
    if bx_tag not in superdict and current_barcodes < max_barcodes:
        superdict[bx_tag] = []
        current_barcodes += 1
    if bx_tag not in superdict:
        continue
    superdict[bx_tag].append((record.reference_start, record.reference_end))

molecule_length = []
molecule_coverage = []
max_diff = 10000
for barcode, vect in superdict.items():
    if len(vect) <= 2:
        continue
    print(barcode)
    #print(vect)
    aligned = 0
    unaligned = 0
    start = vect[0][0]
    superinterval = P.empty()
    for interval in vect:
        superinterval = superinterval | P.closed(interval[0], interval[1])
    superinterval = list(superinterval)
    if len(superinterval) <= 1:
        continue

    for i in range(len(superinterval) - 1):
        if superinterval[i + 1].lower - superinterval[i].upper > max_diff:
            aligned += superinterval[i].upper - superinterval[i].lower
            if superinterval[i].upper - start > 5000:
                molecule_coverage.append(aligned/float(aligned+unaligned))
            if superinterval[i].upper - start > 5000:
                molecule_length.append(superinterval[i].upper - start)
            start = superinterval[i+1].lower
            aligned = 0
            unaligned = 0
        else:
            aligned += superinterval[i].upper - superinterval[i].lower
            unaligned += superinterval[i+1].lower - superinterval[i].upper
    aligned += superinterval[-1].upper - superinterval[-1].lower
    if superinterval[i].upper - start > 5000:
        molecule_coverage.append(aligned/float(aligned+unaligned))
    if superinterval[i].upper - start > 5000:
        molecule_length.append(superinterval[-1].upper - start)

from matplotlib import pyplot as plt
import numpy as np
bins1 = np.arange(0, 1, 0.01)
bins2 = np.arange(0, 200000, 1000)

print(molecule_length)
print(molecule_coverage)

with open(sys.argv[2] + "/mol_length.csv", 'w') as molecule_length_csv:
    molecule_length_csv.write("length")
    for item in molecule_length:
        molecule_length_csv.write(item+"\n")

with open(sys.argv[2] + "/mol_coverage.csv", 'w') as molecule_length_csv:
    molecule_length_csv.write("coverage")
    for item in molecule_length:
        molecule_length_csv.write(item+"\n")

plt.hist(molecule_length, bins=bins2, alpha=0.5)
plt.title('Molecule length distribution')
plt.xlabel('Molecule length observed')
plt.ylabel('Count')
plt.savefig(sys.argv[2] + "/mol_length.png", dpi = 150)
plt.close()

plt.hist(molecule_coverage, bins=bins1, alpha=0.5)
plt.title('Molecule coverage distribution')
plt.xlabel('Molecule coverage observed')
plt.ylabel('Count')
plt.savefig(sys.argv[2] + "/mol_coverage.png", dpi = 150)
plt.close()

fig, ax = plt.subplots()
ax.scatter(molecule_length,molecule_coverage)
plt.savefig(sys.argv[2] + "/scatterplot.png", dpi = 150)

