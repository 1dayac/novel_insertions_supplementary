import  sys
from subprocess import call, Popen

count_bad = 0
count_good = 0
with open(sys.argv[1]) as vcf:
    for r in vcf.readlines():
        if (r.startswith('#')):
            continue
        splitted = r.split("\t")
        chrom = splitted[0]
        pos = int(splitted[1])
        start = pos - 300
        end = pos + 300
        if len(splitted[4]) < 300:
            continue
        process = Popen(['/local/workdir/dmm2017/novel-x-paper/chm1_mappability/bigWigToWig', '-chrom=' + chrom,
                         '-start=' + str(start), 'end=' + str(end), sys.argv[2], chrom + '_' + str(pos) + '.wig'])
        process.wait()
        bad = False
        with open(chrom + '_' + str(pos) + '.wig') as wig:
            for r in wig.readlines():
                try:
                    f = float(r.strip())
                    if (f < 0.9):
                        bad = True
                except:
                    pass

        if bad:
            count_bad += 1
        else:
            count_good += 1

print("No mappability issues " + str(count_good))
print("Mappability issues " + str(count_bad))
