import gffutils

class SV:
    def __init__(self, chrom, pos, length, seq = ""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.seq = seq.upper()

#fn = "../gencode.v24.annotation.gtf"
#db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', disable_infer_transcripts = True, disable_infer_genes=True)

db = gffutils.FeatureDB('test.db', keep_order=True)
overlap = 0
non_overlap = 0
with open("./../results/chm1/novelx.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, len(r.split("\t")[4]), r.split("\t")[4])

        if len(r.split("\t")[4]) >= 300:
            if len(list(db.region(region=(chrom, pos, pos+1)))):
                for f in list(db.region(region=(chrom, pos, pos+1))):
                    print(f)
                overlap += 1
            else:
                non_overlap += 1
print(overlap)
print(non_overlap)
            #    print(feature)