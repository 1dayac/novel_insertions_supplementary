from subprocess import call, Popen

with open("CHM13_180GB_CrG_GRCh38_phased_possorted.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        if len(chrom) > 7:
            continue
        pos = int(r.split("\t")[1])


        number_of_insertion = r.split("\t")[8].split("_")[1]

        start = max(0, pos - 5000)
        end = pos + 500
        with open(number_of_insertion + ".bam", 'w') as bam:
            process = Popen(['samtools', 'view', "-b", 'sample/CHM13_180GB_CrG_GRCh38_phased_possorted.bam', chrom + ":" + str(start) + "-" + str(end)], stdout=bam)
            process.wait()
        with open(number_of_insertion + ".txt", 'w') as stats:
            process = Popen(['/pbtech_mounts/homes048/dmm2017/Novel-X/bxtools/bin/bxtools', 'stats', number_of_insertion + ".bam"], stdout=stats)
            process.wait()

        barcodes = set()
        with open("CHM13_180GB_CrG_GRCh38_phased_possorted_barcodes/" + number_of_insertion + ".txt", "r") as barcode_list:
            for r in barcode_list.readlines():
                barcodes.add(r.strip())

        bam_barcodes = set()
        with open(number_of_insertion + ".txt", "r") as new_barcode_list:
            for r in new_barcode_list.readlines():
                if not r.startswith("N"):
                    bam_barcodes.add(r.split("\t")[0].strip())


        with open("percentage_of_barcodes_in.txt", 'a') as file1:
            with open("percentage_of_barcodes_of_bam.txt", 'a') as file2:

                in_bam = 0
                not_in_bam = 0

                for bx in barcodes:
                    if bx in bam_barcodes:
                        in_bam += 1
                    else:
                        not_in_bam += 1
                try:
                    file1.write(str(float(in_bam)/float(in_bam+not_in_bam)))
                    file1.write("\n")
                    file2.write(str(float(in_bam)/len(bam_barcodes)))
                    file2.write("\n")
                except:
                    print("Zero barcodes")


