# Novel-X supplementary data

This repository contains supplementary data for "Efficient detection and assembly of non-reference DNA sequences
with linked-reads" paper.

**10X_insertions_supplementary.pdf** contains supplementary tables and figures for Novel-X paper. 

### VCF-files

Results folder contains all vcf files that we used for benchmarks that generated Figures 4-5 and Tables S1-S3 from the article.

NA12878, NA19240, CHM1, CHM13, HG002 contains VCFs created using 10X Chromium datasets.

NA12878_tellseq, HG002_tellseq - UST Tell-Seq datasets.

NA12878_stlfr, HG002_stlfr - stLFR datasets.

Folder simulated contains different VCFs:

* No postfix - simple simulated dataset created using LRsim.
* 80, 60, 40, 20 postfix - datasets that were created by downsampling original simulated removing 80/60/40/20 percent of reads
* spades_first, supernova_first, velvel_second, supernova_second - VCFs that were used to create table S2. Spades_first means that in first assembly round we tried to use SPAdes instead of Velvet.

### Important Python scripts

**compare_vcf.py** - scripts that was used to perform event comparison using positions and SV length. In order to run it go to the root repository folder and pass dataset parameter (simulated, HG002, NA12878_tellseq, ...).

### Novel-X flowchart

flowchart.pdf contains a detailed pipeline scheme.


### Read and barcode usage statistics

In read_barcode_statistics.pdf you can find a statistics on read and barcode usage on the key steps of Novel-X for each dataset from the paper.
