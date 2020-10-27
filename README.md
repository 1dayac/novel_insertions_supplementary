# Novel-X supplementary data

This repository contains data that didn't show up in "Novel sequence insertion detection using Linked-Reads" paper. It is not structured at the moment but will be.

**10X_insertions_supplementary.pdf** contains supplementary tables and figures for Novel-X paper. 

### Novel-X flowchart

flowchart.pdf contains a detailed pipeline scheme.

### VCF-files

Results folder contains all vcf files that we used for generating Tables from the article.

### Simulated dataset

Folder simulation contains various scripts we used to obtain results for the "Simulated dataset" section:

* ... - script for creating reference with insertions, given list of insertions and reference.
* ... - commands used to simulated reads and bam-file
* ... - list of insertions used (borrowed from Huddlestone et al, 2015)
* ...vcf - list of insertions found

### ALU and tandem repeat filtering

Folder filtered_statistics contains comparison with SMRT-SV calls without ALU and tandem repeat insertions:

* NA19240_table_filtered.pdf - table for NA19240 dataset
* ....py - script for comparison of vcf-files

### Read and barcode usage statistics

In read_barcode_statistics.pdf you can find a statistics on read and barcode usage on the key steps of Novel-X for each dataset from the paper.
