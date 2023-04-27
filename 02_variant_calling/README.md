# 02_variant_calling

This folder contains scripts and codes which does the following:

1. Remove barcodes from the raw reads
2. Remove adapters from the raw reads
3. Align raw reads to reference genome
4. Perform variant calling using BCFtools
5. Perform variant calling using STACKS
6. Filter the variants called via BCFtools and STACKS

## Remove barcodes from the raw reads
```bash
process_radtags -p /path/to/folder/with/raw/reads -b /path/to/barcodesID.txt -o /path/to/output/folder1 -e pstI -c -q -r
# /path/to/folder/with/raw/reads is the path to the folder containing the raw reads
# /path/to/barcodesID.txt is the path to tab delimited txt file containing the barcode sequence in the first column (e.g. CCGTCTTATGCACT), and the sample ID in the second column (e.g. M0001). Example of the first 5 rows of this file:
# CCGTCTTATGCACT  M0001a
# ACATGACCGACCT   M0002a
# TGTTGCCGCGT     M0003a
# TGCCGACAAGAA    M0004a
# CACCGTCCATG     M0005a
# /path/to/output/folder1 is the path to the folder where the reads with the barcodes removed will be stored

```

Barcodes were also removed from raw reads of the samples from the Australian Museum samples (Ewart et al., 2019) using the process_radtags program in STACKS version 2.58 (Catchen et al., 2013). However, this had to be done on each individual file at a time to avoid mixing issues as some barcodes were subset of other barcodes (e.g. ACATGACCGACCT for sample A, and ACATGACCG for sample B). This means that the process has to be looped over all the raw read file. This was executed using the following command:

```bash
process_radtags -f /path/to/raw/read/file.FASTQ.gz -b /path/to/sample/barcodesID.txt -o /path/to/output/folder1 -e pstI -c -q -r
# /path/to/raw/read/file.FASTQ.gz is the path to the sample raw read file. 
# /path/to/sample/barcodesID.txt is the path to tab delimited txt file containing the barcode sequence in the first column (e.g. CCGTCTTATGCACT), and the sample ID in the second column (e.g. M0001). However, this is only for the particular sample. Example of the one of the file:
# CTGCT	12709_rep_b 
# /path/to/output/folder1 is the path to the folder where the reads with the barcodes removed will be stored

```
## Remove adapters from the raw reads
```bash
fastp -l 40 --adapter_sequence=AGATCGGAAGAG -i /path/to/output/folder1/barcodes_removed.fq.gz -o /path/to/output/folder2/adapters_trimmed.fq.gz
# /path/to/output/folder1/barcodes_removed.fq.gz is the path to the file with the barcodes removed from the previous step. 
# /path/to/output/folder2/adapters_trimmed.fq.gz is the path to the output file (barcodes removed and adapters trimmed).

```
