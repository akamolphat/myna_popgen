# 02_variant_calling

This folder contains scripts and codes which does the following:

1. Remove barcodes from the raw reads
2. Remove adapters from the raw reads
3. Align raw reads to reference genome
4. Perform variant calling using BCFtools
5. Perform variant calling using STACKS
6. Filter the variants called via BCFtools and STACKS

## 1. Remove barcodes from the raw reads

The removal of barcodes and adapters, and the alignment of reads to the reference genome was performed for both the STACKS and the BCFtools pipeline. Barcodes were removed from raw reads of the samples we sent to dart using the process_radtags program in STACKS version 2.58 (Catchen et al., 2013), using the following command:

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
## 2. Remove adapters from the raw reads
Adapters are trimmed using fastp version 0.20.0 (Chen et al., 2018) with the following command. The command has to be looped over all files (barcodes removed from previous step):
```bash
fastp -l 40 --adapter_sequence=AGATCGGAAGAG -i /path/to/output/folder1/barcodes_removed.fq.gz -o /path/to/output/folder2/adapters_trimmed.fq.gz
# /path/to/output/folder1/barcodes_removed.fq.gz is the path to the file with the barcodes removed from the previous step. 
# /path/to/output/folder2/adapters_trimmed.fq.gz is the path to the output file (barcodes removed and adapters trimmed).

```

### Quality control checks
Quality control checks were performed using fastqc version 0.11.9 (Andrews, 2010) and multiqc version 1.9 (Ewels et al., 2016). fastqc was performed on each file as output from the previous step while multiqc was performed on the outputs from fastqc using the following commands:

```bash
fastqc /path/to/output/folder2/adapters_trimmed.fq.gz -o /path/to/fastqc/folder 
multiqc /path/to/fastqc/folder -o /path/to/multiqc/folder
# /path/to/fastqc/folder is the path to the folder containing the fastqc outputs
# /path/to/multiqc/folder is the path to the folder for the multiqc outputs
```

## 3. Align raw reads to reference genome
Genomic DNA was extracted from a male myna collected in Newcastle, Australia using the Qiagen DNeasy blood and tissue kit.  A Chromium Genome library v2 was prepared by Macrogen and sequenced in PE150 mode on a NovaSeq. A total of 110.7 Gb data was generated. A genome assembly was generated from a subsample of 260 million read pairs, targeting ~60x coverage using Supernova version 2.1.1.  The resulting assembly was 1 044,988,094 bp in length, comprising 17,560 contigs with a contig N50 of 3.86 Mb. 
Reads were then aligned to the draft common myna reference genome using BWA version 0.7.17 (Li and Durbin, 2009). This was done using the following command which has to be looped over all files.

```bash
bwa mem -R $(echo "@RG\tID:${i%.fq.gz}\tSM:${i%.fq.gz}\tLB:${i%.fq.gz}\tPL:ILLUMINA") $ref $i > /path/to/output/folder3/${i%.fq.gz}.sam
# This executed in the folder containing the files containing the reads with the barcodes removed and adapters trimmed.
# $i is the name of the individual file containing the barcodes removed and adapters trimmed reads.
# /path/to/output/folder3/ is the path to the output folder where the file with reads aligned to the reference would be stored.
```
### Sort and index files
The reads were then sorted and indexed using SAMTOOLS version 1.12 (Li et al., 2009). This was done using the following command which has to be looped over all files.
```
samtools sort ${i%.fq.gz}.sam > ${i%.fq.gz}.sorted.bam
samtools index ${i%.fq.gz}.sorted.bam
samtools flagstat ${i%.fq.gz}.sorted.bam > ${i%.fq.gz}.flagstat.txt
# As above, $i is the name of the individual file containing the barcodes removed and adapters trimmed reads. These commands are to be executed at /path/to/output/folder3/
```
Variants were called based on three subsets of individuals – the entire dataset (ALL dataset: 814 individuals, 1042 samples), a New Zealand-only dataset (NZ dataset: 226 individuals, 282 samples) and an Indian-only dataset (IND dataset: 78 individuals, 80 samples). Note that the number of samples and individuals are different due to the inclusion/removal of replicates, including DArT-specific proprietary replicates which includes approximately 30% of the samples that were sequenced.

## 4. Variant calling with BCFtools
Variants were called through the mpileup and call commands in BCFtools version 1.13 (Li, 2011). This refers to step 4b on the Figure S2.1, and section 2.2 in the main text and on Figure 2.

### bcftools mpileup
```bash
bcftools mpileup -b list_of_file_name.txt -a "DP,AD,SP" --ignore-RG -Ou -f path/to/ref.fasta -o /path/to/mpileup_output.bcf
# list_of_file_name.txt is the path to file containing all the files with the aligned reads from step 3. Ideally this is executed in the directory contain the file and list_of_file_name.txt will only contain the file name
# path/to/ref.fasta is the path to the file with reference myna genome
# /path/to/mpileup_output.bcf is the path to the output .bcf file from bcftools mpileup
```
### bcftools call
bcftools call can be executed on its own but the commands below rename the samples to remove file suffix before performing bcftools call.

```bash
# Make list of samples for renaming (remove .filtered.sorted.bam file suffix)
bcftools query -l /path/to/mpileup_output.bcf | awk -F'.filtered.sorted.bam' '{print $1}' > path/to/renamed_sample.txt
# Rename samples and pipe this into bcftools call
bcftools reheader -s path/to/renamed_sample.txt /path/to/mpileup_output.bcf | bcftools call -mv -Ob -o /path/to/variant_call_output.bcf -f GQ
# path/to/renamed_sample.txt is the path to the file with the sample name with the file suffix removed
# /path/to/variant_call_output.bcf is the path to the .bcf file containing all the variants as identified by bcftools call
```

## 5. Variant calling with STACKS
Variants were called through the gstacks and populations program in STACKS version 2.58 (Catchen et al., 2013). This refers to step 4a on the Figure S2.1. All SNPs identified in gstacks and populations were first output to a VCF file before further filtering with BCFtools version 1.13 (Li, 2011), and VCFtools version 0.1.15 (Danecek et al., 2011). 
### gstacks
gstacks were executed using the reference-based mode as follows:
```bash
gstacks -I /path/to/output/folder3/ -S .filtered.sorted.bam -O /path/to/gstacks/folder/ -M /path/to/popmap_file.txt
# /path/to/output/folder3/ is the path to the folder containing .filtered.sorted.bam files which were output from step 3 (Align to reference)
# /path/to/gstacks/folder/ is the path to the folder where the outputs from gstacks will be stored
# /path/to/popmap_file.txt is the path to the population map giving the list of samples in the first column, and the population number code in the second column. The population number code can just be one number just to allow the program to work. Example of the first 5 rows of this file:
# 10201_rep_b	1
# 10201	1
# 10202_rep_b	1
# 10202	1
# 10203	1
```
### populations
The populations program reads in the outputs from gstacks and export the results in the specified format. The command below was used to export the variants identified to a VCF format:
```bash
populations -P /path/to/gstacks/folder/ -O /path/to/populations/folder/ --vcf
# /path/to/gstacks/folder/ is the path to the folder where the outputs from gstacks are stored
# /path/to/populations/folder/ is the path to the folder where the outputs from populations will be stored
```
### sort and index VCF files
The commands below were used to sort and index the VCF file for further use. Sorting of the VCF was done using the vcf-sort function in VCFtools version 0.1.15 (Danecek et al., 2011):
```bash
cd path/to/populations/folder/
vcf-sort populations.snps.vcf | bgzip -c > populations.snps.sorted.vcf.gz
bcftools index populations.snps.sorted.vcf.gz
# populations.snps.vcf is the output from the populations program above. This is the standard output name from the populations program.
# populations.snps.sorted.vcf is the output is the indexed, sorted vcf file.

```
## 6. Filter the variants called via BCFtools and STACKS

### Biallelic SNPs
The variants identified by the BCFtools and STACKS pipelines were filtered to only retain biallelic SNPs using the view and filter command in BCFtools. After filtering, the resulting file is indexed using the index command in BCFtools. Indexing is performed following all filtering steps.
bcftools view and filter
```bash
bcftools view -m2 -M2 -v snps /path/to/variant_file | bcftools view -c1 | bcftools filter -e 'AC==0 || AC=AN' | bgzip -c > /path/to/vcf_filtered/variant_call_biallelic_snps.vcf.gz
bcftools index /path/to/vcf_filtered/variant_call_biallelic_snps.vcf.gz
# /path/to/variant_file is the path to the .bcf or .vcf.gz file containing all the variants as identified by the BCFtools and STACKS pipelines (section S2.2 and S2.3 above respectively)
# /path/to/vcf_filtered/variant_call_biallelic_snps.vcf.gz is the path to the output vcf.gz file containing only biallelic SNPs
# bcftools view -m2 -M2 -v snps filters the .bcf file for only biallelic SNPs
# bcftools view -c1 removes monomorphic 0/0 loci. This part is actually redundant as this is covered in the next section of the code
# bcftools filter -e 'AC==0 || AC=AN' filters out all monomorphic loci for reference (AC==0) and alternative genotype (AC=AN).
# bgzip -c compresses the output VCF from bcftools filter to a vcf.gz file. Note that this was used over the -Oz -o flag in bcftools because -Oz -o flag requires the input file to be indexed first.
# bcftools index indexed the resulting vcf.gz file which is required for further downstream processes.
```

### Read depths and genotype quality
Generally, minimum read depth filtering will remove false positive genotype calls, while maximum read depth cutoffs will remove genotype calls that may have potentially been mapped onto paralogs or repetitive parts of the genome.

The threshold for the read depths for the STACKS and BCFtools dataset were based on the distribution of depths of well-defined SNPs. Well-defined SNPs are based on the datasets that were filtered to only retain SNPs with GQ ≥ 30 (and QUAL ≥ 30 for BCFtools) with no missing data (genotyped in all individuals) and the distribution of the read depths were determined in each dataset. Whole number thresholds were guided by the 95-percentile centered around the median: minimum depth = 15, and maximum depth = 100 for the NZ and IND datasets and 125 for the ALL dataset. Figure S3.1-3 show the DP distribution for the ALL, NZ, and IND BCFtools datasets respectively. A similar distribution is observed for the STACKS dataset (results not shown). Note that DArT SNP dataset does not have genotype read depth information.

VCFtools was used to recode genotype calls to NA where genotype quality scores (GQ) were < 30 or where the coverage fell outside the chosen read depth cut off point of minimum and maximum read depth thresholds. The BCFtools dataset was also filtered to only retain SNPs with SNP quality scores (QUAL) ≥ 30. After some genotype calls were recoded to NA by VCFtools, the filter command in BCFtools was also used to remove any monomorphic SNPs.
```bash
vcftools --gzvcf ${inputa} --remove-indels --minQ 30 --minGQ 30 --minDP 15 --maxDP ${maxDP} --recode --recode-INFO-all --stdout | bcftools filter -e 'AC==0 || AC=AN' -Oz -o ${outputa}
# ${inputa} is the input vcf.gz file with only biallelic SNPs from the biallelic SNP filter
# ${maxDP} is the maximum depth cutoff point. This is 100 for the IND and NZ dataset, and 125 for the entire/ALL dataset
# ${outputa} is the output vcf.gz file from this filter. 
# NOTE: --minQ 30 is only used when working with BCFtools dataset as the STACKS datasets do not have SNP quality scores
```

### Admixed Australian samples and replicates
Admixed Australian samples were removed from the dataset as admixed individuals are known to affect some downstream analysis. We also remove replicates at random prior to further SNP filtering. Both these steps were executed using the same command, except for the input text file containing the list of samples to be removed (${removesampletxt} in the command below):
```bash
bcftools view -S ^${removesampletxt} --force-samples ${inputa} | bcftools filter -e 'AC==0 || AC=AN' -Oz -o ${outputa}
# ${removesampletxt} is the path to the .txt file containing the sample ID(s) of the sample(s) to be removed from the dataset. This could be the list of admixed Australian samples or the list of replicates.
# ${inputa} is the resulting vcf.gz file from section S3.2 above.
# ${outputa} is the output vcf.gz file after the removal of the samples in ${removesampletxt}.
```
### Inbreeding coefficient (FIS)
FIS is the proportion of the variance in the subpopulation contained in an individual. A high positive value indicates some levels of inbreeding, and a high negative value indicates heterozygote excess (outbreeding, or potential sample mixing). Figure S3.4 shows the FIS for the BCFtools ALL dataset as calculated by the vcftools function with the --het option in VCFtools, showing the three anomalous samples (sample 12718, and 12719 from Sydney, and M0271 from India) with FIS < -0.25. FIS outliers were also observed in the other datasets in a similar manner (IND dataset, and STACKS dataset, results not shown). In the DART dataset, the FIS outliers were also observed with observed heterozygosity (Ho), but also including sample M0208 from Fiji (results not shown).
The anomalous samples were removed using the view command in BCFtools. This was only executed in the ALL and the IND dataset and was executed alongside the SNP missingness filter (see command in section below).

### SNP missingness
The locus data missingness thresholds were determined based on a balance between the number of SNPs retained and the data missingness. A histogram of the number of SNPs and the data missingness was produced and an arbitrary threshold was used. Figure S3.6 shows this for the ALL dataset as calculated by the vcftools function with the --missing-site option in VCFtools. Similar distributions are observed with the NZ and IND dataset, and also datasets from the STACKS and DART pipelines but not shown here.

A threshold of 20% SNP data missingness was chosen and the view command in BCFtools was used to remove SNPs with more than 20% missingness. This was executed alongside the removal of the samples with anomalous inbreeding coefficients. 
```bash
bcftools view -s ^12718,12719,M0271 ${inputa} | bcftools view -i 'F_MISSING<=0.2' | bcftools filter -e 'AC==0 || AC=AN' -Oz -o ${outputa}
# or
bcftools view -i 'F_MISSING<=0.2' ${inputa} | bcftools filter -e 'AC==0 || AC=AN' -Oz -o ${outputa}
# The first command excludes sample 12718, 12719 and M0271. The list of samples to exclude may change depending on the dataset being filtered.
# The second command does not exclude any samples and is used for the NZ dataset which has no samples with anomalous inbreeding coefficient.
# ${inputa} is the resulting vcf.gz file from section 3.3 above
# ${outputa} is the output vcf.gz file after the removal of samples with anomalous inbreeding coefficients, the removal of SNPs with more than 20% missingness, and the removal of monomorphic SNPs.
```
###	Singletons and doubletons

For all downstream analyses, excluding the construction of the site-frequency-spectrum (SFS), demographic inference, and gene flow analysis (see Figure 2 in the main text and Figure S2.1), singletons and doubletons (SNPs only occurring in one sample) are removed from the dataset using the following VCFtools commands.
```bash
vcftools --gzvcf ${vcfa} --singletons --out ${singleton}
vcftools --gzvcf ${vcfa} --exclude-positions ${singleton}.singletons --recode --recode-INFO-all --stdout | bcftools annotate --set-id "%CHROM\_%POS\_%REF\_%ALT" -O z -o ${vcfb}
# The first command makes a list of singletons and doubletons
# The second command exclude singletons and doubletons from the VCF file using the list created from the first command and pipes this through bcftools annotate and setting the ID to the CHROM_POS_REF_ALT nomenclature before outputting to a vcf.gz file.
# ${vcfa} is the path to the resulting vcf.gz file from section 3.5 (after removing samples with anomalous inbreeding coefficients and after removing SNPs with more than 20 percent data missingness)
# ${singleton} is path to the output text file from the first command, containing a list of singletons and doubletons.
# ${vcfb} is the path to the output vcf.gz file after the removal of singletons and doubletons.
```

### Linkage disequilibrium and SNP thinning

Pairwise LD between all loci within 500kb on the same contig was calculated and plotted against distance in kb using the PopLDdecay program (Zhang et al., 2019). This was performed prior to thinning the SNPs in the VCF files using the following command:
```bash
./PopLDdecay -InVCF ${vcfa} -OutStat ${outpath}
perl Plot_OnePop.pl -inFile ${outpath}.stat.gz -output ${outputPDF}
# The first command runs PopLDdecay and output a file containing all pairwise LD between all loci within 500kb on the same contig
# The second command makes the LD decay plots (e.g. Figure S4.1-3) based on the output from the first command using a perl script that is provided with PopLDdecay program
# ${vcfa} is the path to the vcf.gz file to calculate LD decay. The command was run with vcf.gz files before and after the removal of singletons and doubletons.
# ${outpath} is the path to the output file containing the pairwise LD between all loci within 500kb on the same contig. The resulting file often has the .stat.gz suffix appended to the path.
# ${outputPDF} is the path to the output PDF file containing the LD plot. The output PDF file often has the .pdf suffix appended to the path.
```

In summary, LD decays to near background level after 50kb. Background level differs depending on the dataset. LD decays slower in the NZ dataset (introduced populations) when compared to the IND dataset. A 100kb threshold was conservatively chosen to remove strong linkage disequilibrium in all datasets, retaining a good number of SNPs (about 5000-6000 per dataset). This was done using the following VCFtools command:
```bash
vcftools --gzvcf ${vcfa} --thin 100000 --recode --recode-INFO-all --stdout | bcftools filter -e 'AC==0 || AC=AN' -Oz -o ${vcfthin}
# ${vcfa} is the path to the vcf.gz file to thin. The command was run with vcf.gz files before and after the removal of singletons and doubletons.
# ${vcfthin} is the path to the resulting vcf.gz file after thinning. Note that BCFtools was used to ensure removal of monomorphic SNPs and compress the resulting VCF file.
```

## References
* Andrews S (2010) FASTQC. A quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Catchen J, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) Stacks: an analysis tool set for population genomics. Molecular ecology 22:3124–3140. https://doi.org/10.1111/MEC.12354
* Chen S, Zhou Y, Chen Y, Gu J (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34:i884–i890. https://doi.org/10.1093/BIOINFORMATICS/BTY560
* Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. (2011) The variant call format and VCFtools. Bioinformatics 27:2156–2158. https://doi.org/10.1093/bioinformatics/btr330
* Ewart KM, Griffin AS, Johnson RN, Kark S, Magory Cohen T, Lo N, et al. (2019) Two speed invasion: assisted and intrinsic dispersal of common mynas over 150 years of colonization. Journal of Biogeography 46:45–57. https://doi.org/10.1111/jbi.13473
* Ewels P, Magnusson M, Lundin S, Käller M (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32:3047–3048. https://doi.org/10.1093/bioinformatics/btw354
* Li H (2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics 27:2987–2993. https://doi.org/10.1093/bioinformatics/btr509
* Li H, Durbin R (2009) Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25:1754–1760. https://doi.org/10.1093/BIOINFORMATICS/BTP324
* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. (2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics 25:2078–2079. https://doi.org/10.1093/bioinformatics/btp352


