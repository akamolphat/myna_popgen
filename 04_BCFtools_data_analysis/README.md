# 04_BCFtools_data_analysis

This folder contain the scripts that perform the analyses on the SNP data called by BCFtools.

The script(s) found here generates results for the main text, but also Appendix S8 and S9.

The following scripts can be found inside the scripts folder as follows:

```bash
├── scripts
│   ├── NZ_pop_struc_BCFtools.R
│   ├── IND_pairwiseFST_from_VCF.R
│   ├── IND_pop_struc_BCFtools.R
│   ├── ALL_pairwiseFST_from_VCF.R
│   ├── ALL_pop_struc_BCFtools.R
│   ├── ALL_pop_struc_BCFtools_S8_5_to_S8_7.R
│   ├── NZ_popdiv.R
│   ├── IND_popdiv.R
│   ├── ALL_popdiv.R
│   └── ALL_SFS.R
├── data
│   └── processed
│       ├── HP-Rare
│       │   ├── glNZsub2_nmin6_ngene12_HP-rare_out.txt
│       │   ├── glINDsub2_nmin6_ngene10_HP-rare_out.txt
│       │   ├── glALLmerge2_ngene10_HP-rare_out.txt
│       │   └── glALLsub2_nmin6_ngene10_HP-rare_out.txt
│       └── ALL_BCFtools_subset_nmax20.csv
└── README.md
```
## Additional data
`data/processed/ALL_BCFtools_subset_nmax20.csv` contains the subset of individuals used when subsampling was done. This is to ensure the same individuals as in the paper.
`data/processed/HP-Rare/*_HP-rare_out.txt` contains output of from HP-rare for making Table 1 and Table S9.1-3.

## Population structure 

The population pairwiseFST matrix were to be calculated first as follows:

```bash
Rscript scripts/IND_pairwiseFST_from_VCF.R 
Rscript scripts/ALL_pairwiseFST_from_VCF.R 
```

The next three scripts run population structure analyses on NZ, IND, and ALL dataset, and can be executed in any order:

```bash
Rscript scripts/NZ_pop_struc_BCFtools.R
Rscript scripts/IND_pop_struc_BCFtools.R
Rscript scripts/ALL_pop_struc_BCFtools.R
```

* `NZ_pop_struc_BCFtools.R` generate results corresponding to Figure 2 in the main text, and Figures S8.1-S8.7 from Appendix S8.1.
* `IND_pop_struc_BCFtools.R` generate results corresponding to Figure 3 in the main text, and Figures S8.8-S8.13 from Appendix S8.2.
* `ALL_pop_struc_BCFtools.R` generate results corresponding to Figures 4 and 5 in the main text, and Figures S8.14-S8.22 from Appendix S8.3 and S8.4.

## Genetic diversity indices 
If run from scratch, these scripts are to be executed interactively as it first outputs data for HP-rare to calculate rarefied allelic richness and private allelic richness via a HP-rare GUI. However, the results from HP-rare have been uploaded here (`data/processed/HP-Rare/*_HP-rare_out.txt`) and the scripts can therefore be sourced/ran in one go (e.g. Rscript scripts/NZ_popdiv.R)

* `scripts/NZ_popdiv.R` calculates population diversity metrics and output Table S9.1.
* `scripts/IND_popdiv.R` calculates population diversity metrics and output Table S9.2.
* `scripts/ALL_popdiv.R` calculates population diversity metrics and output Table 1 and Table S9.3.

### Site Frequency Spectrum
Site frequency spectrums were constructed using `ALL_SFS.R`. This script require the input "ALL" dataset here were filtered differently from the other datasets described above. The data are available on request or it can be generated using the filters described in Appendix S3 and S4, or in the 02_variant_calling folder of this repository, except for the filters for singletons and doubletons. This is because, while filtering out singletons and doubletons avoid inclusion of incorrect genotyping error, it also causes biases in the SNPs retained per population. The singleton/doubleton filter will remove more rare SNPs from populations (or group of individuals with similar genetic background) with lower number of individuals.

* `ALL_SFS.R` generates results corresponding to Figure 6 in the main text, and Figure S9.1 from Appendix S9.

## Additional population structure analyses
Additional population structure analyses were performed on ALL dataset to creat Appendix S8.5 - S8.7 using the following script:

* `ALL_pop_struc_BCFtools_S8_5_to_S8_7.R`

However, the input "ALL" dataset here were filtered differently from the other datasets described above. The data are available on request or it can be generated using the filters described in Appendix S3, or in the section **6. Filter the variants called via BCFtools and STACKS** in the 02_variant_calling folder of this repository. 
