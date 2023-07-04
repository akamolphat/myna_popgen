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
│   └── ALL_pop_struc_BCFtools.R
├── scripts
│   └── ALL_BCFtools_subset_nmax20.csv
└── README.md
```

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

`NZ_pop_struc_BCFtools.R` generate results corresponding to Figure 2 in the main text, and Figures S8.1-S8.7 from Appendix S8.1.
`IND_pop_struc_BCFtools.R` generate results corresponding to Figure 3 in the main text, and Figures S8.8-S8.13 from Appendix S8.2.
`ALL_pop_struc_BCFtools.R` generate results corresponding to Figures 4 and 5 in the main text, and Figures S8.14-S8.22 from Appendix S8.3 and S8.4.