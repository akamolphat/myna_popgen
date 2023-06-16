# 06_DArT_data_analysis

This folder contain the scripts that perform the analyses on the SNP data called by DArT's proprietary DArTsoft14 SNP calling pipeline.

The script(s) found here generates results which can be found in Appendix S10.4-6

The following scripts can be found inside the scripts folder as follows:
```bash
├── scripts
│   ├── Create_metadata_file.R 
│   ├── NZ_popstruc_SI_DART.R
│   ├── IND_popstruc_SI_DART.R
│   └── ALL_popstruc_SI_DART.R
└── README.md
```

Create_metadata_file.R has to be executed first. This is to be executed only after the data has been downloaded from Dryad and the DART.zip unzipped in the 01_download_data folder. The next three scripts run population structure analyses on NZ, IND, and ALL dataset, respectively, and can be executed in any order.

```bash
Rscript scripts/Create_metadata_file.R
Rscript NZ_popstruc_SI_DART.R
```

