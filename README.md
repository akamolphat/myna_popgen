# myna_popgen

Scripts and codes related to the following manuscript and data:

* Atsawawaranunt, K., Ewart, K.M., Major, R.E. et al. Tracing the introduction of the invasive common myna using population genomics. Heredity 131, 56–67 (2023). https://doi.org/10.1038/s41437-023-00621-w
* https://doi.org/10.5061/dryad.xsj3tx9m7

This repository is divided into different folders for different tasks. Each folder has a README.md file which describes the tasks being executed in each folder. This repository is structured as follows:

```bash
├── 01_download_data
│   └── README.md
│   
├── 02_variant_calling
│   └── README.md
|
├── 03_plot_map
|   ├── scripts
│   └── README.md
|
├── 04_BCFtools_data_analysis
|   ├── data/processed
|   ├── scripts
│   └── README.md
|
├── 05_STACKS_data_analysis
|   ├── scripts
│   └── README.md
|
├── 06_DART_data_analysis
|   ├── scripts
│   └── README.md
|
├── README.md
└── .gitignore
```

* `01_download_data` is where the data is downloaded. Sequence and sample data are to be downloaded from [Dryad](https://doi.org/10.5061/dryad.xsj3tx9m7), shapefile for the common myna distribution is to be downloaded from the IUCN redlist website.
* `02_variant_calling` documents how the variants were called and filtered.
* `03_plot_map` plots maps of the species and sample distribution (Figure 1).
* `04_BCFtools_data_analysis` performs population structure analyses and calculate population genetic diversity metrics based on the SNPs called via BCFtools. This section consititutes the main analyses and makes Figures 2-6, Figures in Appendix S8, Table 1, and Tables in Appendix S9.
* `05_STACKS_data_analysis` performs population structure analyses and calculate population genetic diversity metrics based on the SNPs called via STACKS. This section makes figures in Appendix S10.1-S10.3.
* `06_DART_data_analysis` performs population structure analyses and calculate population genetic diversity metrics based on the SNPs called via DArT P/L proprietary DArTsoft14 SNP calling pipeline. This section makes figures in Appendix S10.4-S10.6.



