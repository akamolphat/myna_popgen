# IND genetic diversity
# This script calculates population diversity metrics
library(PopGenReport)
library(data.table)
library(dartR)
library(xlsx)
source("../shared_scripts/functions.R")

# 1. Define input file ----------------------------------------------#####
vcffileIND <- "../01_download_data/BCFtools_IND/variant_calls.IND.bialminQ30minGQ30DP15-100.norep.highnegfis.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2v2.csv"

# 2. Read in VCF file -----------------------------------------------#####
glIND <- gl.read.vcf(vcffileIND)
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile, na.strings = "n/a")
# Subset metadt for the genlight object
indiv_namesIND <- gsub('[a-z]$', '', indNames(glIND))
metadtsub <- metadt %>% 
  filter(ID %in% indiv_namesIND) %>%
  arrange(factor(ID, levels = indiv_namesIND)) 
# Attach individual metadata to genlight object ---------------------
glIND@other$ind.metrics <- metadtsub
# Assign individual ID to individual metadata
glIND@other$ind.metrics$ID <- indNames(glIND)
# Attach population information to genlight object
pop(glIND) <- glIND$other$ind.metrics$popdef1

### Calculate observed and expected heterozygousity, allelic richness and private allelic richness and allele frequency differences
glINDsub <- gl.drop.pop(glIND, pop.list = c("Gujarat"), recalc = T, mono.rm = T)
# Calculate obs and exp heterozygousity
dthet_INDsub <- gl.report.heterozygosity(glINDsub)
# Calculate allelic richness
# Needs to make sure that there are genotypes observed within each population
glINDsub2 <- gl.filter.callrate(glINDsub, method='pop', threshold=0.8, verbose=3, recalc = T, mono.rm = T)
## Convert to genind object
giINDsub2 <- gl2gi(glINDsub2)
## Calculate allelic richness
dtINDsub2 <- allel.rich(giINDsub2)

## Export to .gp for HP-rare ----------------------------------------#####
dir.create("data/processed/HP-Rare", showWarnings = F, recursive = T)
# All
gl2genalex(glINDsub2, outfile = 'glINDsub2_nmin6_CR08perpop.csv', outpath = './data/processed/HP-Rare/')
genalex2genepop("data/processed/HP-Rare/glINDsub2_nmin6_CR08perpop.csv", "data/processed/HP-Rare/glINDsub2_nmin6_CR08perpop.gp", "IND dataset (ROM and modern kept separate) after removal of populations with n <= 5. Filter for loci with CR >= 0.8 for each population as pop with missing loci cause error within allel.rich function in R.")
## Calculate rarefied AR and Private Allelic Richness in HP-rare ----#####
# Merge the tables NOTE that H0 = observed Heterozygousity
## All
dthet_INDsub$mean_allele_richness_allele.rich <- dtINDsub2$mean.richness
## Add data from HP-rare
### All
#### Read full dataset and calculate mean allelic richness ourselves
dtINDsub2_HP.rare_AR <- fread("data/processed/HP-Rare/glINDsub2_nmin6_ngene10_HP-rare_out.txt", skip = 9, nrows = 9)
dtINDsub2_HP.rare_allel_rich <- data.frame(pop = dtINDsub2_HP.rare_AR[,c(1)],
                                          mean_allelic_richness_HP.rare = round(rowMeans(dtINDsub2_HP.rare_AR[,-c(1,2,3)]),3))
colnames(dtINDsub2_HP.rare_allel_rich) <- c("pop", "mean_allelic_richness_HP.rare")
dtINDsub2_HP.rare_allel_rich$n_rarefied_nindx2 <- 10
dtINDsub2_HP.rare_PAR <- fread("data/processed/HP-Rare/glINDsub2_nmin6_ngene10_HP-rare_out.txt", skip = 23, nrows = 9)
dtINDsub2_HP.rare_p_allel_rich <- data.frame(pop = dtINDsub2_HP.rare_PAR[,c(1)],
                                            mean_private_allelic_richness_HP.rare = round(rowMeans(dtINDsub2_HP.rare_PAR[,-c(1,2,3)]),3))

colnames(dtINDsub2_HP.rare_p_allel_rich) <- c("pop", "mean_private_allelic_richness_HP.rare")
dtdiv_glINDsub2 <- merge(dthet_INDsub, dtINDsub2_HP.rare_allel_rich, all.x = T)
dtdiv_glINDsub2 <- merge(dtdiv_glINDsub2, dtINDsub2_HP.rare_p_allel_rich, all.x = T)

# Calculate proportion of polymorphic loci
dtpoly_glINDsub2 <- calc_prop_polymorphic_loci_perpop_rarefact(glINDsub2,
                                                              n_per_pop = 6, 
                                                              n_iter = 100)

dtdiv_glINDsub2 <- merge(dtdiv_glINDsub2, dtpoly_glINDsub2$dt_summary)

# Output results to xlsx --------------------------------------------#####
dir.create("results/Pop_div/")
filename <- "results/Pop_div/Pop_div_IND.xlsx"
write.xlsx(dtdiv_glINDsub2, filename, sheetName = "TableS9_2_full", row.names = F)
write.xlsx(dtdiv_glINDsub2[,c("pop", "Ho", "He", "mean_allele_richness_allele.rich", "mean_allelic_richness_HP.rare", "mean_private_allelic_richness_HP.rare", "median_prop_poly_loci")], filename, 
           sheetName = "TableS9_2", row.names = F, append = T)
