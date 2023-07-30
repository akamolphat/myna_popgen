# NZ genetic diversity
# This script calculates population diversity metrics
library(PopGenReport)
library(data.table)
library(dartR)
library(xlsx)
source("../shared_scripts/functions.R")

# 1. Define input file ----------------------------------------------#####
vcffileNZ <- "../01_download_data/BCFtools_NZ/variant_calls.NZ.bialminQ30minGQ30DP15-100.norep.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2v2.csv"

# 2. Read in VCF file -----------------------------------------------#####
glNZ <- gl.read.vcf(vcffileNZ)
# Subset metadt for the genlight object
indiv_namesNZ <- gsub('[a-z]$', '', indNames(glNZ))
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile, na.strings = "n/a")
metadtsubNZ <- metadt %>% 
  filter(ID %in% indiv_namesNZ) %>%
  arrange(factor(ID, levels = indiv_namesNZ)) 
# Attach individual metadata to genlight object ---------------------#####
glNZ@other$ind.metrics <- metadtsubNZ
# Assign individual ID to individual metadata
glNZ@other$ind.metrics$ID <- indNames(glNZ)
# 3. Calculate pop diversity with populations separated --------------#####
# Attach population information to genlight object
pop(glNZ) <- glNZ$other$ind.metrics$popdef1
### Calculate observed and expected heterozygousity, allelic richness and private allelic richness and allele frequency differences
glNZsub <- gl.drop.pop(glNZ, pop.list = c("Taupo (ROM)", "Auckland", "Kaikohe (ROM)"), recalc = T, mono.rm = T)
# Calculate obs and exp heterozygousity
dthet_NZsub <- gl.report.heterozygosity(glNZsub)
# Calculate allelic richness
# Needs to make sure that there are genotypes observed within each population
glNZsub2 <- gl.filter.callrate(glNZsub, method='pop', threshold=0.8, verbose=3, recalc = T, mono.rm = T)
## Convert to genind object
giNZsub2 <- gl2gi(glNZsub2)
## Calculate allelic richness
dtNZsub2 <- allel.rich(giNZsub2)

## Export to .gp for HP-rare ----------------------------------------#####
dir.create("data/processed/HP-Rare", showWarnings = F, recursive = T)
# All
gl2genalex(glNZsub2, outfile = 'glNZsub2_nmin6_CR08perpop.csv', outpath = './data/processed/HP-Rare/')
genalex2genepop("data/processed/HP-Rare/glNZsub2_nmin6_CR08perpop.csv", "data/processed/HP-Rare/glNZsub2_nmin6_CR08perpop.gp", "NZ dataset (ROM and modern kept separate) after removal of populations with n <= 5. Filter for loci with CR >= 0.8 for each population as pop with missing loci cause error within allel.rich function in R.")
## Calculate rarefied AR and Private Allelic Richness in HP-rare ----#####
# Merge the tables NOTE that H0 = observed Heterozygousity
## All
dthet_NZsub$mean_allele_richness_allele.rich <- dtNZsub2$mean.richness
## Add data from HP-rare
### All
#### Read full dataset and calculate mean allelic richness ourselves
dtNZsub2_HP.rare_AR <- fread("data/processed/HP-Rare/glNZsub2_nmin6_ngene12_HP-rare_out.txt", skip = 9, nrows = 12)
dtNZsub2_HP.rare_allel_rich <- data.frame(pop = dtNZsub2_HP.rare_AR[,c(1)],
                                           mean_allelic_richness_HP.rare = round(rowMeans(dtNZsub2_HP.rare_AR[,-c(1,2,3)]),3))
colnames(dtNZsub2_HP.rare_allel_rich) <- c("pop", "mean_allelic_richness_HP.rare")
dtNZsub2_HP.rare_allel_rich$n_rarefied_nindx2 <- 12
dtNZsub2_HP.rare_PAR <- fread("data/processed/HP-Rare/glNZsub2_nmin6_ngene12_HP-rare_out.txt", skip = 25, nrows = 12)
dtNZsub2_HP.rare_p_allel_rich <- data.frame(pop = dtNZsub2_HP.rare_PAR[,c(1)],
                                             mean_private_allelic_richness_HP.rare = round(rowMeans(dtNZsub2_HP.rare_PAR[,-c(1,2,3)]),3))

colnames(dtNZsub2_HP.rare_p_allel_rich) <- c("pop", "mean_private_allelic_richness_HP.rare")
dtdiv_glNZsub2 <- merge(dthet_NZsub, dtNZsub2_HP.rare_allel_rich, all.x = T)
dtdiv_glNZsub2 <- merge(dtdiv_glNZsub2, dtNZsub2_HP.rare_p_allel_rich, all.x = T)

# Calculate proportion of polymorphic loci
dtpoly_glNZsub2 <- calc_prop_polymorphic_loci_perpop_rarefact(glNZsub2,
                                                               n_per_pop = 7, 
                                                               n_iter = 100)

dtdiv_glNZsub2 <- merge(dtdiv_glNZsub2, dtpoly_glNZsub2$dt_summary)

# Output results to xlsx --------------------------------------------#####
dir.create("results/Pop_div/")
filename <- "results/Pop_div/Pop_div_NZ.xlsx"
write.xlsx(dtdiv_glNZsub2, filename, 
           sheetName = "TableS9_1_full", row.names = F)
write.xlsx(dtdiv_glNZsub2[,c("pop", "Ho", "He", "mean_allele_richness_allele.rich", "mean_allelic_richness_HP.rare", "mean_private_allelic_richness_HP.rare", "median_prop_poly_loci")], filename, 
           sheetName = "TableS9_1", row.names = F, append = T)
