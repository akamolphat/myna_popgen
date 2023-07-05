# This script calculates genetic diversity metrics for the 
# ALL BCFtools dataset
#
# ALL genetic diversity
# This script calculates population diversity metrics
# Libraries ---------------------------------------------------------
library(PopGenReport)
library(data.table)
library(dartR)
library(xlsx)
source("../shared_scripts/functions.R")

# Define input file -------------------------------------------------
vcffile <- "../01_download_data/BCFtools_ALL/variant_calls.ALL.bialminQ30minGQ30DP15-125.norep.noadm.highnegfis.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2.csv"
# Read in VCF file --------------------------------------------------
glALL <- gl.read.vcf(vcffileALL)
# Subset metadt for the genlight object
indiv_namesALL <- gsub('[a-z]$', '', indNames(glALL))
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile, na.strings = "n/a")
# Subset metadt for the genlight object
indiv_names <- gsub('[a-z]$', '', indNames(glALL))
metadtsub <- metadt %>% 
  filter(ID %in% indiv_names) %>%
  arrange(factor(ID, levels = indiv_names)) 
# Attach individual metadata to genlight object ---------------------
glALL@other$ind.metrics <- metadtsub
# Assign individual ID to individual metadata
glALL@other$ind.metrics$ID <- indNames(glALL)
# Recode popdef2 with Maharashtra subpop. A separated ---------------
glALL$other$ind.metrics$ALL_pop_recode <- glALL$other$ind.metrics$popdef2
glALL$other$ind.metrics <- glALL$other$ind.metrics %>%
  mutate(ALL_pop_recode = replace(ALL_pop_recode, id %in% paste("M0", seq(210, 215), sep = ""), "Maharashtra subpop. A"))
# Export to HP-rare -------------------------------------------------
## popdef1 ----------------------------------------------------------
### Export to .gp for HP-rare ---------------------------------------
pop(glALL) <- glALL$other$ind.metrics$popdef1
glALLmerge2 <- gl.filter.callrate(glALL, method='pop', threshold=0.8, verbose=3, recalc = T, mono.rm = T)
dir.create("data/processed/HP-Rare", showWarnings = F)
gl2genalex(glALLmerge2, outfile = 'glALLmerge2_CR08perpop.csv', outpath = './data/processed/HP-Rare/', overwrite = T)
genalex2genepop("data/processed/HP-Rare/glALLmerge2_CR08perpop.csv", "data/processed/HP-Rare/glALLmerge2_CR08perpop.gp", "ALL dataset (populations within NZ: Other and India: Other merged). Filter for loci with CR >= 0.8 for each population as pop with missing loci cause error within allel.rich function in R.")

## popdef2 with Maharashtra subpop. A separated ---------------------
### Export to .gp for HP-rare ---------------------------------------
gl2genalex(glALLsub2, outfile = 'glALLsub2_nmin6_CR08perpop.csv', outpath = './data/processed/HP-Rare/')
genalex2genepop("data/processed/HP-Rare/glALLsub2_nmin6_CR08perpop.csv", "data/processed/HP-Rare/glALLsub2_nmin6_CR08perpop.gp", "ALL dataset (ROM and modern kept separate) after removal of populations with n <= 5. Filter for loci with CR >= 0.8 for each population as pop with missing loci cause error within allel.rich function in R.")
## Calculate rarefied AR and PAR in HP-rare -------------------------
## AR = Allelic richness, PAR = Private allelic richness


# Calculate pop. diversity in R -------------------------------------
## popdef1 ----------------------------------------------------------
pop(glALL) <- glALL$other$ind.metrics$popdef1
### Ho and He -------------------------------------------------------
dthet_ALLmerge2 <- gl.report.heterozygosity(glALL)
### PopGenReport AR -------------------------------------------------
# Needs to make sure that there are genotypes observed within each population
glALLmerge2 <- gl.filter.callrate(glALL, method='pop', threshold=0.8, verbose=3, recalc = T, mono.rm = T)
## Convert to genind object
giALLmerge2 <- gl2gi(glALLmerge2)
## Calculate allelic richness
dtALLmerge2 <- allel.rich(giALLmerge2)
### Proportion of polymorphic loci-----------------------------------
dtpoly_glALLmerge2 <- calc_prop_polymorphic_loci_perpop_rarefact(glALLmerge2, n_per_pop = 6, n_iter = 100)
dthet_ALLmerge2$mean_allele_richness_allele.rich <- dtALLmerge2$mean.richness
### Merge the outputs in R with outputs from HP-rare -----------------
#### Read AR from HP-rare
dtALLmerge_HP.rare_AR2 <- fread("data/processed/HP-Rare/glALLmerge2_ngene10_HP-rare_out.txt", skip = 9, nrows = 11)
dtALLmerge_HP.rare_allel_rich2 <- data.frame(pop = dtALLmerge_HP.rare_AR2[,c(1)],
                                             mean_allelic_richness_HP.rare = round(rowMeans(dtALLmerge_HP.rare_AR2[,-c(1,2,3)]),3))
colnames(dtALLmerge_HP.rare_allel_rich2) <- c("pop", "mean_allelic_richness_HP.rare")
dtALLmerge_HP.rare_allel_rich2$n_rarefied_nindx2 <- 10
#### Read PAR from HP-rare
dtALLmerge_HP.rare_PAR2 <- fread("data/processed/HP-Rare/glALLmerge2_ngene10_HP-rare_out.txt", skip = 24, nrows = 11)
dtALLmerge_HP.rare_p_allel_rich2 <- data.frame(pop = dtALLmerge_HP.rare_PAR2[,c(1)],
                                               mean_private_allelic_richness_HP.rare = round(rowMeans(dtALLmerge_HP.rare_PAR2[,-c(1,2,3)]),3))
colnames(dtALLmerge_HP.rare_p_allel_rich2) <- c("pop", "mean_private_allelic_richness_HP.rare")
#### Merge tables
dtdiv_glALLmerge2 <- merge(dthet_ALLmerge2, dtALLmerge_HP.rare_allel_rich2, all.x = T)
dtdiv_glALLmerge2 <- merge(dtdiv_glALLmerge2, dtALLmerge_HP.rare_p_allel_rich2, all.x = T)
dtdiv_glALLmerge2 <- merge(dtdiv_glALLmerge2, dtpoly_glALLmerge2$dt_summary)
## popdef2 with Maharashtra subpop. A separated ---------------------
pop(glALL) <- glALL$other$ind.metrics$ALL_pop_recode
### Remove populations with n < 5 -----------------------------------
glALLsub <- gl.drop.pop(glALL, pop.list = c("Gujarat", "Taupo (ROM)", "Auckland", "Kaikohe (ROM)"), recalc = T, mono.rm = T)
### Ho and He -------------------------------------------------------
dthet_ALLsub <- gl.report.heterozygosity(glALLsub)
# Calculate allelic richness
# Needs to make sure that there are genotypes observed within each population
glALLsub2 <- gl.filter.callrate(glALLsub, method='pop', threshold=0.8, verbose=3, recalc = T, mono.rm = T)
## Convert to genind object
giALLsub2 <- gl2gi(glALLsub2)
## Calculate allelic richness
dtALLsub2 <- allel.rich(giALLsub2)
### Proportion of polymorphic loci-----------------------------------
dtpoly_glALLsub2 <- calc_prop_polymorphic_loci_perpop_rarefact(glALLsub2,
                                                               n_per_pop = 6, 
                                                               n_iter = 100)
## All
dthet_ALLsub$mean_allele_richness_allele.rich <- dtALLsub2$mean.richness

### Merge the outputs in R with outputs from HP-rare -----------------
#### Read AR from HP-rare
dtALLsub2_HP.rare_AR <- fread("data/processed/HP-Rare/glALLsub2_nmin6_ngene10_HP-rare_out.txt", skip = 9, nrows = 28)
dtALLsub2_HP.rare_allel_rich <- data.frame(pop = dtALLsub2_HP.rare_AR[,c(1)],
           mean_allelic_richness_HP.rare = round(rowMeans(dtALLsub2_HP.rare_AR[,-c(1,2,3)]),3))
colnames(dtALLsub2_HP.rare_allel_rich) <- c("pop", "mean_allelic_richness_HP.rare")
dtALLsub2_HP.rare_allel_rich$n_rarefied_nindx2 <- 10
#### Read PAR from HP-rare
dtALLsub2_HP.rare_PAR <- fread("data/processed/HP-Rare/glALLsub2_nmin6_ngene10_HP-rare_out.txt", skip = 42, nrows = 28)
dtALLsub2_HP.rare_p_allel_rich <- data.frame(pop = dtALLsub2_HP.rare_PAR[,c(1)],
                                             mean_private_allelic_richness_HP.rare = round(rowMeans(dtALLsub2_HP.rare_PAR[,-c(1,2,3)]),3))
colnames(dtALLsub2_HP.rare_p_allel_rich) <- c("pop", "mean_private_allelic_richness_HP.rare")
#### Merge tables
dtdiv_glALLsub2 <- merge(dthet_ALLsub, dtALLsub2_HP.rare_allel_rich, all.x = T)
dtdiv_glALLsub2 <- merge(dtdiv_glALLsub2, dtALLsub2_HP.rare_p_allel_rich, all.x = T)
dtdiv_glALLsub2 <- merge(dtdiv_glALLsub2, dtpoly_glALLsub2$dt_summary)

# Output results to xlsx --------------------------------------------
dir.create("results/Pop_div/")
filename <- "results/Pop_div/Pop_div_ALL.xlsx"
write.xlsx(dtdiv_glALLmerge2[,c("pop", "Ho", "He", "mean_allele_richness_allele.rich", "mean_allelic_richness_HP.rare", "mean_private_allelic_richness_HP.rare", "median_prop_poly_loci")], filename, 
           sheetName = "Table1", row.names = F)
write.xlsx(dtdiv_glALLsub2[,c("pop", "Ho", "He", "mean_allele_richness_allele.rich", "mean_allelic_richness_HP.rare", "mean_private_allelic_richness_HP.rare", "median_prop_poly_loci")], filename, 
           sheetName = "Table9_3", row.names = F, append = T)
write.xlsx(dtdiv_glALLsub2, filename, sheetName = "Table1_full", row.names = F, append = T)
write.xlsx(dtdiv_glALLmerge2, filename, sheetName = "Table9_3_full", row.names = F, append = T)

