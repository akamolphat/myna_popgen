# ALL_pairwiseFST_from_VCF.R
# This script calculates pairwiseFST in the Indian dataset
# Define input variables --------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# vcffile <- args[1]
# metadtfile <- args[2]
# outputxlsx <- args[3]
vcffile <- "../01_download_data/BCFtools_ALL/variant_calls.ALL.bialminQ30minGQ30DP15-125.norep.noadm.highnegfis.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2.csv"
outputxlsx <- "results/FST/03_ALL/pairwiseFST_ALL_BCF_hwe_100bs.xlsx"

# Load libraries ----------------------------------------------------
source("../shared_scripts/functions.R")
library(dartR)
library(tidyverse)
library(LEA)
library(pcadapt)
library(vcfR)
library(xlsx)
# Read in VCF file --------------------------------------------------
gl <- gl.read.vcf(vcffile)
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile, na.strings = "n/a")
# Subset metadt for the genlight object
indiv_names <- gsub('[a-z]$', '', indNames(gl))
metadtsub <- metadt %>% 
  filter(ID %in% indiv_names) %>%
  arrange(factor(ID, levels = indiv_names)) 
# Attach individual metadata to genlight object ---------------------
gl@other$ind.metrics <- metadtsub
# Assign individual ID to individual metadata
gl@other$ind.metrics$ID <- indNames(gl)
# Attach population information to genlight object
# Number of bootstraps -------------------------------------------
bootn <- 100
bootcol <- 6+bootn-1 
# FST with popdef1 -----------------------------------------------
pop(gl) <- gl$other$ind.metrics$popdef1
# FST between all popdef1 ----------------------------------------
fstmat <- gl.fst.pop(gl, nboots = bootn)
# Reich FST between all popdef1 ----------------------------------
reichfstmat <- reich.fst(gl, bootstrap=bootn)
reichfstmat$bootstraps$p_value <- rowSums(reichfstmat$bootstraps[,6:bootcol]<0)/bootn 
reichpvalmat <- pivot_wider(reichfstmat$bootstraps[,c("pop1", "pop2", "p_value")], names_from = pop2, values_from = p_value)

# Add spreadsheet with popdef1
write.xlsx(fstmat$Fsts, file = outputxlsx, sheetName = "FST_100bs", append = F)
write.xlsx(fstmat$Pvalues, file = outputxlsx, sheetName = "pval_100bs", append = T)
write.xlsx(fstmat$Bootstraps, file = outputxlsx, sheetName = "FST_100bs_val", append = T)
write.xlsx(reichfstmat$fsts, file = outputxlsx, sheetName = "ReichFST_100bs", append = T)
write.xlsx(reichpvalmat, file = outputxlsx, sheetName = "Reichpval_100bs", append = T)
write.xlsx(reichfstmat$bootstraps, file = outputxlsx, sheetName = "ReichFST_100bs_val", append = T)

# Perform FST with popdef2 ------------------------------------------
gl$other$ind.metrics$ALL_pop_recode <- gl$other$ind.metrics$popdef2
gl$other$ind.metrics <- gl$other$ind.metrics %>%
  mutate(ALL_pop_recode = replace(ALL_pop_recode, ALL_pop_recode == "Central Division", "Fiji")) %>%
  mutate(ALL_pop_recode = replace(ALL_pop_recode, ALL_pop_recode == "Gauteng", "South Africa")) %>%
  mutate(ALL_pop_recode = replace(ALL_pop_recode, id %in% paste("M0", seq(210, 215), sep = ""), "Maharashtra subpop. A")) 
pop(gl) <- gl$other$ind.metrics$ALL_pop_recode
# Standard FST
fstmat <- gl.fst.pop(gl, nboots = bootn)
# Reich FST
reichfstmat <- reich.fst(gl, bootstrap=bootn)
reichfstmat$bootstraps$p_value <- rowSums(reichfstmat$bootstraps[,6:bootcol]<0)/bootn 
reichpvalmat <- pivot_wider(reichfstmat$bootstraps[,c("pop1", "pop2", "p_value")], names_from = pop2, values_from = p_value)

# Add spreadsheet with popdef2
write.xlsx(fstmat$Fsts, file = outputxlsx, sheetName = "FST_100bs_modROMsep", append = T)
write.xlsx(fstmat$Pvalues, file = outputxlsx, sheetName = "pval_100bs_modROMsep", append = T)
write.xlsx(fstmat$Bootstraps, file = outputxlsx, sheetName = "FST_100bs_val_modROMsep", append = T)
write.xlsx(reichfstmat$fsts, file = outputxlsx, sheetName = "ReichFST_100bs_modROMsep", append = T)
write.xlsx(reichpvalmat, file = outputxlsx, sheetName = "Reichpval_100bs_modROMsep", append = T)
write.xlsx(reichfstmat$bootstraps, file = outputxlsx, sheetName = "ReichFST_100bs_val_modROMsep", append = T)

