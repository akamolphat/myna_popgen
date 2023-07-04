# IND_pairwiseFST_from_VCF.R
# This script calculates pairwiseFST in the Indian dataset
# Define input variables --------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# vcffile <- args[1]
# metadtfile <- args[2]
# outputxlsx <- args[3]
vcffile <- "../01_download_data/BCFtools_IND/variant_calls.IND.bialminQ30minGQ30DP15-100.norep.highnegfis.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2.csv"
outputxlsx <- "results/FST/02_IND/pairwiseFST_IND_BCF_hwe_100bs.xlsx"

# Create output folder to store outputs from HWE filters
dir.create(path = "results/FST/02_IND", recursive = T, showWarnings = F)

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
# Define Pop4
gl$other$ind.metrics$IND_pop_recode <- gl$other$ind.metrics$popdef2
gl$other$ind.metrics <- gl$other$ind.metrics %>%
  mutate(IND_pop_recode = replace(IND_pop_recode, ID %in% paste("M0", seq(210, 215), sep = ""), "Maharashtra subpop. A"))
# pop(gl) <- gl$other$ind.metrics$IND_pop_recode
pop(gl) <- gl@other$ind.metrics$IND_pop_recode

# FST between all pop_INDrecode ----------------------------------
print("calc pairwise FST")
bootn <- 100
fstmat <- gl.fst.pop(gl, nboots = bootn)
# Reich FST between all pop_INDrecode ----------------------------
reichfstmat <- reich.fst(gl, bootstrap=bootn)
bootcol <- 6+bootn-1 
reichfstmat$bootstraps$p_value <- rowSums(reichfstmat$bootstraps[,6:bootcol]<0)/bootn 
reichpvalmat <- pivot_wider(reichfstmat$bootstraps[,c("pop1", "pop2", "p_value")], names_from = pop2, values_from = p_value)

# Attach population information to genlight object
pop(gl) <- gl$other$ind.metrics$popdef2
# FST between all popdef2 ----------------------------------
fstmat2 <- gl.fst.pop(gl, nboots = bootn)
# Reich FST between all popdef2 ----------------------------
reichfstmat2 <- reich.fst(gl, bootstrap=bootn)
bootcol <- 6+bootn-1 
reichfstmat2$bootstraps$p_value <- rowSums(reichfstmat2$bootstraps[,6:bootcol]<0)/bootn 
reichpvalmat2 <- pivot_wider(reichfstmat2$bootstraps[,c("pop1", "pop2", "p_value")], names_from = pop2, values_from = p_value)


# Output to xlsx file -----------------------------------------------
print("output to xlsx")
write.xlsx(fstmat$Fsts, file = outputxlsx, sheetName = "FST_100bs")
write.xlsx(fstmat$Pvalues, file = outputxlsx, sheetName = "pval_100bs", append = T)
write.xlsx(fstmat$Bootstraps, file = outputxlsx, sheetName = "FST_100bs_val", append = T)
write.xlsx(reichfstmat$fsts, file = outputxlsx, sheetName = "ReichFST_100bs", append = T)
write.xlsx(reichpvalmat, file = outputxlsx, sheetName = "Reichpval_100bs", append = T)
write.xlsx(reichfstmat$bootstraps, file = outputxlsx, sheetName = "ReichFST_100bs_val", append = T)
# Add spreadsheet with popdef2 matrix instead of pop_INDrecode
write.xlsx(fstmat2$Fsts, file = outputxlsx, sheetName = "FST_100bs_popdef2", append = T)
write.xlsx(fstmat2$Pvalues, file = outputxlsx, sheetName = "pval_100bs_popdef2", append = T)
write.xlsx(fstmat2$Bootstraps, file = outputxlsx, sheetName = "FST_100bs_val_popdef2", append = T)
write.xlsx(reichfstmat2$fsts, file = outputxlsx, sheetName = "ReichFST_100bs_popdef2", append = T)
write.xlsx(reichpvalmat2, file = outputxlsx, sheetName = "Reichpval_100bs_popdef2", append = T)
write.xlsx(reichfstmat2$bootstraps, file = outputxlsx, sheetName = "ReichFST_100bs_val_popdef2", append = T)