# This script makes plots for the figures in Appendix S8.5-S8.7
#
# This script is for the BCFtools dataset, on differently subset 
# datasets
#
# Load libraries ----------------------------------------------------
library(dartR)
library(LEA)
library(tidyverse)
library(cowplot)
library(xlsx)
library(pophelper)
source("../shared_scripts/functions.R")
# Define input files -----------------------------------------------
FULLVCF <- "../01_download_data/BCFtools_ALL/variant_calls.ALL.bialminQ30minGQ30DP15-125.norep.noadm.highnegfis.lmiss20.nosingledoubletons.snps.vcf.gz"
# FULLVCF is the path to the vcf file that has not been thinned. This path will have to be changed to match the file path 
metadtfile <- "../01_download_data/TableS1.2v2.csv"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/PCA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)

# Read in VCF -------------------------------------------------------
glFULL <- gl.read.vcf(FULLVCF)
# Get fixed section from VCF for contig and position 
vcf_in <- read.vcfR(FULLVCF, verbose = FALSE )
vcfFIX <- getFIX(vcf_in)
vcfFIXmerged <- cbind(vcfFIX, glFULL$other$loc.metrics)
glFULL$other$loc.metrics <- vcfFIXmerged
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile, na.strings = "n/a")
# Subset metadt for the genlight object
indiv_names <- gsub('[a-z]$', '', indNames(glFULL))
metadtsub <- metadt %>% 
  filter(id %in% indiv_names) %>%
  arrange(factor(id, levels = indiv_names)) 
# Attach individual metadata to genlight object ---------------------
gl@other$ind.metrics <- metadtsub
# Assign individual ID to individual metadata
gl@other$ind.metrics$ID <- indNames(gl)
# Attach population information to genlight object
pop(gl) <- gl$other$ind.metrics$popdef2
# SI S8.5: No SNP thinning ------------------------------------------
## PCA on no vcfthin ------------------------------------------------
pcFULL <- gl.pcoa(glFULL, nfactors=10, parallel = T, n.cores = 8)
# save(pcFULL, file = "results/PCA/03_ALL/PCA_ALL_FULL_NOTHIN.Rdata")
# load(file = "results/PCA/03_ALL/PCA_ALL_FULL_NOTHIN.Rdata")
### Figure S8.23: scree plot ----------------------------------------
png("results/Figs_paper/SI/ALL_FULL_NOTHIN_PCA_BCF_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
PCA_screeplot(pcFULL, nPC = 15)
dev.off()

# Quick plot
gl.pcoa.plot(pcFULL, glFULL, xaxis = 1, yaxis = 2) 
# Create a table with variance explained
dt_var_expl_FULL <- PCA_var_explained(pcFULL)
round(dt_var_expl_FULL$variance_explained, 1)
# Create table with PC scores and populations
dtpcFULL <- data.frame(pcFULL$scores)
dtpcFULL$pop <- factor(glFULL$other$ind.metrics$popdef1,
                       levels = dtcolshapes$pop)
### Figure S8.24: PCA plot ------------------------------------------
pFULL_1v2 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 1, yaxis = 2) + PCA_theme

pFULL_2v3 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 2, yaxis = 3) + PCA_theme

pFULL_3v4 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 3, yaxis = 4) + PCA_theme

pFULL_4v5 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 4, yaxis = 5) + PCA_theme

pFULL_5v6 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 5, yaxis = 6) + PCA_theme

pFULL_6v7 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 6, yaxis = 7) + PCA_theme

pFULL_7v8 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 7, yaxis = 8) + PCA_theme

pFULL_8v9 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                       dtcolshapes, xaxis = 8, yaxis = 9) + PCA_theme

pFULL_9v10 <- plotPCASI(dtpcFULL, varls = round(dt_var_expl_FULL$variance_explained, 1),
                        dtcolshapes, xaxis = 9, yaxis = 10) + PCA_theme

legend_b2 <- get_legend(
  pFULL_1v2  +
    guides(colour=guide_legend(nrow=2,byrow=T), 
           shape=guide_legend(nrow=2,byrow=T)) 
)

pca_SIFULL <- plot_grid(pFULL_1v2 + theme(legend.position = "none"), 
                        pFULL_2v3 + theme(legend.position = "none"), 
                        pFULL_3v4 + theme(legend.position = "none"),
                        pFULL_4v5 + theme(legend.position = "none"), 
                        pFULL_5v6 + theme(legend.position = "none"), 
                        pFULL_6v7 + theme(legend.position = "none"), 
                        pFULL_7v8 + theme(legend.position = "none"), 
                        pFULL_8v9 + theme(legend.position = "none"), 
                        pFULL_9v10 + theme(legend.position = "none"), 
                        labels = "AUTO", ncol = 3, nrow = 3)

png("results/Figs_paper/SI/ALL_PCA_BCF_FULL_NOTHIN.png", width = 8.3, height = 10, units = "in", res = 600)
plot_grid(pca_SIFULL, legend_b2, ncol = 1, rel_heights = c(2, .15))
dev.off()

### Figure S8.25: PC loadings ---------------------------------------
dtpcloadFULL <- data.frame(pcFULL$loadings)
dtpcloadFULL$contig <- as.numeric(glFULL$other$loc.metrics$CHROM)
dtpcloadFULL$position <- as.numeric(glFULL$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcloadFULL <- dtpcloadFULL %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/ALL_FULL_NOTHIN_PCA_BCF_loadings.png", width = 8.3, height = 11.7, units = "in", res = 600)

op <- par(mfrow = c(10,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

for (i in 1:10){
  plot(dtpcloadFULL[,i], ylab = paste("PC", i, sep = ""), xaxt = "n", xlab = "", cex = 0.25)
  axis(side = 2, labels = T)
  box(which = "plot", bty = "l")
}

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5)
dev.off()


## sNMF on no vcfthin -----------------------------------------------
genofileFULL <- "data/processed/LEA/03_ALL/ALL_FULL_NOTHIN.geno"
# Convert genlight to matrix of 0,1,2 
mat_glFULL <- t(as.matrix(glFULL))
# Write table to .geno file which is essentially a matrix 
write.table(mat_glFULL, file = genofileFULL,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)

Krange <- 1:20
# snmf_ALLFULL <- snmf(genofileFULL, K = Krange, repetitions = 10, ploidy = 2,
#                     entropy = T, alpha = 100, project = "new")
snmf_ALLFULL <- load.snmfProject("data/processed/LEA/03_ALL/ALL_FULL_NOTHIN.snmfProject")
plot(snmf_ALLFULL, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), ALL dataset, after HWE filter, only relevant pops.")
# Calculate mean entropy
crossentropy.matALLFULL <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_ALLFULL, K = x)})))
mean.entropyALLFULL   <- apply(crossentropy.matALLFULL, MARGIN=1, FUN=mean, na.rm=TRUE)
if(any(diff(mean.entropyALLFULL)>0)){
  bestK <- unname(which(diff(mean.entropyALLFULL)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
bestK

### Figure S8.26: Mean cross-entropy plot----------------------------
png("results/Figs_paper/SI/ALL_FULL_NOTHIN_sNMF_meanxentropy_k1-20.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropyALLFULL, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()

q.ALLFULL.ls <- merge_sNMF_qmatrix_multiK(snmf_ALLFULL, glFULL, k2plot = 1:10)
pop_label_ALLFULL <- glFULL$other$ind.metrics[, c("id","popdef1", "pop")]

# Attach other information that might be more useful
pop_label_ALLFULL <- pop_label_ALLFULL %>% 
  mutate(snmfcluster = replace(popdef1, popdef1 %in% c("AUS: Gold Coast"), "Group 1")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("AUS: Sydney"), "Group 2")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("NZ: Other"), "Group 3")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Maharashtra subpop. A", "Fiji", "NZ: Napier", "AUS: Melbourne"), "Group 4")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("India: Other"), "Group 5")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Hawaii"), "Group 6")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("South Africa"), "Group 7")) %>%
  mutate(pcacluster = replace(snmfcluster, snmfcluster %in% c("Group 1"), "Cluster 1")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 2", "Group 3", "Group 4"), "Cluster 2")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 5"), "Cluster 3")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 6"), "Cluster 4")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 7"), "Cluster 5")) 

pop_labels <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
                "NZ: Napier", "Fiji", "AUS: Melbourne", 
                "Maharashtra subpop. A", "India: Other", 
                "Hawaii", "South Africa")
pop_abb <- c("GOL", "SYD", "NZO", "NZN", "FIJ", "MEL", 
             "MAH", "IND", "HAW", "SAF")
pop_trans <- data.frame(popdef1 = pop_labels,
                        pop_codes = pop_abb)
pop_label_ALLFULL <- merge(pop_label_ALLFULL, pop_trans)
# Make sure that the ID is sorted properly
pop_label_ALLFULL <- pop_label_ALLFULL[match(indNames(glFULL), pop_label_ALLFULL$id),]

### Figure S8.27: sNMF plot -----------------------------------------
plotQ(q.ALLFULL.ls[2:10],imgoutput="join",
      # clustercol = c("#cab2d6", "#a6cee3",
      #                "#6a3d9a", "#b2df8a",
      #                "yellow", "#fb9a99",
      #                "#ff7f00", "#e31a1c"),
      showdiv = T,
      showindlab=F,
      grplab=pop_label_ALLFULL[, c("pop_codes", "snmfcluster", "pcacluster")],
      selgrp="snmfcluster",
      # subsetgrp=poporder,
      grplabangle = 90,
      grplabheight = 5,
      panelratio = c(4,4),
      grplabsize = 3,
      grplabpos = 0.65,
      splab = c("K=2", "K=3", "K=4", "K=5", "K=6", "K=7", "K=8", "K=9", "K=10"),
      ordergrp=T,
      indlabspacer = 10,
      # showlegend=T,
      # showtitle=T,showsubtitle=T,
      height=2.3,
      width = 16.25,
      barbordercolour=NA,barbordersize=0,
      legendrow = 2,
      legendkeysize = 14,
      legendtextsize = 12,
      legendspacing = 10,
      outputfilename="results/Figs_paper/SI/ALL_FULL_NOTHIN_sNMF_BCF_k2-10",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on ALL dataset (No thinning of SNPs), K = 2-10, 10 repetitions each"
)

# SI: PCA on different MAF thresholds -----------------------------------
## SI S8.6: nmax = 20, MAF > 0.05, thinned ------------------------------
glmaf0_05_thin <- gl.prune.SNP.dist.FROM.VCF(gl = gl.filter.maf(glFULLsub, threshold = 0.05), 
                                             bp_threshold = 100000)

pcmaf0_05_thin <- gl.pcoa(glmaf0_05_thin, nfactors=10, parallel = T, n.cores = 8)
# save(pcmaf0_05_thin, file = "results/PCA/03_ALL/PCA_ALL_FULL_NOTHINnmax20_maf0_05_thin.Rdata")
# load(file = "results/PCA/03_ALL/PCA_ALL_FULL_NOTHINnmax20_maf0_05_thin.Rdata")
### Figure S8.28: scree plot ----------------------------------------
png("results/Figs_paper/SI/ALL_FULL_NOTHINnmax20_maf0_05_thin_PCA_BCF_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
PCA_screeplot(pcmaf0_05_thin, nPC = 15)
dev.off()

# Quick plot
gl.pcoa.plot(pcmaf0_05_thin, glmaf0_05_thin, xaxis = 1, yaxis = 2) 
# Create a table with variance explained
dt_var_expl_FULLsubmaf0_05 <- PCA_var_explained(pcmaf0_05_thin)
round(dt_var_expl_FULLsub$variance_explained, 1)
# Create table with PC scores and populations
dtpcmaf0_05_thin <- data.frame(pcmaf0_05_thin$scores)
dtpcmaf0_05_thin$pop <- factor(glmaf0_05_thin$other$ind.metrics$popdef1,
                               levels = dtcolshapes$pop)
### Figure S8.29: PCA plot ------------------------------------------
pFULLsubmaf0_05_1v2 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 1, yaxis = 2) + PCA_theme

pFULLsubmaf0_05_2v3 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 2, yaxis = 3) + PCA_theme

pFULLsubmaf0_05_3v4 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 3, yaxis = 4) + PCA_theme

pFULLsubmaf0_05_4v5 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 4, yaxis = 5) + PCA_theme

pFULLsubmaf0_05_5v6 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 5, yaxis = 6) + PCA_theme

pFULLsubmaf0_05_6v7 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 6, yaxis = 7) + PCA_theme

pFULLsubmaf0_05_7v8 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 7, yaxis = 8) + PCA_theme

pFULLsubmaf0_05_8v9 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                 dtcolshapes, xaxis = 8, yaxis = 9) + PCA_theme

pFULLsubmaf0_05_9v10 <- plotPCASI(dtpcmaf0_05_thin, varls = round(dt_var_expl_FULLsubmaf0_05$variance_explained, 1),
                                  dtcolshapes, xaxis = 9, yaxis = 10) + PCA_theme

pca_SIFULLsubmaf0_05 <- plot_grid(pFULLsubmaf0_05_1v2 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_2v3 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_3v4 + theme(legend.position = "none"),
                                  pFULLsubmaf0_05_4v5 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_5v6 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_6v7 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_7v8 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_8v9 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_05_9v10 + theme(legend.position = "none"), 
                                  labels = "AUTO", ncol = 3, nrow = 3)

png("results/Figs_paper/SI/ALL_PCA_BCF_FULL_NOTHINnmax20_maf0_05_thin.png", width = 8.3, height = 10, units = "in", res = 600)
plot_grid(pca_SIFULLsubmaf0_05, legend_b2, ncol = 1, rel_heights = c(2, .15))
dev.off()

# PC loadings PCA FULL NOTHIN nmax20 mafthinned
dtpcloadFULLsubmaf0_05 <- data.frame(pcmaf0_05_thin$loadings)
dtpcloadFULLsubmaf0_05$contig <- as.numeric(glmaf0_05_thin$other$loc.metrics$CHROM)
dtpcloadFULLsubmaf0_05$position <- as.numeric(glmaf0_05_thin$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcloadFULLsubmaf0_05 <- dtpcloadFULLsubmaf0_05 %>% arrange(contig, position)

### Figure S8.30: PC loadings ---------------------------------------
png("results/Figs_paper/SI/ALL_FULL_NOTHINnmax20_maf0_05_thin_PCA_BCF_loadings.png", width = 8.3, height = 11.7, units = "in", res = 600)

op <- par(mfrow = c(10,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

for (i in 1:10){
  plot(dtpcloadFULLsubmaf0_05[,i], ylab = paste("PC", i, sep = ""), xaxt = "n", xlab = "", cex = 0.25)
  axis(side = 2, labels = T)
  box(which = "plot", bty = "l")
}

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5)
dev.off()

## sNMF: nmax = 20, MAF > 0.05, thinned -----------------------------
genofileFULLsubmaf0_05 <- "data/processed/LEA/03_ALL/ALL_FULL_nmax20_maf0_05_thin.geno"
# Convert genlight to matrix of 0,1,2 
mat_glmaf0_05_thin <- t(as.matrix(glmaf0_05_thin))
# Write table to .geno file which is essentially a matrix 
write.table(mat_glmaf0_05_thin, file = genofileFULLsubmaf0_05,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)

Krange <- 1:10
snmf_ALLFULLsubmaf0_05 <- snmf(genofileFULLsubmaf0_05, K = Krange, repetitions = 10, ploidy = 2,
                               entropy = T, alpha = 100, project = "new")
# snmf_ALLFULLsubmaf0_05 <- load.snmfProject("data/processed/LEA/03_ALL/ALL_FULL_nmax20_maf0_05_thin.snmfProject")
plot(snmf_ALLFULLsubmaf0_05, col = "blue4", cex = 1.4, pch = 19)
# Calculate mean entropy
crossentropy.matALLFULLsubmaf0_05 <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_ALLFULLsubmaf0_05, K = x)})))
mean.entropyALLFULLsubmaf0_05   <- apply(crossentropy.matALLFULLsubmaf0_05, MARGIN=1, FUN=mean, na.rm=TRUE)
if(any(diff(mean.entropyALLFULLsubmaf0_05)>0)){
  bestK <- unname(which(diff(mean.entropyALLFULLsubmaf0_05)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
bestK

### Figure S8.31: mean cross-entropy plot----------------------------
png("results/Figs_paper/SI/ALL_FULL_nmax20_maf0_05_thin_sNMF_meanxentropy_k1-10.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropyALLFULLsubmaf0_05, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()

q.ALLFULLsubmaf0_05.ls <- merge_sNMF_qmatrix_multiK(snmf_ALLFULLsubmaf0_05, glmaf0_05_thin, k2plot = 1:10)
pop_label_ALLFULLsubmaf0_05 <- glmaf0_05_thin$other$ind.metrics[, c("id","popdef1", "pop")]

# Attach other information that might be more useful
pop_label_ALLFULLsubmaf0_05 <- pop_label_ALLFULLsubmaf0_05 %>% 
  mutate(snmfcluster = replace(popdef1, popdef1 %in% c("AUS: Gold Coast"), "Group 1")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("AUS: Sydney"), "Group 2")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("NZ: Other"), "Group 3")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Maharashtra subpop. A", "Fiji", "NZ: Napier", "AUS: Melbourne"), "Group 4")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("India: Other"), "Group 5")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Hawaii"), "Group 6")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("South Africa"), "Group 7")) %>%
  mutate(pcacluster = replace(snmfcluster, snmfcluster %in% c("Group 1"), "Cluster 1")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 2", "Group 3", "Group 4"), "Cluster 2")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 5"), "Cluster 3")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 6"), "Cluster 4")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 7"), "Cluster 5")) 

pop_labels <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
                "NZ: Napier", "Fiji", "AUS: Melbourne", 
                "Maharashtra subpop. A", "India: Other", 
                "Hawaii", "South Africa")
pop_abb <- c("GOL", "SYD", "NZO", "NZN", "FIJ", "MEL", 
             "MAH", "IND", "HAW", "SAF")
pop_trans <- data.frame(popdef1 = pop_labels,
                        pop_codes = pop_abb)
pop_label_ALLFULLsubmaf0_05 <- merge(pop_label_ALLFULLsubmaf0_05, pop_trans)
# Make sure that the ID is sorted properly
pop_label_ALLFULLsubmaf0_05 <- pop_label_ALLFULLsubmaf0_05[match(indNames(glmaf0_05_thin), pop_label_ALLFULLsubmaf0_05$id),]
### Figure S8.32: sNMF plot -----------------------------------------
plotQ(q.ALLFULLsubmaf0_05.ls[2:10],imgoutput="join",
      # clustercol = c("#cab2d6", "#a6cee3",
      #                "#6a3d9a", "#b2df8a",
      #                "yellow", "#fb9a99",
      #                "#ff7f00", "#e31a1c"),
      showdiv = T,
      showindlab=F,
      grplab=pop_label_ALLFULLsubmaf0_05[, c("pop_codes", "snmfcluster", "pcacluster")],
      selgrp="snmfcluster",
      # subsetgrp=poporder,
      grplabangle = 90,
      grplabheight = 5,
      panelratio = c(4,4),
      grplabsize = 3,
      grplabpos = 0.65,
      splab = c("K=2", "K=3", "K=4", "K=5", "K=6", "K=7", "K=8", "K=9", "K=10"),
      ordergrp=T,
      indlabspacer = 10,
      # showlegend=T,
      # showtitle=T,showsubtitle=T,
      height=2.3,
      width = 16.25,
      barbordercolour=NA,barbordersize=0,
      legendrow = 2,
      legendkeysize = 14,
      legendtextsize = 12,
      legendspacing = 10,
      outputfilename="results/Figs_paper/SI/ALL_FULL_nmax20_maf0_05_thin_sNMF_BCF_k2-10",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on ALL dataset (nmax = 20, MAF > 0.05, thin 100kbp), K = 2-10, 10 repetitions each"
)


## SI S8.7: nmax = 20, MAF > 0.1, thinned ---------------------------
glmaf0_10_thin <- gl.prune.SNP.dist.FROM.VCF(gl = gl.filter.maf(glFULLsub, threshold = 0.1), 
                                             bp_threshold = 100000)

pcmaf0_10_thin <- gl.pcoa(glmaf0_10_thin, nfactors=10, parallel = T, n.cores = 8)

# save(pcmaf0_10_thin, file = "results/PCA/03_ALL/PCA_ALL_FULL_NOTHINnmax20_maf0_10_thin.Rdata")
# load(file = "results/PCA/03_ALL/PCA_ALL_FULL_NOTHINnmax20_maf0_10_thin.Rdata")

### Figure S8.33: screeplot -----------------------------------------
png("results/Figs_paper/SI/ALL_FULL_NOTHINnmax20_maf0_10_thin_PCA_BCF_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
PCA_screeplot(pcmaf0_10_thin, nPC = 15)
dev.off()

# Quick plot
gl.pcoa.plot(pcmaf0_10_thin, glmaf0_10_thin, xaxis = 1, yaxis = 2) 
# Create a table with variance explained
dt_var_expl_FULLsubmaf0_10 <- PCA_var_explained(pcmaf0_10_thin)
round(dt_var_expl_FULLsub$variance_explained, 1)
# Create table with PC scores and populations
dtpcmaf0_10_thin <- data.frame(pcmaf0_10_thin$scores)
dtpcmaf0_10_thin$pop <- factor(glmaf0_10_thin$other$ind.metrics$popdef1,
                               levels = dtcolshapes$pop)
### Figure S8.34: PCA plot ------------------------------------------
pFULLsubmaf0_10_1v2 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 1, yaxis = 2) + PCA_theme

pFULLsubmaf0_10_2v3 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 2, yaxis = 3) + PCA_theme

pFULLsubmaf0_10_3v4 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 3, yaxis = 4) + PCA_theme

pFULLsubmaf0_10_4v5 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 4, yaxis = 5) + PCA_theme

pFULLsubmaf0_10_5v6 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 5, yaxis = 6) + PCA_theme

pFULLsubmaf0_10_6v7 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 6, yaxis = 7) + PCA_theme

pFULLsubmaf0_10_7v8 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 7, yaxis = 8) + PCA_theme

pFULLsubmaf0_10_8v9 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                 dtcolshapes, xaxis = 8, yaxis = 9) + PCA_theme

pFULLsubmaf0_10_9v10 <- plotPCASI(dtpcmaf0_10_thin, varls = round(dt_var_expl_FULLsubmaf0_10$variance_explained, 1),
                                  dtcolshapes, xaxis = 9, yaxis = 10) + PCA_theme

pca_SIFULLsubmaf0_10 <- plot_grid(pFULLsubmaf0_10_1v2 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_2v3 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_3v4 + theme(legend.position = "none"),
                                  pFULLsubmaf0_10_4v5 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_5v6 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_6v7 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_7v8 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_8v9 + theme(legend.position = "none"), 
                                  pFULLsubmaf0_10_9v10 + theme(legend.position = "none"), 
                                  labels = "AUTO", ncol = 3, nrow = 3)

png("results/Figs_paper/SI/ALL_PCA_BCF_FULL_NOTHINnmax20_maf0_10_thin.png", width = 8.3, height = 10, units = "in", res = 600)
plot_grid(pca_SIFULLsubmaf0_10, legend_b2, ncol = 1, rel_heights = c(2, .15))
dev.off()

# PC loadings PCA FULL NOTHIN nmax20 mafthinned
dtpcloadFULLsubmaf0_10 <- data.frame(pcmaf0_10_thin$loadings)
dtpcloadFULLsubmaf0_10$contig <- as.numeric(glmaf0_10_thin$other$loc.metrics$CHROM)
dtpcloadFULLsubmaf0_10$position <- as.numeric(glmaf0_10_thin$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcloadFULLsubmaf0_10 <- dtpcloadFULLsubmaf0_10 %>% arrange(contig, position)

### Figure S8.35: PC loadings ---------------------------------------
png("results/Figs_paper/SI/ALL_FULL_NOTHINnmax20_maf0_10_thin_PCA_BCF_loadings.png", width = 8.3, height = 11.7, units = "in", res = 600)

op <- par(mfrow = c(10,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

for (i in 1:10){
  plot(dtpcloadFULLsubmaf0_10[,i], ylab = paste("PC", i, sep = ""), xaxt = "n", xlab = "", cex = 0.25)
  axis(side = 2, labels = T)
  box(which = "plot", bty = "l")
}

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5)
dev.off()

## sNMF: nmax = 20, MAF > 0.1, thinned ---------------------------------------
genofileFULLsubmaf0_10 <- "data/processed/LEA/03_ALL/ALL_FULL_nmax20_maf0_10_thin.geno"
# Convert genlight to matrix of 0,1,2 
mat_glmaf0_10_thin <- t(as.matrix(glmaf0_10_thin))
# Write table to .geno file which is essentially a matrix 
write.table(mat_glmaf0_10_thin, file = genofileFULLsubmaf0_10,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)

Krange <- 1:10
snmf_ALLFULLsubmaf0_10 <- snmf(genofileFULLsubmaf0_10, K = Krange, repetitions = 10, ploidy = 2,
                               entropy = T, alpha = 100, project = "new")
snmf_ALLFULLsubmaf0_10 <- load.snmfProject("data/processed/LEA/03_ALL/ALL_FULL_nmax20_maf0_10_thin.snmfProject")
plot(snmf_ALLFULLsubmaf0_10, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), ALL dataset, after HWE filter, only relevant pops.")
# Calculate mean entropy
crossentropy.matALLFULLsubmaf0_10 <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_ALLFULLsubmaf0_10, K = x)})))
mean.entropyALLFULLsubmaf0_10   <- apply(crossentropy.matALLFULLsubmaf0_10, MARGIN=1, FUN=mean, na.rm=TRUE)
if(any(diff(mean.entropyALLFULLsubmaf0_10)>0)){
  bestK <- unname(which(diff(mean.entropyALLFULLsubmaf0_10)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
bestK

### Figure S8.36: mean cross-entropy plot----------------------------
png("results/Figs_paper/SI/ALL_FULL_nmax20_maf0_10_thin_sNMF_meanxentropy_k1-10.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropyALLFULLsubmaf0_10, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()

q.ALLFULLsubmaf0_10.ls <- merge_sNMF_qmatrix_multiK(snmf_ALLFULLsubmaf0_10, glmaf0_10_thin, k2plot = 1:10)
pop_label_ALLFULLsubmaf0_10 <- glmaf0_10_thin$other$ind.metrics[, c("id","popdef1", "pop")]

# Attach other information that might be more useful
pop_label_ALLFULLsubmaf0_10 <- pop_label_ALLFULLsubmaf0_10 %>% 
  mutate(snmfcluster = replace(popdef1, popdef1 %in% c("AUS: Gold Coast"), "Group 1")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("AUS: Sydney"), "Group 2")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("NZ: Other"), "Group 3")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Maharashtra subpop. A", "Fiji", "NZ: Napier", "AUS: Melbourne"), "Group 4")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("India: Other"), "Group 5")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Hawaii"), "Group 6")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("South Africa"), "Group 7")) %>%
  mutate(pcacluster = replace(snmfcluster, snmfcluster %in% c("Group 1"), "Cluster 1")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 2", "Group 3", "Group 4"), "Cluster 2")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 5"), "Cluster 3")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 6"), "Cluster 4")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 7"), "Cluster 5")) 

pop_labels <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
                "NZ: Napier", "Fiji", "AUS: Melbourne", 
                "Maharashtra subpop. A", "India: Other", 
                "Hawaii", "South Africa")
pop_abb <- c("GOL", "SYD", "NZO", "NZN", "FIJ", "MEL", 
             "MAH", "IND", "HAW", "SAF")
pop_trans <- data.frame(popdef1 = pop_labels,
                        pop_codes = pop_abb)
pop_label_ALLFULLsubmaf0_10 <- merge(pop_label_ALLFULLsubmaf0_10, pop_trans)
# Make sure that the ID is sorted properly
pop_label_ALLFULLsubmaf0_10 <- pop_label_ALLFULLsubmaf0_10[match(indNames(glmaf0_10_thin), pop_label_ALLFULLsubmaf0_10$id),]
### Figure S8.37: sNMF plot -----------------------------------------
plotQ(q.ALLFULLsubmaf0_10.ls[2:10],imgoutput="join",
      # clustercol = c("#cab2d6", "#a6cee3",
      #                "#6a3d9a", "#b2df8a",
      #                "yellow", "#fb9a99",
      #                "#ff7f00", "#e31a1c"),
      showdiv = T,
      showindlab=F,
      grplab=pop_label_ALLFULLsubmaf0_10[, c("pop_codes", "snmfcluster", "pcacluster")],
      selgrp="snmfcluster",
      # subsetgrp=poporder,
      grplabangle = 90,
      grplabheight = 5,
      panelratio = c(4,4),
      grplabsize = 3,
      grplabpos = 0.65,
      splab = c("K=2", "K=3", "K=4", "K=5", "K=6", "K=7", "K=8", "K=9", "K=10"),
      ordergrp=T,
      indlabspacer = 10,
      # showlegend=T,
      # showtitle=T,showsubtitle=T,
      height=2.3,
      width = 16.25,
      barbordercolour=NA,barbordersize=0,
      legendrow = 2,
      legendkeysize = 14,
      legendtextsize = 12,
      legendspacing = 10,
      outputfilename="results/Figs_paper/SI/ALL_FULL_nmax20_maf0_10_thin_sNMF_BCF_k2-10",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on ALL dataset (nmax = 20, maf > 0.1, thin 100kbp), K = 2-10, 10 repetitions each"
)

