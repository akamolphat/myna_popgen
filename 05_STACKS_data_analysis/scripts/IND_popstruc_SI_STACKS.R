# This script makes plots for the figure for publication for the 
# IND dataset
#
# This script is for the STACKS dataset
#
# This is based on the dataset prior to HWE filtering
#
# This is just to show that the results are relatively
# consistent with BCFtools 
#
# Define input files
#
# Define input ------------------------------------------------------
vcffile <- "../01_download_data/STACKS_IND/populations.IND.bialminGQ30DP15-100.norep.highnegfis.lmiss20.nosingledoubletons.vcfthin.snps.vcf.gz"
metadtfile <- '../01_download_data/TableS1.2.csv'
# Define location to store .geno file for sNMF analysis -------------
genofile <- "data/processed/LEA/02_IND/IND_STACKS_thinnedonly.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/02_IND/", recursive = T, showWarnings = F)
dir.create("results/PCA/02_IND/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
# Load libraries ----------------------------------------------------
library(dartR)
library(tidyverse)
library(cowplot)
# library(ggforce)
library(xlsx)
# library(pheatmap)
library(pophelper)
library(LEA)
source("../shared_scripts/functions.R")
# Read in VCF file --------------------------------------------------
gl <- gl.read.vcf(vcffile)
# Get fixed section from VCF for contig and position ----------------
vcf_in <- read.vcfR(vcffile, verbose = FALSE )
vcfFIX <- getFIX(vcf_in)
vcfFIXmerged <- cbind(vcfFIX, gl$other$loc.metrics)
gl$other$loc.metrics <- vcfFIXmerged
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
# This can also be done by mutating the values from popdef1, but this 
# is based on the ID numbers instead.
# 
# This was only included for making the plots for the sNMF
# population structure plot 
#
pop(gl) <- gl@other$ind.metrics$popdef2
# Count samples per location ----------------------------------------
gl@other$ind.metrics %>% dplyr::count(popdef2)
# Re-perform PCA on the full dataset --------------------------------
pc <- gl.pcoa(gl, nfactors=10, parallel = T, n.cores = 8)
# Save pc object for reproducibility
save(pc, file = "results/PCA/02_IND/PCA_IND_STACKS_thinnedonly.Rdata")
load("results/PCA/02_IND/PCA_IND_STACKS_thinnedonly.Rdata")
# Quick plot to see variance explained ------------------------------
gl.pcoa.plot(pc, gl, pop.labels= "pop") 
# Create dataframe for pricipal component axes ----------------------
dtpc <- data.frame(pc$scores)
# Assign population information -------------------------------------
dtpc$popdef2 <- gl$other$ind.metrics$popdef2
dtpc$IND_pop_recode <- gl$other$ind.metrics$IND_pop_recode
write.csv(dtpc, file = "results/PCA/02_IND/PCA_IND_STACKS_thinnedonly.csv")
dtpc <- read.csv("results/PCA/02_IND/PCA_IND_STACKS_thinnedonly.csv")
# Get lat and long for each location to sort legends by
dtlatlon <- gl@other$ind.metrics %>% 
  group_by(popdef2) %>%
  summarise(lat = mean(latitude),
            lon = mean(longitude))
dtlatlon$popdef2[order(dtlatlon$lon)]
# lab_order <- dtlatlon$popdef2[order(dtlatlon$lon)]
lab_order <- dtlatlon$popdef2
shape_order <- c(16, 11, 17, 15, 3, 12, 10, 8, 5)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_order <- gg_color_hue(9)
dtcolshapes <- data.frame(pop = as.character(lab_order),
                          pop_lab = as.character(lab_order),
                          col = col_order,
                          shape = shape_order)
# Create dataframe with the variance explained for each PC ---------
dtvar <- PCA_var_explained(pc)
# Converrt dtpc population column to factor and in order ------------
dtpc$popdef2 <- factor(dtpc$popdef2,
                       levels = lab_order)
# Create PC1 vs PC2 plot --------------------------------------------
pl <- plotPCASI(dtpc, 
                varls = round(dtvar$variance_explained, 1),
                dtcolshapes,
                xaxis = 1,
                yaxis = 2, popcol = "popdef2") + 
  PCA_theme

pl
pl1 <- pl

# Make screeplot for SI -------------------------------------------------------
screep <- PCA_screeplot(pc, nPC = 15)
png("results/Figs_paper/SI/IND_PCA_STACKS_thinnedonly_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep
dev.off()

# SI plots ----------------------------------------------------------
## PCA --------------------------------------------------------------
pc2v3 <- plotPCASI(dtpc, 
                   varls = round(dtvar$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 2,
                   yaxis = 3, popcol = "popdef2") + 
  PCA_theme

pc3v4 <- plotPCASI(dtpc, 
                   varls = round(dtvar$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 3,
                   yaxis = 4, popcol = "popdef2") + 
  PCA_theme

pc4v5 <- plotPCASI(dtpc, 
                   varls = round(dtvar$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 4,
                   yaxis = 5, popcol = "popdef2") + 
  PCA_theme

legend_b <- get_legend(
  pc2v3
)

pca_SI <- plot_grid(pl + theme(legend.position = "none"), 
                    pc2v3 + theme(legend.position = "none"), 
                    pc3v4 + theme(legend.position = "none"),
                    pc4v5 + theme(legend.position = "none"), 
                    labels = "AUTO", ncol = 2, nrow = 2)
png("results/Figs_paper/SI/IND_PCA_STACKS_thinnedonly.png", width = 8.3, height = 8.5, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(2, .15))
dev.off()

## PC loadings -----------------------------------------------------
dtpcload <- data.frame(pc$loadings)
dtpcload$contig <- as.numeric(gl$other$loc.metrics$CHROM)
dtpcload$position <- as.numeric(gl$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/IND_PCA_STACKS_thinnedonly_loadings.png", width = 8.3, height = 8.4, units = "in", res = 600)

op <- par(mfrow = c(5,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

for (i in 1:5){
  plot(dtpcload[,i], ylab = paste("PC", i, sep = ""), xaxt = "n", xlab = "", cex = 0.25)
  axis(side = 2, labels = T)
  box(which = "plot", bty = "l")
}

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5)
dev.off()

## sNMF -------------------------------------------------------------
# Plot cross entropy plot
# Plot k=2-5
# Convert genlight to matrix of 0,1,2 -------------------------------
mat_gl <- t(as.matrix(gl))
# Write table to .geno file which is essentially a matrix -----------
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
Krange <- 1:5
snmf_IND_thinnedonly <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2, 
                     entropy = T, alpha = 100, project = "new",)
# snmf_IND_thinnedonly <- load.snmfProject("data/processed/LEA/02_IND/IND_STACKS_thinnedonly.snmfProject")
plot(snmf_IND_thinnedonly, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), IND dataset.")
# Calculate mean entropy
crossentropy.matIND <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_IND_thinnedonly, K = x)})))
mean.entropyIND   <- apply(crossentropy.matIND, MARGIN=1, FUN=mean, na.rm=TRUE)
if(any(diff(mean.entropyIND)>0)){
  bestK <- unname(which(diff(mean.entropyIND)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
bestK
# SI: Export mean cross-entropy plot --------------------------------
png("results/Figs_paper/SI/IND_sNMF_meanxentropy_STACKS_thinnedonly_k1-5.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropyIND, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()


q.IND.ls <- merge_sNMF_qmatrix_multiK(snmf_IND_thinnedonly, gl, k2plot = 1:5)
pop_label_IND <- gl$other$ind.metrics[, c("ID","IND_pop_recode")]
dtlatlon$IND_pop_recode[order(dtlatlon$lon)]
poporder <- c("Gujarat", "Maharashtra subpop. A", "Maharashtra", "Madhya Pradesh", "Karnataka",
              "Andhra Pradesh", "Tamil Nadu", 
              "Uttar Pradesh", "Odisha", "West Bengal")
plotQ(q.IND.ls[2:5],imgoutput="join",
      showdiv = T,
      showindlab=F,
      grplab=pop_label_IND[, c("pop"), drop = F],
      subsetgrp=poporder,
      grplabangle = 90,
      grplabheight = 5,
      panelratio = c(4,9),
      grplabsize = 3,
      grplabpos = 0.65,
      splab = c("K=2", "K=3", "K=4", "K=5"),
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
      outputfilename="results/Figs_paper/SI/IND_sNMF_STACKS_thinnedonly_k1-5",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on IND dataset, K = 2-5, 10 repetitions each"
)
