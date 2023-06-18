# IND_popstruc_SI_DArT.R --------------------------------------------
# This script is to show that results from the DArT dataset
# is comparable to BCFtools datasets. Appendix S10.5
# 
# Define input ------------------------------------------------------
dartcsv <- "../01_download_data/DArT/DArT_IND/Report_DImy20-5737_SNP_2.csv"
metacsv <- '../01_download_data/DArT/metadata_all_samples.csv'
# Define location to store .geno file for sNMF analysis -------------
genofile <- "data/processed/LEA/02_IND/IND_DART_thinnedonly.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/02_IND/", recursive = T, showWarnings = F)
dir.create("results/PCA/02_IND/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
# Load libraries ----------------------------------------------------
library(tidyverse)
library(dartR)
library(cowplot)
library(LEA) # Version 3.2.0. Some function names change in newer versions.
library(pophelper)
source("../shared_scripts/functions.R")
# Read in data ------------------------------------------------------
gl_IND <- gl.read.dart(filename = dartcsv,
                      ind.metafile = metacsv)
# Define Maharashtra subpop. A within Indian data -------------------
gl_IND$other$ind.metrics$IND_pop_recode <- as.character(gl_IND$other$ind.metrics$popdef2)
gl_IND$other$ind.metrics <- gl_IND$other$ind.metrics %>%
  mutate(IND_pop_recode = replace(IND_pop_recode, id %in% paste("M0", seq(210, 215), sep = ""), "Maharashtra subpop. A"))
# This can also be done by mutating the values from popdef1, but this 
# is based on the ID numbers instead.
# 
# This was only included for making the plots for the sNMF
# population structure plot 
#
pop(gl_IND) <- gl_IND@other$ind.metrics$popdef2
# save(gl_IND, file = "data/gl_IND.Rdata")
# load("data/gl_IND.Rdata")
# Remove replicats --------------------------------------------------
# Kyle replicates
ind_rep <- grep("^[0-9]*[a-zA-Z]$", indNames(gl_IND), value = T)
# Append our in-house replicates
ind_rep <- c(ind_rep, "M0092", "M0093", "M0094", "M0186", "M0187", "M0370", "M0188", "M0280", "M0281", "M0282", "M0368", "M0369", "M0237")
# Drop replicates 
gl_IND_sub1 <- gl.drop.ind(gl_IND, ind.list = ind_rep, recalc = T)

# Report heterozygousity after removal of replicates ----------------
png("results//DArT_Ho_indiv.png", height = 5, width = 8, units = "in", res = 600)
dt_IND_sub1 <- gl.report.heterozygosity(gl_IND_sub1, method = "ind")
dev.off()

# Sample M0271 should have been removed but was not in the DArT 
# dataset. This was removed in the ALL_DArT dataset and also in
# STACKS and BCFtools ALL and IND datasets. This is a mistake but 
# did not alter the results
#
# gl_IND_sub1 <- gl.filter.heterozygosity(gl_IND_sub1, t.upper = 0.175)


# Filter for loci that has only been called once on the reference ---
loc_alncnt_not1 <- locNames(gl_IND_sub1)[gl_IND_sub1$other$loc.metrics$AlnCnt_Indian_Myna_v2.1 != 1]
gl_IND_sub2 <- gl.drop.loc(gl_IND_sub1, loc.list = loc_alncnt_not1)

# Filter loci based on call rates and reproducibility ---------------
# This also removes single and doubletons. 
gl_IND_sub3 <- gl_filt_CR_rep_MAC(gl = gl_IND_sub2, CR = 0.8, 
                                 rep = 0.95, MAC_threshold = 1, 
                                 SNP_in_single_samp_rm = T)

gl_IND_sub4 <- gl.prune.SNP.dist(gl_IND_sub3, bp_threshold = 100000)

# Count samples per location ----------------------------------------
gl_IND_sub4@other$ind.metrics %>% dplyr::count(popdef2)
pop(gl_IND_sub4) <- gl_IND_sub4$other$ind.metrics$popdef2
# Re-perform PCA on the full dataset --------------------------------
pc <- gl.pcoa(gl_IND_sub4, nfactors=10, parallel = T, n.cores = 8)
# Quick plot to see variance explained ------------------------------
gl.pcoa.plot(pc, gl_IND_sub4, pop.labels= "pop") 
# Create dataframe for pricipal component axes ----------------------
dtpc <- data.frame(pc$scores)
save(pc, file = "results/PCA/02_IND/PCA_IND_DART_thinnedonly.Rdata")
load("results/PCA/02_IND/PCA_IND_DART_thinnedonly.Rdata")
# Assign population information -------------------------------------
dtpc$popdef2 <- gl_IND_sub4$other$ind.metrics$popdef2
dtpc$IND_pop_recode <- gl_IND_sub4$other$ind.metrics$IND_pop_recode
write.csv(dtpc, file = "results/PCA/02_IND/PCA_IND_DART_thinnedonly.csv")
dtpc <- read.csv("results/PCA/02_IND/PCA_IND_DART_thinnedonly.csv")

# Get lat and long for each location to sort legends by
dtlatlon <- gl_IND_sub4@other$ind.metrics %>% 
  dplyr::group_by(popdef2) %>%
  dplyr::summarise(lat = mean(latitude),
            lon = mean(longitude))
dtlatlon$popdef2[order(dtlatlon$lon)]
lab_order <- dtlatlon$popdef2[order(dtlatlon$lon)]
shape_order <- c(16, 11, 17, 15, 3, 12, 10, 8, 5)
col_order <- gg_color_hue(9)
dtcolshapes <- data.frame(pop = as.character(lab_order),
                          pop_lab = as.character(lab_order),
                          col = col_order,
                          shape = shape_order)
# Create dataframe with the variance explained for each PC ---------
dtvar <- PCA_var_explained(pc)
round(dtvar$variance_explained, 1)
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
pl <- ggplot() + 
  geom_point(data = dtpc, 
             aes(x = PC1, y = PC2, 
                 shape = popdef2, 
                 col = popdef2)) + xlab("PCA 1 (3%)") + ylab ("PCA 2 (2.3%)") + theme_classic() +
  scale_shape_manual(values = shape_order) +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, face = "bold"),
        legend.key.size = unit(0.75, 'lines'),
        legend.box.margin = margin(-10,0,-10,-10),
        legend.position = c(0.05,0.975),
        legend.justification = c(0,1),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_rect(colour = NA, fill='transparent'),
        legend.margin = margin(t = -0.2, r = 0.1, b = 0, l = 0, unit = "cm"),
        legend.background = element_rect(color = "black", size = 0.01, linetype = "solid"))
pl1 <- pl 
pl1

# Make screeplot for SI -------------------------------------------------------
screep <- PCA_screeplot(pc, nPC = 15)
png("results/Figs_paper/SI/IND_PCA_DART_thinnedonly_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
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
                    pc4v5 + theme(legend.position = "none"), labels = "AUTO", ncol = 2, nrow = 2)
png("results/Figs_paper/SI/IND_PCA_DART_thinnedonly.png", width = 8.3, height = 8.5, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(2, .15))
dev.off()

## PC loadings -----------------------------------------------------
dtpcload <- data.frame(pc$loadings)
dtpcload$contig <- as.numeric(gl_IND_sub4$other$loc.metrics$Chrom_Indian_Myna_v2.1)
dtpcload$position <- as.numeric(gl_IND_sub4$other$loc.metrics$ChromPos_Indian_Myna_v2.1)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/IND_PCA_DART_thinnedonly_loadings.png", width = 8.3, height = 8.4, units = "in", res = 600)

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
mat_gl <- t(as.matrix(gl_IND_sub4))
# Write table to .geno file which is essentially a matrix -----------
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
Krange <- 1:5
snmf_IND_thinnedonly <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2, 
                             entropy = T, alpha = 100, project = "new",)
snmf_IND_thinnedonly <- load.snmfProject("data/processed/LEA/02_IND/IND_DART_thinnedonly.snmfProject")
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
png("results/Figs_paper/SI/IND_sNMF_meanxentropy_DART_thinnedonly_k1-5.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropyIND, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()


q.IND.ls <- merge_sNMF_qmatrix_multiK(snmf_IND_thinnedonly, gl_IND_sub4, k2plot = 1:5)
pop_label_IND <- gl_IND_sub4$other$ind.metrics[, c("id","IND_pop_recode")]
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
      outputfilename="results/Figs_paper/SI/IND_sNMF_DART_thinnedonly_k2-5",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on IND dataset, K = 2-5, 10 repetitions each"
)
