# This script makes plots for the figure for publication for the ALL dataset
#
# This script is for the BCFtools dataset
#
# Define input files ------------------------------------------------
#
vcffile <- "../01_download_data/BCFtools_ALL/variant_calls.ALL.bialminQ30minGQ30DP15-125.norep.noadm.highnegfis.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2v2.csv"
# Define location to store .geno filefor sNMF analysis --------------
genofile <- "data/processed/LEA/03_ALL/snmf_ALL_hwe.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/PCA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
# Load libraries ----------------------------------------------------
library(dartR)
library(LEA)
library(tidyverse)
library(cowplot)
library(xlsx)
library(pophelper)
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
# Attach population information to genlight object
pop(gl) <- gl$other$ind.metrics$popdef2
# Count samples per location ----------------------------------------
gl@other$ind.metrics %>% dplyr::count(popdef2)
# Re-perform PCA on the ALL dataset, nmax = 20 ----------------------
dtsub <- read.csv("data/processed/ALL_BCFtools_subset_nmax20.csv")
gl_sub <- gl.keep.ind(gl, ind.list = dtsub$subsampled_ID, recalc = T, mono.rm = T)
pc_sub <- gl.pcoa(gl_sub, nfactors=10, parallel = T, n.cores = 8)
# Save pc object to .Rdata for later access
save(pc_sub, file = "results/PCA/03_ALL/PCA_ALL_hwe_nmax20.Rdata")
load("results/PCA/03_ALL/PCA_ALL_hwe_nmax20.Rdata")
# Quick plot to see variance explained
gl.pcoa.plot(pc_sub, gl_sub, pop.labels= "pop") 
# Create dataframe of PC scores
dtpc <- data.frame(pc_sub$scores)
dtpc$pop <- pop(gl_sub)
write.csv(dtpc, file = "results/PCA/03_ALL/PCA_ALL_hwe_nmax20.csv")
dtpc <- read.csv("results/PCA/03_ALL/PCA_ALL_hwe_nmax20.csv")
## Make tables for colours ------------------------------------------
lab_order <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
               "NZ: Napier", "Fiji", "AUS: Melbourne", 
               "IND: Maharashtra subpopulation A", "IND: Other", 
               "Hawaii", "South Africa")
lab_order_abb <- c("AUS: Gold Coast (GOL)", "AUS: Sydney (SYD)", "NZ: Other (NZO)",
                   "NZ: Napier (NZN)", "Fiji (FIJ)", "AUS: Melbourne (MEL)",
                   "IND: Maharashtra subpop. A (MAH)", "IND: Other (IND)",
                   "Hawaii (HAW)", "South Africa (SAF)")

dtpc$pop <- factor(dtpc$pop,
               levels = lab_order)
dtcolshapes <- data.frame(pop = lab_order,
                          pop_lab = lab_order_abb,
                          col = c("#a6cee3", "#b2df8a", "#cab2d6", 
                                  "#ff7f00", "#33a02c", "#1f78b4",
                                  "#e31a1c", "#fdbf6f", "#fb9a99",
                                  "#6a3d9a"),
                          shape = c(16, 3, 12, 18, 8, 14, 11, 17,15,7))

# Main text plots ---------------------------------------------------
## Make PCA plots ---------------------------------------------------
dt_var_expl <- PCA_var_explained(pc_sub)
### Define overall theme for PCA plots 
PCA_theme <- theme_classic() + theme(axis.title = element_text(size = 12, face = "bold"),
                   axis.text = element_text(size = 10, face = "bold"),
                   legend.title = element_blank(),
                   legend.key.size = unit(0.75, 'lines'),
                   legend.text = element_text(size = 8, face = "bold"),
                   legend.background = element_rect(color = "black", size = 0.01, linetype = "solid"),
                   legend.position = "bottom",
                   aspect.ratio=1)
### Make PC1 vs PC2 plot, ALL dataset -------------------------------
pl <- plotPCASI(dtpc, 
                varls = round(dt_var_expl$variance_explained, 1),
                dtcolshapes,
                xaxis = 1,
                yaxis = 2) + 
  PCA_theme +
  guides(colour=guide_legend(nrow=2,byrow=T), 
         shape=guide_legend(nrow=2,byrow=T)) 
pl
### Make PC1 vs PC2 plot, relavent pops -----------------------------
gl_sub2 <- gl.keep.pop(gl_sub, 
                       pop.list = c("AUS: Sydney", "AUS: Melbourne", "NZ: Napier", "NZ: Other", "Fiji", "IND: Maharashtra subpopulation A"), 
                       recalc = T, mono.rm = T)
pc_sub2 <- gl.pcoa(gl_sub2, nfactors=10, parallel = T, n.cores = 8)
save(pc_sub2, file = "results/PCA/03_ALL/PCA_ALLrelaventpops_hwe_nmax20.Rdata")
load("results/PCA/03_ALL/PCA_ALLrelaventpops_hwe_nmax20.Rdata")
gl.pcoa.plot(pc_sub2, gl_sub2, pop.labels= "pop") 
dtpc2 <- data.frame(pc_sub2$scores)
dtpc2$pop <- pop(gl_sub2)
write.csv(dtpc2, file = "results/PCA/03_ALL/PCA_ALLrelaventpops_hwe_nmax20.csv")
dtpc2 <- read.csv("results/PCA/03_ALL/PCA_ALLrelaventpops_hwe_nmax20.csv")
dtpc2$pop <- factor(dtpc2$pop,
                    levels = lab_order)

dt_var_expl2 <- PCA_var_explained(pc_sub2)

pl2 <- plotPCASI(dtpc2, 
                 varls = round(dt_var_expl2$variance_explained, 1),
                 dtcolshapes,
                 xaxis = 1,
                 yaxis = 2) + 
  PCA_theme
### Merge plots into panel ------------------------------------------
prow <- plot_grid(pl + theme(legend.position = "none"), 
                  pl2 + theme(legend.position = "none"), ncol = 2, rel_widths = c(4,4), labels = NULL)
legend_b <- get_legend(
  pl +
    guides(color = guide_legend(nrow = 4),
           shape = guide_legend(nrow = 4)) +
    theme(legend.position = "bottom",
          legend.key.size = unit(1, 'lines'),
          legend.text = element_text(size = 12, face = "bold")
    )
)
#### Add legend 
prow2 <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .25))
# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
png("results/Figs_paper/Main/ALL_pop_str_PCA.png", width = 8.3, height = 5, units = "in", res = 600)
prow2
dev.off()

# sNMF --------------------------------------------------------------
# Previous iterations of the analysis found that optimal K 
# is likely ~ 7 or 8
# Convert genlight to matrix of 0,1,2 -------------------------------
mat_gl <- t(as.matrix(gl))
# Write table to .geno file which is essentially a matrix -----------
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
Krange <- 1:20
snmf_ALL <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2,
                    entropy = T, alpha = 100, project = "new")
# snmf_ALL <- load.snmfProject("data/processed/LEA/03_ALL/snmf_ALL_hwe.snmfProject")

plot(snmf_ALL, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), ALL, after HWE.")

# Calculate mean entropy
crossentropy.mat <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_ALL, K = x)})))
mean.entropy   <- apply(crossentropy.mat, MARGIN=1, FUN=mean, na.rm=TRUE)
min.entropy   <- apply(crossentropy.mat, MARGIN=1, FUN=min, na.rm=TRUE)

# Find k value with lowest mean entropy
if(any(diff(mean.entropy)>0)){
  bestK <- unname(which(diff(mean.entropy)>0)[1])
} else {
  bestK <- unname(Krange[1])
}

# SI: Export cross-entropy plot --------------------------------

png("results/Figs_paper/SI/ALL_sNMF_min_meanxentropy_k1-20.png", width = 8.3, height = 7.5, units = "in", res = 600)
op <- par(mfrow = c(2,1),
          oma = c(4,0,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)
plot(mean.entropy, xaxt = "n", xlab = "", ylab = "Mean cross-entropy")
plot(min.entropy, xlab = "", ylab = "Minimum cross-entropy")
title(xlab = "Number of ancestral populations (K)",
      outer = TRUE, line = 2,
      cex.lab = 1)
dev.off()


rep <- 1:10
# Create labels for samples 
# Note that the gl object here is the same object used to run sNMF that was
# called earlier.
samplenames <- indNames(gl)
numind <- length(samplenames)
# Make table of population information
pop_label_sub <- gl$other$ind.metrics[, c("ID","popdef2", "popdef1")]
# Attach other information that might be more useful
pop_label_sub <- pop_label_sub %>% 
  mutate(snmfcluster = replace(popdef2, popdef2 %in% c("AUS: Gold Coast"), "1")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef2 %in% c("AUS: Sydney"), "2A")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef2 %in% c("NZ: Other"), "2B")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef2 %in% c("IND: Maharashtra subpopulation A", "Fiji", "NZ: Napier", "AUS: Melbourne"), "2C")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef2 %in% c("IND: Other"), "3")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef2 %in% c("Hawaii"), "4")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef2 %in% c("South Africa"), "5")) %>%
  mutate(pcacluster = replace(snmfcluster, snmfcluster %in% c("Group 1"), "Cluster 1")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 2", "Group 3", "Group 4"), "Cluster 2")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 5"), "Cluster 3")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 6"), "Cluster 4")) %>%
  mutate(pcacluster = replace(pcacluster, snmfcluster %in% c("Group 7"), "Cluster 5")) 


pop_labels <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
                "NZ: Napier", "Fiji", "AUS: Melbourne", 
                "IND: Maharashtra subpopulation A", "IND: Other", 
                "Hawaii", "South Africa")
pop_abb <- c("GOL", "SYD", "NZO", "NZN", "FIJ", "MEL", 
             "MAH", "IND", "HAW", "SAF")
pop_trans <- data.frame(popdef2 = pop_labels,
                        pop_codes = pop_abb)
pop_label_sub <- merge(pop_label_sub, pop_trans)
# Make sure that the ID is sorted properly
pop_label_sub <- pop_label_sub[match(indNames(gl), pop_label_sub$ID),]

# Create list of qmatrix from the same K number 
q.merged <- merge_sNMF_qmatrix_multiK(snmf_obj = snmf_ALL,
                                      gl,
                                      k2plot = 8)
pp <- plotQ(q.merged,
            clustercol = c("#cab2d6", "#a6cee3",
                           "#6a3d9a", "#b2df8a",
                           "yellow", "#fb9a99",
                           "#ff7f00", "#e31a1c"),
            # imgoutput="join",
            showdiv = T,
            showindlab=F,
            grplab=pop_label_sub[, c("pop_codes", "snmfcluster")],
            # subsetgrp=c("NZ: Other", "NZ: Napier"),
            grplabangle = 90,
            grplabheight = 3,
            grplabsize = 4,
            grplabpos = 0.65,
            grplabcol = "black",
            selgrp="snmfcluster",
            splab = c("K=8"),
            # panelratio = c(1,1),
            # grplabspacer = -2,
            ordergrp=T,
            showlegend=F,
            showtitle=F,
            showsubtitle=F,
            height=4,
            width =21,
            panelratio = c(9,6),
            linesize = 0.5,
            # indlabspacer=-1,
            barbordercolour=NA,barbordersize=0,
            legendrow = 1,
            showsp = F,
            # legendlab = c("K1", "K2"),
            outputfilename="results/Figs_paper/Main/ALL_BCFtools_hwe_sNMF_k8_merged",
            imgtype="png",
            exportpath=getwd(),
            titlelab="Common myna population structure",
            returnplot = T#,
            # subtitlelab="sNMF on subsetted dataset, n = 10 per population."
)

prow3 <- plot_grid(prow2, pp$plot[[1]], nrow = 2, 
                   rel_heights = c(5.44, 3), 
                   labels = NULL)
## Figure 4: Output to png ------------------------------------------
# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
png("results/Figs_paper/Main/ALL_pop_str.png", width = 8.3, height = 8.44, units = "in", res = 600)
prow3
dev.off()

## Figure 5: FST matrix popdef2 -------------------------------------
FSTpoporder <- c("South Africa", "Hawaii", "AUS: Gold Coast",
                 "AUS: Sydney", "NZ: Other", 
                 "AUS: Melbourne", "NZ: Napier", "Fiji",
                 "IND: Maharashtra subpopulation A", "IND: Other")
pFST <- plot_lowertri_FST_mat_from_xlsx(FSTxlsx = "results/FST/03_ALL/pairwiseFST_ALL_BCF_hwe_100bs.xlsx", 
                                        poporder = FSTpoporder,
                                        dt1spreadsheetname = "FST_100bs",
                                        dt2spreadsheetname = "pval_100bs")
# Abbreviate IND: Maharashtra subpopulation A
FSTpoprename <- FSTpoporder
FSTpoprename[9] <- "IND: Maharashtra subpop. A"

pFST <- pFST + scale_fill_distiller(palette = "Spectral",
                            direction = 1,
                            limits = c(-0.001, 0.32),
                            na.value = "white") + 
  scale_y_discrete(breaks = FSTpoporder,
                   labels = FSTpoprename,
                   # position = "left",
                   position = "right",
                   expand = c(0,0)) +
  scale_x_discrete(breaks = FSTpoporder,
                   labels = FSTpoprename,
                   limits = rev,
                   expand = c(0,0)) +
  guides(fill = guide_colourbar(title = expression(~"Population pairwise -"~italic(F[ST])),
                                barwidth = 20, barheight = 1, title.position = "top")) +
  theme(axis.text.x = element_text(size = 12, face = "bold", angle = 315, vjust = 0.1, hjust=0.02),
        axis.text.y = element_text(size = 12, face = "bold", angle = 315, vjust = 0),
        legend.position = c(0.05,0.95),
        legend.direction = "horizontal",
        legend.justification = c(0,1),
        legend.title = element_text(size = 12, face = "bold"))
pFST
# Export plot to png
png("results/Figs_paper/Main/ALL_FST.png", width = 8.3, height = 10, units = "in", res = 600)
pFST
dev.off()


# SI plots ----------------------------------------------------------
## PCA ALL ----------------------------------------------------------
### Figure S8.14: PCA scree plot ------------------------------------
screep <- PCA_screeplot(pc_sub, nPC = 15)
png("results/Figs_paper/SI/ALL_PCA_BCF_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep
dev.off()
### Figure S8.15: PCA plot ------------------------------------------
# dt_var_expl <- PCA_var_explained(pc_sub)
pc2v3 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 2,
                   yaxis = 3) + PCA_theme
pc3v4 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 3,
                   yaxis = 4) + PCA_theme

pc4v5 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 4,
                   yaxis = 5) + PCA_theme

pc5v6 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 5,
                   yaxis = 6) + PCA_theme

pc6v7 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 6,
                   yaxis = 7) + PCA_theme

pc7v8 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 7,
                   yaxis = 8) + PCA_theme


pc8v9 <- plotPCASI(dtpc, 
                   varls = round(dt_var_expl$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 8,
                   yaxis = 9) + PCA_theme

pc9v10 <- plotPCASI(dtpc, 
                    varls = round(dt_var_expl$variance_explained, 1),
                    dtcolshapes,
                    xaxis = 9,
                    yaxis = 10) + PCA_theme
legend_b <- get_legend(
  pc2v3
)

pca_SI <- plot_grid(pl + theme(legend.position = "none"), 
                    pc2v3 + theme(legend.position = "none"), 
                    pc3v4 + theme(legend.position = "none"),
                    pc4v5 + theme(legend.position = "none"),
                    pc5v6 + theme(legend.position = "none"),
                    pc6v7 + theme(legend.position = "none"),
                    pc7v8 + theme(legend.position = "none"),
                    pc8v9 + theme(legend.position = "none"),
                    pc9v10 + theme(legend.position = "none"),
                    labels = "AUTO", ncol = 3, nrow = 3)

png("results/Figs_paper/SI/ALL_PCA_BCF.png", width = 8.3, height = 10, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(2, .15))
dev.off()

### Figure S8.16: PC loadings ---------------------------------------
dtpcload <- data.frame(pc_sub$loadings)
dtpcload$contig <- as.numeric(gl_sub$other$loc.metrics$CHROM)
dtpcload$position <- as.numeric(gl_sub$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/ALL_PCA_BCF_loadings.png", width = 8.3, height = 11.7, units = "in", res = 600)

op <- par(mfrow = c(10,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

for (i in 1:10){
  plot(dtpcload[,i], ylab = paste("PC", i, sep = ""), xaxt = "n", xlab = "", cex = 0.25)
  axis(side = 2, labels = T)
  box(which = "plot", bty = "l")
}

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5)
dev.off()

## PCA relavent pop -------------------------------------------------
### Figure S8.19: screeplot -----------------------------------------
screep2 <- PCA_screeplot(pc_sub2, nPC = 15)
png("results/Figs_paper/SI/ALL_relevantpop_PCA_BCF_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep2
dev.off()
### Figure S8.20: PCA plot ------------------------------------------
dt_var_expl2 <- PCA_var_explained(pc_sub2)

p2_2v3 <- plotPCASI(dtpc2, 
          varls = round(dt_var_expl2$variance_explained, 1),
          dtcolshapes,
          xaxis = 2,
          yaxis = 3) + PCA_theme

p2_3v4 <- plotPCASI(dtpc2, 
          varls = round(dt_var_expl2$variance_explained, 1),
          dtcolshapes,
          xaxis = 3,
          yaxis = 4) + PCA_theme

p2_4v5 <- plotPCASI(dtpc2, 
          varls = round(dt_var_expl2$variance_explained, 1),
          dtcolshapes,
          xaxis = 4,
          yaxis = 5) + PCA_theme


legend_b2 <- get_legend(
  pl2
)

pca_SI2 <- plot_grid(pl2 + theme(legend.position = "none"), 
                    p2_2v3 + theme(legend.position = "none"), 
                    p2_3v4 + theme(legend.position = "none"),
                    p2_4v5 + theme(legend.position = "none"), 
                    labels = "AUTO", ncol = 2, nrow = 2)
png("results/Figs_paper/SI/ALL_relevantpop_PCA_BCF.png", width = 8.3, height = 9, units = "in", res = 600)
plot_grid(pca_SI2, legend_b2, ncol = 1, rel_heights = c(2, .15))
dev.off()

### Figure S8.21: PC loadings ---------------------------------------
dtpcload2 <- data.frame(pc_sub2$loadings)
dtpcload2$contig <- as.numeric(gl_sub2$other$loc.metrics$CHROM)
dtpcload2$position <- as.numeric(gl_sub2$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload2 <- dtpcload2 %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/ALL_relevantpop_PCA_BCF_loadings.png", width = 8.3, height = 8.4, units = "in", res = 600)

op <- par(mfrow = c(5,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

for (i in 1:5){
  plot(dtpcload2[,i], ylab = paste("PC", i, sep = ""), xaxt = "n", xlab = "", cex = 0.25)
  axis(side = 2, labels = T)
  box(which = "plot", bty = "l")
}

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5)
dev.off()


## Figure S8.22: sNMF on only relevant pops -------------------------
genofile2 <- "data/processed/LEA/03_ALL/ALL_relevantpops_hwe.geno"
# Convert genlight to matrix of 0,1,2 
mat_gl2 <- t(as.matrix(gl_sub2))
# Write table to .geno file which is essentially a matrix 
write.table(mat_gl2, file = genofile2,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)

Krange <- 1:20
snmf_ALLsub <- snmf(genofile2, K = Krange, repetitions = 10, ploidy = 2,
                    entropy = T, alpha = 100, project = "new")
# snmf_ALLsub <- load.snmfProject(file = "data/processed/LEA/03_ALL/ALL_relevantpops_hwe.snmfProject")

plot(snmf_ALLsub, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), ALL dataset, after HWE filter, only relevant pops.")
# Calculate mean entropy
crossentropy.matALLrel <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_ALLsub, K = x)})))
mean.entropyALLrel   <- apply(crossentropy.matALLrel, MARGIN=1, FUN=mean, na.rm=TRUE)
if(any(diff(mean.entropyALLrel)>0)){
  bestK <- unname(which(diff(mean.entropyALLrel)>0)[1])
} else {
  bestK <- unname(Krange[1])
}

q.ALLrel.ls <- merge_sNMF_qmatrix_multiK(snmf_ALLsub, gl_sub2, k2plot = 1:10)
pop_label_ALLrel<- gl_sub2$other$ind.metrics[, c("ID","popdef2")]
poporder <- c("IND: Maharashtra subpopulation A", "AUS: Melbourne", "Fiji", 
              "NZ: Napier", "NZ: Other", "AUS: Sydney")
plotQ(q.ALLrel.ls[2:10],imgoutput="join",
      showdiv = T,
      showindlab=F,
      grplab=pop_label_ALLrel[, c("popdef2"), drop = F],
      subsetgrp=poporder,
      grplabangle = 90,
      grplabheight = 5,
      panelratio = c(4,9),
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
      outputfilename="results/Figs_paper/SI/ALL_relevantpops_sNMF_BCF_k2-5",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on ALL dataset (subset relevant pops) (after HWE filter), K = 2-10, 10 repetitions each"
)



## Figures S8.18: FST ALL subpopulations ----------------------------
# This is the FST population pairwise FST matrix between all populations
# as defined by popdef1
pALL_FSTmatALL <- plot_FST_mat_from_xlsx(FSTxlsx = "results/FST/03_ALL/pairwiseFST_ALL_BCF_hwe_100bs.xlsx", 
                                         poporder = c("South Africa", "Hawaii", "Gold Coast",
                                                      "Sydney", "Sydney (ROM)", "Helensville", "Ngunguru",
                                                      "Waitakeres", "Leigh", "Auckland (ROM)", "Waiheke",
                                                      "Hamilton (ROM)", "Great Barrier Island", "Thames", 
                                                      "Melbourne", "Melbourne (ROM)", "Napier", "Napier (ROM)", "Fiji",
                                                      "Maharashtra subpop. A", "Maharashtra", "Madhya Pradesh", "Karnataka",
                                                      "Andhra Pradesh", "Tamil Nadu", 
                                                      "Uttar Pradesh", "Odisha", "West Bengal"),
                                         dt1spreadsheetname = "FST_100bs_modROMsep",
                                         dt2spreadsheetname = "FST_100bs_val_modROMsep",
                                         dt3spreadsheetname = "pval_100bs_modROMsep",
                                         fontsize = 1.5)

pALL_FSTmatALL <- pALL_FSTmatALL + 
  scale_fill_distiller(palette = "Spectral",
                       direction = 1,
                       limits = c(-0.001, 0.32)) + 
  scale_y_discrete(labels = c("South Africa", "Hawaii", "Gold Coast",
                              "Sydney", "Sydney (ROM)", "Helensville", "Ngunguru",
                              "Waitakeres", "Leigh", "Auckland (ROM)", "Waiheke",
                              "Hamilton (ROM)", "Great Barrier Island", "Thames",
                              "Melbourne", "Melbourne (ROM)", "Napier", "Napier (ROM)", "Fiji",
                              "Maharashtra subpop. A", "Maharashtra", "Madhya Pradesh", "Karnataka",
                              "Andhra Pradesh", "Tamil Nadu",
                              "Uttar Pradesh", "Odisha", "West Bengal"),
                              position = "right") +
  scale_x_discrete(labels = c("South Africa", "Hawaii", "Gold Coast",
                              "Sydney", "Sydney (ROM)", "Helensville", "Ngunguru",
                              "Waitakeres", "Leigh", "Auckland (ROM)", "Waiheke",
                              "Hamilton (ROM)", "Great Barrier Island", "Thames",
                              "Melbourne", "Melbourne (ROM)", "Napier", "Napier (ROM)", "Fiji",
                              "Maharashtra subpop. A", "Maharashtra", "Madhya Pradesh", "Karnataka",
                              "Andhra Pradesh", "Tamil Nadu",
                              "Uttar Pradesh", "Odisha", "West Bengal"))

# Export plot to png
png("results/Figs_paper/SI/ALL_FST_modROMsep.png", width = 8.3, height = 10, units = "in", res = 600)
pALL_FSTmatALL
dev.off()

