# All_popstruc_SI_DArT.R --------------------------------------------
# This script is to show that results from the DArT dataset
# is comparable to BCFtools datasets
# 
# This script does the following:
# 1. Remove AUS: Other and AUS: Other Victoria + Townsville + Cairns
#    (only keeping AUS: Melbourne, Sydney, and Gold Coast)
# 2. Remove samples with outlier observed heterozygousity 
#    These are
# 2. Remove loci with call rates < 0.8, reproducibility < 0.95 and MAC = 1 and loci where SNPs only present in one individual
# 3. Export loci and individual list for filtering on VCF, for
#    plink LD clumping and also for other procedures.
#    The reason to do the LD clumping in PLINK is so that it is consistent
#    with how the data from the other pipeline is processed. 
# 4. Report call rates of individuals and loci. Report loci reproducibility, 
#    output to PDF
#
# Define input ------------------------------------------------------
dartcsv <- "../01_download_data/DArT/DArT_ALL/Report_DImy20-5737_2_moreOrders_SNP_1.csv"
metacsv <- '../01_download_data/DArT/metadata_all_samples.csv'
# Define location to store .geno file for sNMF analysis -------------
genofile <- "data/processed/LEA/03_ALL/ALL_DART_thinnedonly.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/PCA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
# Define output folder to store subsampled dataset ------------------
dir.create("data/processed/sample_subset", recursive = T, showWarnings = F)

# Load libraries ----------------------------------------------------
library(tidyverse)
library(dartR)
library(cowplot)
library(LEA) # Version 3.2.0. Some function names change in newer versions.
library(pophelper)
source("../shared_scripts/functions.R")
# Read in data ------------------------------------------------------
gl_all <- gl.read.dart(filename = dartcsv,
                       ind.metafile = metacsv)
pop(gl_all) <- gl_all$other$ind.metrics$popdef1
# save(gl_all, file = "data/gl_all.Rdata")
# load("data/gl_all.Rdata")
# Remove replicats --------------------------------------------------
# Kyle replicates
ind_rep <- grep("^[0-9]*[a-zA-Z]$", indNames(gl_all), value = T)
# Append our in-house replicates
ind_rep <- c(ind_rep, "M0092", "M0093", "M0094", "M0186", "M0187", "M0370", "M0188", "M0280", "M0281", "M0282", "M0368", "M0369", "M0237")
# Drop replicates 
gl_all_sub1 <- gl.drop.ind(gl_all, ind.list = ind_rep, recalc = T)

# Remove admixed samples --------------------------------------------
# The admixed samples from Australia had no population assignment.
gl_all_sub1a <- gl.drop.pop(gl_all_sub1, pop.list = c(NA), 
                            mono.rm = T, recalc = T)
# Report heterozygousity after removal of admixed samples -----------
png("results/Figs_paper/SI/DArT_Ho_indiv.png", height = 5, width = 8, units = "in", res = 600)
dt_all_sub1 <- gl.report.heterozygosity(gl_all_sub1a, method = "ind")
dev.off()

gl_all_sub1b <- gl.filter.heterozygosity(gl_all_sub1a, t.upper = 0.175)

# Filter for loci that has only been called once on the reference ---
loc_alncnt_not1 <- locNames(gl_all_sub1b)[gl_all_sub1b$other$loc.metrics$AlnCnt_Indian_Myna_v2.1 != 1]
gl_all_sub2 <- gl.drop.loc(gl_all_sub1b, loc.list = loc_alncnt_not1)

# Filter loci based on call rates and reproducibility ---------------
# This also removes single and doubletons. 
gl_all_sub3 <- gl_filt_CR_rep_MAC(gl_all_sub2, CR = 0.8, 
                                  rep = 0.95, MAC_threshold = 1, 
                                  SNP_in_single_samp_rm = T)

gl_all_sub4 <- gl.prune.SNP.dist(gl_all_sub3, bp_threshold = 100000)

# Re-perform PCA on the ALL dataset, nmax = 20 ----------------------
dtsub <- read.csv("../04_BCFtools_data_analysis/data/processed/ALL_BCFtools_subset_nmax20.csv")
# Remove a/b at the end of name which only exists in STACKS/BCFtools 
# dataset as there is data from DART replicates
idsub <- gsub(pattern = "[a-z]$", "", dtsub$subsampled_ID)
gl_sub <- gl.keep.ind(gl_all_sub4, ind.list = idsub, recalc = T, mono.rm = T)
pc_sub <- gl.pcoa(gl_sub, nfactors=10, parallel = T, n.cores = 8)
# Save pc object to .Rdata for later access
# save(pc_sub, file = "results/PCA/03_ALL/PCA_ALL_DART_thinnedonly_nmax20.Rdata")
# load("results/PCA/03_ALL/PCA_ALL_DART_thinnedonly_nmax20.Rdata")
# Quick plot to see variance explained
gl.pcoa.plot(pc_sub, gl_sub, pop.labels= "pop") 
# Create dataframe of PC scores
dtpc <- data.frame(pc_sub$scores)
dtpc$pop <- pop(gl_sub)
# write.csv(dtpc, file = "results/PCA/03_ALL/PCA_ALL_DART_thinnedonly_nmax20.csv")
# dtpc <- read.csv("results/PCA/03_ALL/PCA_ALL_DART_thinnedonly_nmax20.csv")
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
# SI ----------------------------------------------------------------
## PCA ALL ----------------------------------------------------------
### PCA scree plot --------------------------------------------------
screep <- PCA_screeplot(pc_sub, nPC = 15)
png("results/Figs_paper/SI/ALL_PCA_DART_thinnedonly_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep
dev.off()
### PCA plot -------------------------------------------------------------
dt_var_expl <- PCA_var_explained(pc_sub)
round(dt_var_expl$variance_explained, 1)

pl <- plotPCASI(dtpc, 
                varls = round(dt_var_expl$variance_explained, 1),
                dtcolshapes,
                xaxis = 1,
                yaxis = 2) + PCA_theme


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

png("results/Figs_paper/SI/ALL_PCA_DART_thinnedonly.png", width = 8.3, height = 10, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(2, .15))
dev.off()

### PC loadings -----------------------------------------------------
dtpcload <- data.frame(pc_sub$loadings)
dtpcload$contig <- as.numeric(gl_sub$other$loc.metrics$Chrom_Indian_Myna_v2.1)
dtpcload$position <- as.numeric(gl_sub$other$loc.metrics$ChromPos_Indian_Myna_v2.1)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/ALL_PCA_DART_thinnedonly_loadings.png", width = 8.3, height = 11.7, units = "in", res = 600)

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
# sNMF --------------------------------------------------------------
# Previous iterations of the analysis found that optimal K 
# is likely ~ 7 or 8
# genofile <- "data/processed/LEA/03_ALL/ALL_DART_thinnedonly.geno"
# Convert genlight to matrix of 0,1,2 -------------------------------
mat_gl <- t(as.matrix(gl_all_sub4))
# Write table to .geno file which is essentially a matrix -----------
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
Krange <- 1:20
snmf_ALL <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2,
                 entropy = T, alpha = 100, project = "new")
snmf_ALL <- load.snmfProject("data/processed/LEA/03_ALL/ALL_DART_thinnedonly.snmfProject")

plot(snmf_ALL, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), ALL")

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
png("results/Figs_paper/SI/ALL_sNMF_meanxentropy_DART_thinnedonly_k1-20.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropy, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()

rep <- 1:10
# Create labels for samples 
# Note that the gl object here is the same object used to run sNMF that was
# called earlier.
samplenames <- indNames(gl_all_sub4)
numind <- length(samplenames)
# Make table of population information
gl_all_sub4$other$ind.metrics$popdef1 <- as.character(gl_all_sub4$other$ind.metrics$popdef1) # convert to character from factor
pop_label_sub <- gl_all_sub4$other$ind.metrics[, c("id","popdef1")]
# Attach other information that might be more useful
pop_label_sub <- pop_label_sub %>% 
  mutate(snmfcluster = replace(popdef1, popdef1 %in% c("AUS: Gold Coast"), "Group 1")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("AUS: Sydney"), "Group 2")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("NZ: Other"), "Group 3")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("IND: Maharashtra subpopulation A", "Fiji", "NZ: Napier", "AUS: Melbourne"), "Group 4")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("IND: Other"), "Group 5")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("Hawaii"), "Group 6")) %>%
  mutate(snmfcluster = replace(snmfcluster, popdef1 %in% c("South Africa"), "Group 7")) %>%
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
pop_trans <- data.frame(popdef1 = pop_labels,
                        pop_codes = pop_abb)
pop_label_sub <- merge(pop_label_sub, pop_trans)
# Make sure that the ID is sorted properly
pop_label_sub <- pop_label_sub[match(indNames(gl_all_sub4), pop_label_sub$id),]

q.ALL.ls <- merge_sNMF_qmatrix_multiK(snmf_ALL, gl_all_sub4, k2plot = 1:10)


plotQ(q.ALL.ls[2:10],imgoutput="join",
      # clustercol = c("#cab2d6", "#a6cee3",
      #                "#6a3d9a", "#b2df8a",
      #                "yellow", "#fb9a99",
      #                "#ff7f00", "#e31a1c"),
      showdiv = T,
      showindlab=F,
      grplab=pop_label_sub[, c("pop_codes", "snmfcluster", "pcacluster")],
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
      outputfilename="results/Figs_paper/SI/ALL_DART_thinnedonly_sNMF_BCF_k2-10",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on STACKS ALL dataset (nmax = 20), K = 2-10, 10 repetitions each"
)
