# NZ_popstruc_SI_DArT.R ---------------------------------------------
# This script is to show that results from the DArT dataset
# is comparable to BCFtools datasets. Appendix S10.4
# 
# Define input ------------------------------------------------------
dartcsv <- "../01_download_data/DArT/DArT_NZ/Report_DImy20-5737_SNP_2.csv"
metacsv <- '../01_download_data/DArT/metadata_all_samples.csv'
# Define location to store .geno file for sNMF analysis -------------
genofile <- "data/processed/LEA/01_NZ/NZ_DART_thinnedonly.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/01_NZ/", recursive = T, showWarnings = F)
dir.create("results/PCA/01_NZ/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)

# Load libraries ----------------------------------------------------
library(tidyverse)
library(dartR)
library(cowplot)
library(LEA) # Version 3.2.0. Some function names change in newer versions.
library(pophelper)
source("../shared_scripts/functions.R")
# Read in data ------------------------------------------------------
gl_NZ <- gl.read.dart(filename = dartcsv,
ind.metafile = metacsv)
# save(gl_NZ, file = "data/gl_NZ.Rdata")
# load("data/gl_NZ.Rdata")
# Remove replicates -------------------------------------------------
# Kyle replicates
ind_rep <- grep("^[0-9]*[a-zA-Z]$", indNames(gl_NZ), value = T)
# Append our in-house replicates
ind_rep <- c(ind_rep, "M0092", "M0093", "M0094", "M0186", "M0187", "M0370", "M0188", "M0280", "M0281", "M0282", "M0368", "M0369", "M0237")
# Drop replicates 
gl_NZ_sub1 <- gl.drop.ind(gl_NZ, ind.list = ind_rep, recalc = T)

# Report heterozygousity after removal of replicates ----------------
#
# This is not presented in paper. Used in ALL dataset to identify
# samples with outlier heterozygosity levels.
#
png("results/NZ_DArT_Ho_indiv.png", height = 5, width = 8, units = "in", res = 600)
dt_NZ_sub1 <- gl.report.heterozygosity(gl_NZ_sub1, method = "ind")
dev.off()

# Filter for loci that has only been called once on the reference ---
loc_alncnt_not1 <- locNames(gl_NZ_sub1)[gl_NZ_sub1$other$loc.metrics$AlnCnt_Indian_Myna_v2.1 != 1]
gl_NZ_sub2 <- gl.drop.loc(gl_NZ_sub1, loc.list = loc_alncnt_not1)

# Filter loci based on call rates and reproducibility ---------------
# This also removes single and doubletons. 
gl_NZ_sub3 <- gl_filt_CR_rep_MAC(gl = gl_NZ_sub2, CR = 0.8, 
                                  rep = 0.95, MAC_threshold = 1, 
                                  SNP_in_single_samp_rm = T)

# Prune SNPs every 100kb --------------------------------------------
gl_NZ_sub4 <- gl.prune.SNP.dist(gl_NZ_sub3, bp_threshold = 100000)

# Count samples per location ----------------------------------------
gl_NZ_sub4@other$ind.metrics %>% dplyr::count(popdef2)
pop(gl_NZ_sub4) <- gl_NZ_sub4$other$ind.metrics$popdef2
# Re-perform PCA on the full dataset --------------------------------
pc <- gl.pcoa(gl_NZ_sub4, nfactors=10, parallel = T, n.cores = 8)
# Quick plot to see variance explained ------------------------------
gl.pcoa.plot(pc, gl_NZ_sub4, pop.labels= "pop") 
# Create dataframe for pricipal component axes ----------------------
dtpc <- data.frame(pc$scores)
save(pc, file = "results/PCA/01_NZ/PCA_NZ_DART_thinnedonly.Rdata")
load("results/PCA/01_NZ/PCA_NZ_DART_thinnedonly.Rdata")
# Assign population information -------------------------------------
dtpc$popdef2 <- gl_NZ_sub4$other$ind.metrics$popdef2
# Assign whether it is ROM or Modern
dtpc <- dtpc %>% 
  mutate(ROM = ifelse(grepl("(ROM)", popdef2), 
                      "ROM", 
                      "Modern")) %>%
  mutate(ROM = replace(ROM, "(ROM)" %in% popdef2, "ROM samples")) %>%
  # Assign location regardless of whether it is ROM or not
  # Note that Napier (ROM) is changed to "Napier " with a space at the end
  # This is for plotting purpuses
  mutate(location = replace(as.character(popdef2), popdef2 == "Kaikohe (ROM)", "Kaikohe")) %>%
  mutate(location = replace(location, location == "Auckland (ROM)", "Auckland")) %>%
  mutate(location = replace(location, location == "Hamilton (ROM)", "Hamilton")) %>%
  mutate(location = replace(location, location == "Taupo (ROM)", "Taupo")) %>%
  mutate(location = replace(location, location == "Napier (ROM)", "Napier ")) 

write.csv(dtpc, file = "results/PCA/01_NZ/PCA_NZ_DART_thinnedonly.csv")
# Get lat and long for each location to sort legends by
dtlatlon <- gl_NZ_sub4@other$ind.metrics %>% 
  dplyr::group_by(popdef2) %>%
  dplyr::summarise(lat = mean(latitude),
            lon = mean(longitude))
dtlatlon$popdef2[order(dtlatlon$lon)]
lab_order <- c("Helensville", "Ngunguru", "Waitakeres", 
               "Onehunga", "Leigh", "Waiheke",
               "Great Barrier Island", "Thames", "Napier", 
               "Kaikohe","Auckland", 
               "Hamilton", "Taupo", "Napier ")
lab_order_abb <- c("Helensville (HEL)", "Ngunguru (NGU)", "Waitakeres (WTK)", 
                   "Onehunga (ONE)", "Leigh (LEI)", "Waiheke (WHK)",
                   "Great Barrier Island (GBI)", "Thames (THA)", "Napier (NAP)", 
                   "Kaikohe (KKH)","Auckland (AKL)", 
                   "Hamilton (HAM)", "Taupo (TAU)", "Napier (NAR)")
shape_order <- c(1, 13, 3, 12, 10, 8, 14, 11, 24,15,16,18,19, 17)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_order <- c(gg_color_hue(13), gg_color_hue(13)[9])
dtcolshapes <- data.frame(pop = lab_order,
                          pop_lab = lab_order_abb,
                          col = col_order,
                          shape = shape_order)
dtpc <- read.csv("results/PCA/01_NZ/PCA_NZ_DART_thinnedonly.csv")
dtpc$location <- factor(dtpc$location,
                        levels = lab_order)

dtvar <- PCA_var_explained(pc)
round(dtvar$variance_explained, 1)
# Create PC1 vs PC2 plot --------------------------------------------
pl <- plotPCASI(dtpc, 
                varls = round(dtvar$variance_explained, 1),
                dtcolshapes,
                xaxis = 1,
                yaxis = 2, popcol = "location") + 
  PCA_theme

# grid.text("ROM", x = unit(0.48, "npc"), y = unit(0.97, "npc"), gp = gpar(fontsize = 8, fontface = "bold"))
# grid.text("Modern", x = unit(0.36, "npc"), y = unit(0.97, "npc"), gp = gpar(fontsize = 8, fontface = "bold"))

# Create PCA related SI plots ---------------------------------------
##  Screeplot 
screep <- PCA_screeplot(pc, nPC = 15)
png("results/Figs_paper/SI/NZ_PCA_DART_thinnedonly_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep
dev.off()
## PC1 vs PC2
## PC2 vs PC3
pc2v3 <- plotPCASI(dtpc, 
                   varls = round(dtvar$variance_explained, 1),
                   dtcolshapes,
                   xaxis = 2,
                   yaxis = 3, popcol = "location") + 
  PCA_theme

legend_b <- get_legend(
  pc2v3
)

pca_SI <- plot_grid(pl + theme(legend.position = "none",
                               aspect.ratio = 1), 
                    pc2v3 + theme(legend.position = "none",
                                  aspect.ratio = 1), 
                    labels = "AUTO", ncol = 2)
png("results/Figs_paper/SI/NZ_PCA_DART_thinnedonly.png", width = 8.3, height = 4.2, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(1, .15))
dev.off()

## Plot PC loadings 
dtpcload <- data.frame(pc$loadings)
gl_NZ_sub4$other$loc.metrics$CHROM <- gl_NZ_sub4$other$loc.metrics$Chrom_Indian_Myna_v2.1
gl_NZ_sub4$other$loc.metrics$POS <- gl_NZ_sub4$other$loc.metrics$ChromPos_Indian_Myna_v2.1
dtpcload$contig <- as.numeric(gl_NZ_sub4$other$loc.metrics$CHROM)
dtpcload$position <- as.numeric(gl_NZ_sub4$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/NZ_PCA_DART_thinnedonly_loadings.png", width = 8.3, height = 5, units = "in", res = 600)

op <- par(mfrow = c(3,1),
          oma = c(2,2,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)

plot(dtpcload[,1], ylab = "PC1", xaxt = "n", xlab = "", cex = 0.25)
# axis(side = 1,
#      at=1:nlevels(x),
#      labels = F)
axis(side = 2, labels = T)
box(which = "plot", bty = "l")

plot(dtpcload[,2], ylab = "PC2", xaxt = "n", xlab = "", cex = 0.25)
axis(side = 2, labels = T)
box(which = "plot", bty = "l")

plot(dtpcload[,3], ylab = "PC3", xaxt = "n", xlab = "", cex = 0.25)
axis(side = 2, labels = T)
box(which = "plot", bty = "l")

title(xlab = "SNP ordered by CONTIG number and position",
      ylab = "Loadings",
      outer = TRUE, line = 0.5,
      cex.lab = 1.5
)
dev.off()

# sNMF --------------------------------------------------------------
## Convert genlight to matrix of 0,1,2 ------------------------------
mat_gl <- t(as.matrix(gl_NZ_sub4))
## Write table to .geno file which is essentially a matrix ----------
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
## Perform sNMF -----------------------------------------------------
Krange <- 1:5
snmf_NZ_thinnedonly <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2, 
                            entropy = T, alpha = 100, project = "new")
snmf_NZ_thinnedonly <- load.snmfproject("data/processed/LEA/01_NZ/NZ_DART_thinnedonly.snmfProject")
plot(snmf_NZ_thinnedonly, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), NZ, thinned dataset only")

# Calculate mean entrop
crossentropy.mat <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_NZ_thinnedonly, K = x)})))
mean.entropy   <- apply(crossentropy.mat, MARGIN=1, FUN=mean, na.rm=TRUE)
# SI: Export mean cross-entropy plot---------------------------------
png("results/Figs_paper/SI/NZ_sNMF_meanxentropy_DART_thinnedonly_k1-5.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropy, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()

if(any(diff(mean.entropy)>0)){
  bestK <- unname(which(diff(mean.entropy)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
bestK
# plot k = 3
rep <- 1:10
samplenames <- indNames(gl_NZ_sub4)
numind <- length(samplenames)
pop_label_sub <- gl_NZ_sub4$other$ind.metrics[, c("id","popdef2", "popdef1")]
pop_labels <- c("Auckland (ROM)", "Great Barrier Island",
                "Hamilton (ROM)", "Helensville", 
                "Kaikohe (ROM)", "Leigh",
                "Ngunguru", "Auckland",
                "Taupo (ROM)", "Thames",
                "Waiheke", "Waitakeres",
                "Napier", "Napier (ROM)")
dtlatlon$popdef2[order(dtlatlon$lon)]
pop_trans <- data.frame(popdef2 = pop_labels,
                        loc_time = LETTERS[1:length(pop_labels)],
                        loc_time_abb = c("AKR", "GBI", "HAM", "HEL", "KKH", "LEI",
                                         "NGU", "AKL", 
                                         "TAU", "THA", "WHK", "WTK", "NAP", "NAR"))
pop_label_sub <- merge(pop_label_sub, pop_trans)
colnames(pop_label_sub) <- c("popdef2", "id", "population", "location/sample", "loc_time_abb")
pop_label_sub <- pop_label_sub[match(indNames(gl_NZ_sub4), pop_label_sub$id),]
pop_label_sub$population <- as.character(pop_label_sub$population)
# SI plot for K = 1-5
q.df.ls <- merge_sNMF_qmatrix_multiK(snmf_NZ_thinnedonly, gl_NZ_sub4, k2plot = 1:5)
plotQ(q.df.ls[2:5],imgoutput="join",
      showdiv = T,
      showindlab=F,
      grplab=pop_label_sub[, c("loc_time_abb", "population"), drop = F],
      # subsetgrp=c("NZ: Other", "NZ: Napier"),
      subsetgrp=c("KKH", "HEL", "NGU", "WTK", "ONE", "LEI", 
                  "AKL", "WHK", "HAM", "GBI", "THA", "TAU", "NAP", "NAR"),
      grplabangle = 90,
      # grplabsize = 2.5,
      # splabsize = 10,
      selgrp="loc_time_abb",
      splab = c("K=2", "K=3", "K=4", "K=5"),
      ordergrp=T,
      showlegend=T,
      showtitle=T,showsubtitle=T,
      height=2.3,
      width = 16.25,
      indlabspacer=-1,
      barbordercolour=NA,barbordersize=0,
      legendrow = 2,
      legendkeysize = 14,
      legendtextsize = 12,
      legendspacing = 10,
      outputfilename="results/Figs_paper/SI/NZ_sNMF_DART_thinnedonly_k2-5",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on NZ dataset, K = 2-5, 10 repetitions each"
)

