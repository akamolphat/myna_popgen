# This script makes plots for the figure for publication for the IND dataset
#
# This script is for the BCFtools dataset
#
# Define input files ------------------------------------------------
#
vcffile <- "../01_download_data/BCFtools_IND/variant_calls.IND.bialminQ30minGQ30DP15-100.norep.highnegfis.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2.csv"
# Define .geno file to store data for performing LEA ----------------
genofile <- "data/processed/LEA/02_IND/IND_hwe.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/02_IND/", recursive = T, showWarnings = F)
dir.create("results/PCA/02_IND/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/Main/", recursive = T, showWarnings = F)
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
# Define Pop4
gl$other$ind.metrics$IND_pop_recode <- gl$other$ind.metrics$popdef2
gl$other$ind.metrics <- gl$other$ind.metrics %>%
  mutate(IND_pop_recode = replace(IND_pop_recode, ID %in% paste("M0", seq(210, 215), sep = ""), "Maharashtra subpop. A"))
# pop(gl) <- gl$other$ind.metrics$IND_pop_recode
pop(gl) <- gl@other$ind.metrics$popdef2
# Count samples per location ----------------------------------------
gl@other$ind.metrics %>% dplyr::count(popdef2)
# Re-perform PCA on the full dataset --------------------------------
pc <- gl.pcoa(gl, nfactors=10, parallel = T, n.cores = 8)
# Save pc object for reproducibility
save(pc, file = "results/PCA/02_IND/PCA_IND_hwe.Rdata")
load("results/PCA/02_IND/PCA_IND_hwe.Rdata")
# Quick plot to see variance explained ------------------------------
gl.pcoa.plot(pc, gl, pop.labels= "pop") 
# Create dataframe for pricipal component axes ----------------------
dtpc <- data.frame(pc$scores)
# Assign population information -------------------------------------
dtpc$popdef2 <- gl$other$ind.metrics$popdef2
dtpc$IND_pop_recode <- gl$other$ind.metrics$IND_pop_recode
write.csv(dtpc, file = "results/PCA/02_IND/PCA_IND_hwe.csv")
dtpc <- read.csv("results/PCA/02_IND/PCA_IND_hwe.csv")
# Get lat and long for each location to sort legends by
dtlatlon <- gl@other$ind.metrics %>% 
  group_by(popdef2) %>%
  summarise(lat = mean(latitude),
            lon = mean(longitude))
# dtlatlon$popdef2[order(dtlatlon$lon)]
lab_order <- dtlatlon$popdef2[order(dtlatlon$lon)]
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
# Convert dtpc population column to factor and in order ------------
dtpc$popdef2 <- factor(dtpc$popdef2,
                       levels = lab_order)
# Create PC1 vs PC2 plot --------------------------------------------
pl <- plotPCASI(dtpc, 
                varls = round(dtvar$variance_explained, 1),
                dtcolshapes,
                xaxis = 1,
                yaxis = 2, popcol = "popdef2") + 
  PCA_theme
pl <- pl +
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

# Create Isolation by distance plot -------------------------------------------
# Attach latlong data to gl$other$latlong
latlon_mat <- as.matrix(gl$other$ind.metrics[,c("latitude", "longitude")])
colnames(latlon_mat) <- c("lat", "lon")
rownames(latlon_mat) <- indNames(gl)
gl$other$latlon <- latlon_mat

# Performs IBD manually instead of the default ----------------------
# This is because the genetic distance matrix can be changed
# The default does not use bootstrapping to calculate FST
#
# Generally the results are the same, but this done for completeness
# purpose for the final plot
#
# Retain populations with sufficient n (n > 5)
#
#
gl$other$ind.metrics %>% dplyr::count(popdef2)
# Drop populations and samples 
glsub <- gl.drop.pop(gl, pop.list = c("Gujarat"))
# Drop samples from Maharashtra subpop. A
glsub <- gl.drop.ind(glsub, ind.list = paste("M0", seq(210,215), sep = ""))
# Covert coordinates to Mercator projection
coords <- dismo::Mercator(glsub@other$latlon[, c("lon", "lat")])
# Create distances between each population, using the mean 
# coordinate of each population. All samples from the same
# population have the same coordinates
pop.xy <- apply(coords, 2, function(a)
  tapply(a, pop(glsub), mean, na.rm = T))
## Create distance matrix in metres ----------------------------------
Dgeo <- as.dist(dist(pop.xy))
## Create FST matrix ------------------------------------------------
## Standard
fst <- gl.fst.pop(glsub, nboots = 100)
Dgen <- as.dist(fst$Fsts)
# Create FST matrix using the Reich et al. (2009) method
fstreich <- reich.fst(glsub, bootstrap = 100)
Dgen2 <- as.dist(fstreich$fsts)
# Has to reorder the matrix to make sure that the 
# column and row names are the same
oo <- order(colnames(as.matrix(Dgen)))
Dgen <- as.dist(as.matrix(Dgen)[oo, oo])
oo <- order(colnames(as.matrix(Dgeo)))
Dgeo <- as.dist(as.matrix(Dgeo)[oo, oo])
oo <- order(colnames(as.matrix(Dgen2)))
Dgen2 <- as.dist(as.matrix(Dgen2)[oo, oo])

Dgeo <- as.dist(Dgeo)
Dgen <- as.dist(Dgen)
Dgen2 <- as.dist(Dgen2)

# save(Dgen, Dgen2, Dgeo, file = "results/Figs_paper/IND_IBD_distance_matrix.Rdata")
# load("results/Figs_paper/IND_IBD_distance_matrix.Rdata")

# Checking if the results from doing it manually is similar to
# the default gl.ibd function
gl.ibd(Dgen = Dgen, Dgeo = Dgeo, Dgeo_trans = "log((Dgeo + 1)/1000)")
gl.ibd(Dgen = Dgen2, Dgeo = Dgeo)
gl.ibd(glsub)

## Test out different transformations -------------------------------
# The IBD is significant in all transformations
gl.ibd(Dgen = Dgen2, Dgeo = Dgeo, Dgen_trans='Dgen/(1-Dgen)')
gl.ibd(Dgen = Dgen, Dgeo = Dgeo, Dgen_trans='Dgen/(1-Dgen)')
gl.ibd(Dgen = Dgen, Dgeo = Dgeo, Dgen_trans='Dgen/(1-Dgen)', Dgeo_trans = "Dgeo")
gl.ibd(Dgen = Dgen, Dgeo = Dgeo)

## IBD using FST/(1-FST) vs log(distance) ---------------------------
# Rousset (1997) (see gl.ibd manual) suggests FST/(1-FST) vs log(distance)
# for IBD analysis
#
# Use FST/(1-FST) vs Log(Distance in km)
# Run analysis with the distance matrices. Note that the plot is 
# not saved and therefore needs to be replotted
p <- gl.ibd(Dgen = Dgen, Dgeo = Dgeo, Dgen_trans='Dgen/(1-Dgen)', Dgeo_trans = "log(Dgeo/1000)")
manteltest <- p$mantel 
lm_eqn <- function(df,
                   r = manteltest$statistic,
                   pp = manteltest$signif) {
  # Function from the source code of gl.ibd
  m <- lm(Dgen ~ Dgeo, df)
  eq <-
    substitute(
      italic(y) == a + b %.% italic(x) * "," ~ ~ italic(R) ^ 2 ~ "=" ~ r2 * "," ~ ~
        italic(p) ~ "=" ~ pp,
      list(
        a = format(unname(coef(m)[1]),
                   digits = 2),
        b = format(unname(coef(m)[2]), digits = 2),
        r2 = format(summary(m)$r.squared, digits = 3),
        pp = format(pp, digits = 3)
      )
    )
  as.character(as.expression(eq))
}
res <- data.frame(Dgen = as.numeric(p$Dgen), Dgeo = as.numeric(p$Dgeo))
p3 <- ggplot(res, aes(x = Dgeo, y = Dgen)) + geom_point() + 
  geom_smooth(method = "lm", se = TRUE, colour = "black") + 
  ylab("Dgen_trans") + xlab("Dgeo_trans") +
  annotate(
    "text",
    label = lm_eqn(res),
    x = Inf,
    y = -Inf,
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.5,
    size = 2.5
  ) + 
  theme_classic()
p3 <- p3 + ylab(expression(italic(F[ST])/(1-italic(F[ST])))) + 
  xlab("ln(distance in km)") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8, face = "bold"))
# Make plot with PCA and mantel test --------------------------------
prow <- plot_grid(pl, p3, ncol = 2, rel_widths = c(4,4), labels = "AUTO")

prow
png("results/Figs_paper/Main/IND_pop_str_PCA_mantel.png", width = 8.3, height = 4, units = "in", res = 600)
prow
dev.off()
# Appendix S8.2, SI plots -------------------------------------------
## PCA --------------------------------------------------------------
### Figure S8.8: screeplot ------------------------------------------
screep <- PCA_screeplot(pc, nPC = 15)
png("results/Figs_paper/SI/IND_PCA_BCF_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep
dev.off()
### Figure S8.9: PC1 vs 2-5 -----------------------------------------
# dtvar <- PCA_var_explained(pc)  # Already performed earlier
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
  pc2v3 + guides(colour=guide_legend(nrow=2,byrow=T), 
               shape=guide_legend(nrow=2,byrow=T))
)

pca_SI <- plot_grid(pl + theme(legend.position = "none"), 
                    pc2v3 + theme(legend.position = "none"), 
                    pc3v4 + theme(legend.position = "none"),
                    pc4v5 + theme(legend.position = "none"), labels = "AUTO", ncol = 2, nrow = 2)
png("results/Figs_paper/SI/IND_PCA_BCF.png", width = 8.3, height = 8.5, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(2, .15))
dev.off()

### Figure S8.10: PC loadings ---------------------------------------
dtpcload <- data.frame(pc$loadings)
dtpcload$contig <- as.numeric(gl$other$loc.metrics$CHROM)
dtpcload$position <- as.numeric(gl$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

# Base plot method seems nice
png("results/Figs_paper/SI/IND_PCA_BCF_loadings.png", width = 8.3, height = 8.4, units = "in", res = 600)

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

## FST matrix -------------------------------------------------------
poporder1 <- c("Maharashtra","Madhya Pradesh","Karnataka",
              "Andhra Pradesh","Tamil Nadu",
              "Uttar Pradesh","Odisha","West Bengal")

pIND_FSTmatmerge <- plot_FST_mat_from_xlsx(FSTxlsx = "results/FST/02_IND/pairwiseFST_IND_BCF_hwe_100bs.xlsx",
                                                    poporder = poporder1,
                                          dt1spreadsheetname = "FST_100bs_popdef2",
                                          dt2spreadsheetname = "FST_100bs_val_popdef2",
                                          dt3spreadsheetname = "pval_100bs")


poporder2 <- c("Maharashtra subpop. A", "Maharashtra","Madhya Pradesh","Karnataka",
               "Andhra Pradesh","Tamil Nadu",
               "Uttar Pradesh","Odisha","West Bengal")

pIND_FSTmatALL <- plot_FST_mat_from_xlsx(FSTxlsx = "results/FST/02_IND/pairwiseFST_IND_BCF_hwe_100bs.xlsx",
                                         poporder = poporder2)
poporder2lab <- c("Maharashtra subpop. A", "Maharashtra","Madhya Pradesh","Karnataka",
                  "Andhra Pradesh","Tamil Nadu",
                  "Uttar Pradesh","Odisha","West Bengal")
pIND_FSTmatALL <- pIND_FSTmatALL + 
  scale_fill_distiller(palette = "Spectral",
  direction = 1,
  limits = c(-0.001, 0.32)) + 
  scale_y_discrete(labels = poporder2lab,
                   breaks = poporder2,
                   position = "right", 
                   expand = c(0,0)) +
  scale_x_discrete(labels = poporder2lab,
                   breaks = poporder2,
                   expand = c(0,0))


# Export plot to png
png("results/Figs_paper/SI/IND_FST.png", width = 8.3, height = 10, units = "in", res = 600)
pIND_FSTmatALL
dev.off()

## sNMF -------------------------------------------------------------
# Plot cross entropy plot
# Plot k=2-5
# Convert genlight to matrix of 0,1,2 -------------------------------
mat_gl <- t(as.matrix(gl))
# Write table to .geno file which is essentially a matrix -----------
genofile <- "data/processed/LEA/02_IND/IND_hwe.geno"
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
Krange <- 1:5
snmf_IND_hwe <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2, 
                     entropy = T, alpha = 100, project = "new")

plot(snmf_IND_hwe, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), IND dataset, after HWE filter.")
# Calculate mean entropy
crossentropy.matIND <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_IND_hwe, K = x)})))
mean.entropyIND   <- apply(crossentropy.matIND, MARGIN=1, FUN=mean, na.rm=TRUE)
if(any(diff(mean.entropyIND)>0)){
  bestK <- unname(which(diff(mean.entropyIND)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
# SI: Export mean cross-entropy plot --------------------------------
png("results/Figs_paper/SI/IND_sNMF_meanxentropy_k1-5.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropyIND, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()


q.IND.ls <- merge_sNMF_qmatrix_multiK(snmf_IND_hwe, gl, k2plot = 1:5)
pop_label_IND <- gl$other$ind.metrics[, c("ID","IND_pop_recode")]
poporder <- c("Gujarat", "Maharashtra subpop. A", "Maharashtra", "Madhya Pradesh", "Karnataka",
              "Andhra Pradesh", "Tamil Nadu", 
              "Uttar Pradesh", "Odisha", "West Bengal")
plotQ(q.IND.ls[2:5],imgoutput="join",
      showdiv = T,
      showindlab=F,
      grplab=pop_label_IND[, c("IND_pop_recode"), drop = F],
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
      outputfilename="results/Figs_paper/SI/IND_sNMF_BCF_k1-5",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on IND dataset (after HWE filter), K = 2-5, 10 repetitions each"
)

