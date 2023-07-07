# This script makes plots for the figure for publication for the NZ dataset
#
# This script is for the BCFtools dataset
#
# Please note that the output plots are very slightly different but 
# showed the same results.
#
# Likely due to a mistake of re-reading an older file containing the
# principal compenents made with a different input file
# 
# Define input files ------------------------------------------------
#
vcffile <- "../01_download_data/BCFtools_NZ/variant_calls.NZ.bialminQ30minGQ30DP15-100.norep.lmiss20.nosingledoubletons.vcfthin.hwe.snps.vcf.gz"
metadtfile <- "../01_download_data/TableS1.2v2.csv"
# Define .geno file to store data for performing LEA ----------------
genofile <- "data/processed/LEA/01_NZ/NZ_BCFtools_hwe.geno"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/01_NZ/", recursive = T, showWarnings = F)
dir.create("results/PCA/01_NZ/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/Main/", recursive = T, showWarnings = F)



# Load libraries ----------------------------------------------------
library(dartR)
library(LEA)
library(tidyverse)
library(cowplot)
library(vcfR)
library(xlsx)
library(pheatmap)
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
pop(gl) <- gl@other$ind.metrics$popdef1
# Count samples per location ----------------------------------------
gl@other$ind.metrics %>% dplyr::count(popdef1)

# Re-perform PCA on the full dataset --------------------------------
pc <- gl.pcoa(gl, nfactors=10, parallel = T, n.cores = 8)
# Quick plot to see variance explained ------------------------------
gl.pcoa.plot(pc, gl, pop.labels= "pop") 
# Create dataframe for pricipal component axes ----------------------
dtpc <- data.frame(pc$scores)
save(pc, file = "results/PCA/01_NZ/PCA_NZ_BCFtools_hwe.Rdata")
load("results/PCA/01_NZ/PCA_NZ_BCFtools_hwe.Rdata")
# Assign population information -------------------------------------
dtpc$popdef1 <- gl$other$ind.metrics$popdef1
# Assign whether it is ROM or Modern
dtpc <- dtpc %>% 
  mutate(ROM = ifelse(grepl("(ROM)", popdef1), 
                                "ROM", 
                                "Modern")) %>%
  mutate(ROM = replace(ROM, "(ROM)" %in% popdef1, "ROM samples")) %>%
# Assign location regardless of whether it is ROM or not
# Note that Napier (ROM) is changed to "Napier " with a space at the end
# This is for plotting purpuses
  mutate(location = replace(popdef1, popdef1 == "Kaikohe (ROM)", "Kaikohe")) %>%
  mutate(location = replace(location, location == "Auckland", "Auckland (AKL)")) %>%  # NOTE that in Appendix figure, this is Onehunga
  mutate(location = replace(location, location == "Auckland (ROM)", "Auckland (AKR)")) %>%
  mutate(location = replace(location, location == "Hamilton (ROM)", "Hamilton")) %>%
  mutate(location = replace(location, location == "Taupo (ROM)", "Taupo")) %>%
  mutate(location = replace(location, location == "Napier (ROM)", "Napier ")) 
  
  
write.csv(dtpc, file = "results/PCA/01_NZ/PCA_NZ_BCFtools_hwe.csv")
# Get lat and long for each location to sort legends by
# This was not used in the end
dtlatlon <- gl@other$ind.metrics %>% 
  group_by(popdef1) %>%
  summarise(lat = mean(latitude),
            lon = mean(longitude))
# dtlatlon$popdef1[order(dtlatlon$lon)]
lab_order <- c("Helensville", "Ngunguru", "Waitakeres", 
               "Auckland (AKL)", "Leigh", "Waiheke",
               "Great Barrier Island", "Thames", "Napier", 
               "Kaikohe","Auckland (AKR)", 
               "Hamilton", "Taupo", "Napier ")
lab_order_abb <- c("Helensville (HEL)", "Ngunguru (NGU)", "Waitakeres (WTK)", 
                   "Auckland (AKL)", "Leigh (LEI)", "Waiheke (WHK)",
                   "Great Barrier Island (GBI)", "Thames (THA)", "Napier (NAP)", 
                   "Kaikohe (KKH)","Auckland (AKR)", 
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
# dtpc <- read.csv("results/PCA/01_NZ/PCA_NZ_BCFtools_hwe.csv")
dtpc$location <- factor(dtpc$location,
                        levels = lab_order)

# Create PCA related SI plots ---------------------------------------
## Figure S8.1: Screeplot -------------------------------------------
screep <- PCA_screeplot(pc, nPC = 15)
png("results/Figs_paper/SI/NZ_PCA_BCFtools_hwe_screeplot.png", width = 8.3, height = 3, units = "in", res = 600)
screep
dev.off()

## Figure S8.2: PC1v2-3 ---------------------------------------------
# Calculate variance explained for each PC 
dtvar <- PCA_var_explained(pc)
## PC1 vs PC2
pl <- plotPCASI(dtpc, 
                varls = round(dtvar$variance_explained, 1),
                dtcolshapes,
                xaxis = 1,
                yaxis = 2, popcol = "location") + 
  PCA_theme

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
# Output to PNG
png("results/Figs_paper/SI/NZ_PCA_BCFtools_hwe.png", width = 8.3, height = 4.2, units = "in", res = 600)
plot_grid(pca_SI, legend_b, ncol = 1, rel_heights = c(1, .15))
dev.off()

## Figure S8.3: PC loadings -----------------------------------------
dtpcload <- data.frame(pc$loadings)
dtpcload$contig <- as.numeric(gl$other$loc.metrics$CHROM)
dtpcload$position <- as.numeric(gl$other$loc.metrics$POS)
## Order by contig and position. It is already ordered by CONTIG
# and POSITION but doing this again to make sure
dtpcload <- dtpcload %>% arrange(contig, position)

png("results/Figs_paper/SI/NZ_PCA_BCF_loadings.png", width = 8.3, height = 5, units = "in", res = 600)

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


# Create Isolation by distance plot ---------------------------------
# Attach latlong data to gl$other$latlong
latlon_mat <- as.matrix(gl$other$ind.metrics[,c("latitude", "longitude")])
colnames(latlon_mat) <- c("lat", "lon")
rownames(latlon_mat) <- indNames(gl)
gl$other$latlon <- latlon_mat

# Performs IBD using our own distance matrix ------------------------
# This is because the genetic distance matrix can be changed
# The default does not use bootstrapping to calculate FST
#
# Generally the results are the same, but this done for completeness
# purpose for the final plot
#
# Retain populations with sufficient n (n > 5)
#
# Note that there is one "Auckland" sample which I have grouped with
# the Auckland Whitford landfill samples. However, we will also exclude 
# this one sample due to the difference in coordinates. This is sample
# M0367
#
gl$other$ind.metrics %>% dplyr::count(popdef1)
# Drop populations and samples 
glsub <- gl.drop.pop(gl, pop.list = c("Auckland", "Taupo (ROM)", "Kaikohe (ROM)"))#, "Great Barrier Island"))
glsub <- gl.drop.ind(glsub, ind.list = c("M0367"))
glsub <- gl.drop.pop(glsub, pop.list = c("Napier", "Napier (ROM)"))
# Covert coordinates to Mercator projection
coords <- dismo::Mercator(glsub@other$latlon[, c("lon", "lat")])
# Create distances between each population, using the mean 
# coordinate of each population. All samples from the same
# population have the same coordinates
pop.xy <- apply(coords, 2, function(a)
    tapply(a, pop(glsub), mean, na.rm = T))
## Create distance matrix in metres ---------------------------------
Dgeo <- as.dist(dist(pop.xy))
## Create FST matrix ------------------------------------------------
## Standard
fst <- gl.fst.pop(glsub, nboots = 100)
Dgen <- as.dist(fst$Fsts)
# UNUSED: Reich's FST 
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

save(Dgen, Dgen2, Dgeo, file = "results/Figs_paper/NZ_IBD_distance_matrix.Rdata")
load("results/Figs_paper/NZ_IBD_distance_matrix.Rdata")

# Checking if the results from doing it manually is similar to
# the default gl.ibd function
gl.ibd(Dgen = Dgen, Dgeo = Dgeo)
gl.ibd(Dgen = Dgen2, Dgeo = Dgeo)
gl.ibd(glsub)

## Test out different transformations -------------------------------
# The IBD is insignificant in all transformations
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

p <- gl.ibd(Dgen = Dgen, Dgeo = Dgeo, Dgen_trans='Dgen/(1-Dgen)', 
            Dgeo_trans = "log(Dgeo/1000)")
manteltest <- p$mantel 
lm_eqn <- function(df,
           r = manteltest$statistic,
           pp = manteltest$signif) {
  # Function from the source code of gl.ibd
  # This is to extract the equation from gl.ibd outputs
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
    label = lm_eqn(res,
                   p$mantel$statistic,
                   p$mantel$signif),
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
        axis.text = element_text(size = 8, face = "bold"))#,
        # aspect.ratio = 1)

## Figure 2A+B: PCA and mantel test ---------------------------------
#
# Note that Figure 2A is slightly different from that of the main text
# This is due to a slight mistake in input file. The axis have also 
# been inverted as often happens with PCA.
pl <- pl + 
  guides(colour=guide_legend(nrow=9,byrow=F), 
         shape=guide_legend(nrow=9,byrow=F)) +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, face = "bold"),
        legend.key.size = unit(0.75, 'lines'),
        legend.box.margin = margin(-10,-10,-10,-10),
        # legend.position = c(0.975,0.05),  
        # legend.justification = c(1,0),
        legend.position = c(0.975,0.975),
        legend.justification = c(1,1),
        legend.background = element_rect(color = "black", size = 0.01, linetype = "solid"))#,
        # aspect.ratio = 1)

prow <- plot_grid(pl, p3, ncol = 2, rel_widths = c(4,4), labels = NULL)

prow
png("results/Figs_paper/Main/NZ_pop_str_PCA_mantel.png", width = 8.3, height = 4, units = "in", res = 600)
prow
dev.off()

## Figure S8.4: IBD without GBI -------------------------------------
Dgenmat <- as.matrix(Dgen)
Dgen_noGBI <- Dgenmat[(rownames(Dgenmat) != "Great Barrier Island"),
                      (colnames(Dgenmat) != "Great Barrier Island")]
Dgeomat <- as.matrix(Dgeo)
Dgeo_noGBI <- Dgeomat[(rownames(Dgeomat) != "Great Barrier Island"),
                      (colnames(Dgeomat) != "Great Barrier Island")]
Dgen_noGBI <- as.dist(Dgen_noGBI)
Dgeo_noGBI <- as.dist(Dgeo_noGBI)
pIBD_noGBI <- gl.ibd(Dgen = Dgen_noGBI, Dgeo = Dgeo_noGBI, Dgen_trans='Dgen/(1-Dgen)', Dgeo_trans = "log(Dgeo/1000)")
manteltest_noGBI <- pIBD_noGBI$mantel 
res2 <- data.frame(Dgen = as.numeric(pIBD_noGBI$Dgen), Dgeo = as.numeric(pIBD_noGBI$Dgeo))
p_noGBI <- ggplot(res2, aes(x = Dgeo, y = Dgen)) + geom_point() + 
  geom_smooth(method = "lm", se = TRUE, colour = "black") + 
  ylab("Dgen_trans") + xlab("Dgeo_trans") +
  annotate(
    "text",
    label = lm_eqn(res2,
                   pIBD_noGBI$mantel$statistic,
                   pIBD_noGBI$mantel$signif),
    x = Inf,
    y = -Inf,
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.5,
    size = 2.5
  ) + 
  theme_classic() + ylab(expression(italic(F[ST])/(1-italic(F[ST])))) + 
  xlab("ln(distance in km)") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8, face = "bold"))

png("results/Figs_paper/SI/NZ_IBD_noGBI.png", width = 8.3, height = 8.3, units = "in", res = 600)
p_noGBI
dev.off()

# sNMF --------------------------------------------------------------
## Convert genlight to matrix of 0,1,2 ------------------------------
mat_gl <- t(as.matrix(gl))
## Write table to .geno file which is essentially a matrix ----------
write.table(mat_gl, file = genofile,
            quote = F, sep = "", na = "9", 
            row.names = F, col.names = F)
## Perform sNMF -----------------------------------------------------
Krange <- 1:5
snmf_NZ_hwe <- snmf(genofile, K = Krange, repetitions = 10, ploidy = 2, 
                         entropy = T, alpha = 100, project = "new")
# snmf_NZ_hwe <- load.snmfproject("data/processed/LEA/01_NZ/NZ_BCFtools_hwe.snmfProject")
## Plot cross-entropy
plot(snmf_NZ_hwe, col = "blue4", cex = 1.4, pch = 19, main = "Cross entropy plot (snmf), NZ, after HWE filter.")

## Figure S8.6: Mean cross-entropy plot-------------------------------
# Calculate mean entropy 
crossentropy.mat <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf_NZ_hwe, K = x)})))
mean.entropy   <- apply(crossentropy.mat, MARGIN=1, FUN=mean, na.rm=TRUE)
# Make plot
png("results/Figs_paper/SI/NZ_sNMF_meanxentropy_k1-5.png", width = 8.3, height = 4.5, units = "in", res = 600)
plot(mean.entropy, xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy (10 repetitions)")
dev.off()

# Find best K value with minimum mean cross-entropy
if(any(diff(mean.entropy)>0)){
  bestK <- unname(which(diff(mean.entropy)>0)[1])
} else {
  bestK <- unname(Krange[1])
}
bestK
# plot k = 2
rep <- 1:10
samplenames <- indNames(gl)
numind <- length(samplenames)
pop_label_sub <- gl$other$ind.metrics[, c("ID","popdef1", "popdef2")]
pop_labels <- c("Auckland (ROM)", "Great Barrier Island",
                "Hamilton (ROM)", "Helensville", 
                "Kaikohe (ROM)", "Leigh",
                "Ngunguru", "Auckland", 
                "Taupo (ROM)", "Thames",
                "Waiheke", "Waitakeres",
                "Napier", "Napier (ROM)")
dtlatlon$popdef1[order(dtlatlon$lon)]
pop_trans <- data.frame(popdef1 = pop_labels,
                        loc_time = LETTERS[1:length(pop_labels)],
                        loc_time_abb = c("AKR", "GBI", "HAM", "HEL", "KKH", "LEI",
                                         "NGU", "AKL", "TAU", "THA", "WHK", "WTK", "NAP", "NAR"))
pop_label_sub <- merge(pop_label_sub, pop_trans)
colnames(pop_label_sub) <- c("popdef1", "ID", "population", "location/sample", "loc_time_abb")
pop_label_sub <- pop_label_sub[match(indNames(gl), pop_label_sub$ID),]

q.df.ls <- merge_sNMF_qmatrix_multiK(snmf_NZ_hwe, gl, k2plot = 2)

pp <- plotQ(q.df.ls,
            clustercol = c("#cab2d6", "#ff7f00"),
            # imgoutput="join",
            showdiv = T,
            showindlab=F,
            grplab=pop_label_sub[, c("loc_time_abb", "population"), drop = F],
            # subsetgrp=c("NZ: Other", "NZ: Napier"),
            subsetgrp=c("KKH", "HEL", "NGU", "WTK", "AKL", "LEI", 
                        "AKR", "WHK", "HAM", "GBI", "THA", "TAU", "NAP", "NAR"),
            grplabangle = 90,
            grplabheight = 3,
            grplabsize = 3,
            grplabpos = 0.65,
            grplabcol = "black",
            selgrp="loc_time_abb",
            splab = c("K=2"),
            # panelratio = c(1,1),
            # grplabspacer = -2,
            ordergrp=T,
            showlegend=F,
            showtitle=F,
            showsubtitle=F,
            height=4,
            width =21,
            linesize = 0.5,
            # indlabspacer=-1,
            barbordercolour=NA,barbordersize=0,
            legendrow = 1,
            showsp = F,
            legendlab = c("K1", "K2"),
            outputfilename="results/Figs_paper/Main/LEA_NZ_hwe_k2_merged",
            imgtype="png",
            exportpath=getwd(),
            titlelab="New Zealand Common myna population structure",
            returnplot = T#,
            # subtitlelab="sNMF on subsetted dataset, n = 10 per population."
)

## Figure 2: PCA, Mantel's test and sNMF to png ---------------------
prow2 <- plot_grid(prow, pp$plot[[1]], nrow = 2, rel_heights = c(4,2), labels = NULL)
prow2
png("results/Figs_paper/Main/NZ_pop_str.png", width = 8.3, height = 6, units = "in", res = 600)
prow2
dev.off()

## Figure S8.7: sNMF K = 2-5 -----------------------------------------
q.df.ls <- merge_sNMF_qmatrix_multiK(snmf_NZ_hwe, gl, k2plot = 1:5)
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
      outputfilename="results/Figs_paper/SI/NZ_sNMF_BCF_k1-5",imgtype="png",
      exportpath=getwd(),
      titlelab="Common myna population structure",
      subtitlelab="sNMF on NZ dataset, K = 2-5, 10 repetitions each"
)
