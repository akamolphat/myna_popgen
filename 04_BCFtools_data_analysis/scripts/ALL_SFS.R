# This script makes plots SFS for publication for the ALL dataset
#
# This script is for the BCFtools dataset
#
# Libraries ---------------------------------------------------------
library(dartR)
library(tidyverse)
source("../shared_scripts/functions.R")

# Define input files ------------------------------------------------
vcffileALL <- "../01_download_data/BCFtools_ALL/variant_calls.ALL.bialminQ30minGQ30DP15-125.norep.noadm.highnegfis.lmiss20.vcfthin.snps.vcf.gz"
# vcffileALL is the path to the vcf file that has not been filtered for singletons and doubeltons. This path will have to be changed to match the file path 
metadtfile <- "../01_download_data/TableS1.2v2.csv"
# Create folders to store some data/outputs -------------------------
dir.create("data/processed/LEA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/PCA/03_ALL/", recursive = T, showWarnings = F)
dir.create("results/Figs_paper/SI/", recursive = T, showWarnings = F)
# Read in VCF file --------------------------------------------------
glALL <- gl.read.vcf(vcffileALL)
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile, na.strings = "n/a")
# Subset metadt for the genlight object
indiv_names <- gsub('[a-z]$', '', indNames(glALL))
metadtsubALL <- metadt %>% 
  filter(ID %in% indiv_names) %>%
  arrange(factor(ID, levels = indiv_names)) 
# Attach individual metadata to genlight object ---------------------
glALL@other$ind.metrics <- metadtsubALL
# Assign individual ID to individual metadata
glALL@other$ind.metrics$ID <- indNames(glALL)
# Attach population information to genlight object
pop(glALL) <- glALL$other$ind.metrics$popdef2
# Order of SFS plots for main text ----------------------------------
poporder <- c("IND: Other", "IND: Maharashtra subpopulation A",
              "AUS: Melbourne", "Fiji", "NZ: Napier", "NZ: Other", 
              "AUS: Sydney", "AUS: Gold Coast",
              "Hawaii", "South Africa")

sfsdt <- plot_dt_SFS_compare_same_n_conf(glALL, 
                                         popname_ls = poporder,
                                         iterations = 100)

sfsdt$median_quantile_plot_nLoc

# Define colours, same as PCA plot

lab_order <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
               "NZ: Napier", "Fiji", "AUS: Melbourne", 
               "IND: Maharashtra subpopulation A", "IND: Other", 
               "Hawaii", "South Africa")
lab_order_renamed <- c("AUS: Gold Coast", "AUS: Sydney", "NZ: Other", 
                       "NZ: Napier", "Fiji", "AUS: Melbourne", 
                       "IND: Maharashtra subpop. A", "IND: Other", 
                       "Hawaii", "South Africa")
dtcolshapes <- data.frame(pop = lab_order,
                          pop_lab = lab_order_renamed,
                          col = c("#a6cee3", "#b2df8a", "#cab2d6", 
                                  "#ff7f00", "#33a02c", "#1f78b4",
                                  "#e31a1c", "#fdbf6f", "#fb9a99",
                                  "#6a3d9a"),
                          shape = c(16, 3, 12, 18, 8, 14, 11, 17,15,7))
## Figure 6 ---------------------------------------------------------
png("results/Figs_paper/Main/SFS_ALL.png", width = 8.3, height = 3.5, units = "in", res = 600)
sfsdt$median_quantile_plot_nLoc + 
  facet_grid(.~popn) +
  theme(axis.text.x = element_text(angle = 90, colour = "black", size = 8, vjust = 0.5),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        strip.text.x = element_blank()) +
  scale_fill_manual(values = dtcolshapes$col,
                    breaks = dtcolshapes$pop,
                    labels = dtcolshapes$pop_lab) +
  scale_x_continuous(breaks = c(0, 25, 50),
                     labels = c(0, 25, 50)) +
  xlab("Minor allele frequency (%)") +
  ylab("Median number of SNPs") +
  ggtitle("SFS plot, subset n = 6, 100 iterations")
dev.off()

# Figure S9.1: SFS on popdef1 ---------------------------------------
pop(glALL) <- glALL$other$ind.metrics$popdef1
glALLsub <- gl.drop.pop(glALL, pop.list = c("Gujarat", "Auckland", "Taupo (ROM)", "Kaikohe (ROM)"))

sfsdt2 <- plot_dt_SFS_compare_same_n_conf(glALLsub, 
                                          popname_ls = popNames(glALLsub),
                                          iterations = 100)
sfspoporder <- sfsdt2$datasum %>% 
  group_by(popn) %>%
  dplyr::summarise(loc_max = max(nLoc_median)) %>%
  arrange(desc(loc_max))

dtsumall <- sfsdt2$datasum 
dtsumall$popn <- factor(dtsumall$popn, 
                        levels = sfspoporder$popn)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dtcolall <- data.frame(values = gg_color_hue(dim(sfspoporder)[1]),
                       label_order = sfspoporder$popn)
plotSFS <- function(dtcombsum, labelorder, colourorder){
  dtcombsum <- dtcombsum %>% 
    filter(popn %in% labelorder)
  dtcombsum$popn <- factor(dtcombsum$popn, 
                          levels = labelorder)
  p3 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = nLoc_median, ymin = nLoc_min, 
                                               ymax = nLoc_max, lower = nLoc_quantile25,
                                               middle = nLoc_median, upper = nLoc_quantile75, 
                                               fill = popn))+
    scale_fill_manual(values = colourorder) +
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median number of SNPs (25%, 75% and min/max quantiles)") +
    # ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  return(p3)
}

p1 <- plotSFS(dtsumall, sfspoporder$popn[1:7], colourorder = dtcolall$values[1:7])
p2 <- plotSFS(dtsumall, sfspoporder$popn[8:14], colourorder = dtcolall$values[8:14])
p3 <- plotSFS(dtsumall, sfspoporder$popn[15:21], colourorder = dtcolall$values[15:21])
p4 <- plotSFS(dtsumall, sfspoporder$popn[22:28], colourorder = dtcolall$values[22:28])

sfstheme <-  theme(axis.title = element_blank(),
                     axis.text.x = element_blank(), 
                     legend.position = "none")

png("results/Figs_paper/SI/ALL_SFS.png", width = 8.3, height = 8.3, res = 600, units = "in")
cowplot::plot_grid(p1 + facet_grid(.~popn) + ylim(0, max(dtsumall$nLoc_max)) + sfstheme, 
                   p2 + facet_grid(.~popn) + ylim(0, max(dtsumall$nLoc_max)) + sfstheme,
                   p3 + facet_grid(.~popn) + ylim(0, max(dtsumall$nLoc_max)) + sfstheme,
                   p4 + facet_grid(.~popn) + ylim(0, max(dtsumall$nLoc_max)) + sfstheme,
                   ncol = 1,
                   nrow = 4)
dev.off()

