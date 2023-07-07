# This script plot maps for the manuscript --------------------------#####
# Libraries ---------------------------------------------------------#####
library(ggplot2)
library(rgdal)
library(plyr)
library(tidyverse)
library(gridExtra)
# library(gtable)
# library(grid)
library(cowplot)  
library(ggrepel)  # Not used in final version of making the plots as 
                  # labels were added manually
# Define inputs -----------------------------------------------------#####
# Paths to files may change depending on where it is stored in the computer
path2metadatafile <- "../01_download_data/TableS1.2v2.csv"  # This path may change
path2shapefile <- "../01_download_data/A_tristis_distribution_map/data_0.shp"  # This path will depend on where the shapefile is stored.

# Read metadata to get lat/lon --------------------------------------#####
meta_dt <- read.csv(path2metadatafile, header = T, na.strings = "n/a")

# Add country of origin ---------------------------------------------#####
meta_dt$COUNTRY.OF.ORIGIN <- meta_dt$popdef2
meta_dt$popdef2 %>% unique()
meta_dt <- meta_dt %>%
  mutate(COUNTRY.OF.ORIGIN = case_when(
    popdef2 %in% c("NZ: Other", "NZ: Napier") ~ "New Zealand",
    popdef2 %in% c("IND: Other", "IND: Maharashtra subpopulation A") ~ "India",
    popdef2 %in% c("AUS: Sydney", "AUS: Gold Coast", "AUS: Melbourne", NA) ~ "Australia",
    popdef2 == "Fiji" ~ "Fiji",
    popdef2 == "Hawaii" ~ "USA",
    popdef2 == "South Africa" ~ "South Africa"))

# Add more precise location for plotting ----------------------------#####
# These are pop_def2 locations 
# when defined and original Australian locations as defined by 
# Ewart et al. (2019)
meta_dt$pop_clean <- meta_dt$Precise.location
meta_dt$pop_clean[meta_dt$COUNTRY.OF.ORIGIN %in% c("Australia") & is.na(meta_dt$popdef1)] <- meta_dt$popdef1[meta_dt$COUNTRY.OF.ORIGIN %in% c("Australia") & is.na(meta_dt$popdef1)]
meta_dt$pop_clean[meta_dt$COUNTRY.OF.ORIGIN %in% c("New Zealand")] <- meta_dt$popdef1[meta_dt$COUNTRY.OF.ORIGIN %in% c("New Zealand")]

meta_dt$pop_clean[meta_dt$COUNTRY.OF.ORIGIN %in% c("USA", "South Africa", "Fiji")] <- meta_dt$popdef2[meta_dt$COUNTRY.OF.ORIGIN %in% c("USA", "South Africa", "Fiji")]
meta_dt$pop_clean[meta_dt$COUNTRY.OF.ORIGIN %in% c("India")] <- meta_dt$popdef1[meta_dt$COUNTRY.OF.ORIGIN %in% c("India")]
meta_dt$pop_clean[meta_dt$pop_clean == "Melbourne, Nunawading Recycling Centre"] <- "Melbourne"
meta_dt$pop_clean[is.na(meta_dt$latitude)] <- NA

# Remove duplicates -------------------------------------------------#####
# Keeping only one of our inhouse replicates and one of Ewart et al. 
# inhouse replicates
inhouse_rm_ls <- c("M0092", "M0093", "M0111", "M0186", "M0187", "M0370", "M0188", "M0280", "M0281", "M0340", "M0368", "M0369", "M0236")
kyle_inhouse_ls <- grep("^[0-9]*[a-zA-Z]$", meta_dt$ID, value = T)
inhouse_ls <- c(inhouse_rm_ls, kyle_inhouse_ls)
dim(meta_dt[meta_dt$ID %in% inhouse_ls,])
meta_dt <- meta_dt[!meta_dt$ID %in% inhouse_ls,]

# Calculate average lat/lon for plotting ----------------------------#####
meta_dt_pop <- meta_dt %>% 
  dplyr::group_by(COUNTRY.OF.ORIGIN, pop_clean) %>%
  dplyr::summarise(latitude = mean(latitude, na.rm = T),
            longitude = mean(longitude, na.rm = T),
            n = n())
meta_dt_pop$longitude_360 <- meta_dt_pop$longitude
meta_dt_pop$longitude_360[meta_dt_pop$longitude < 0 & !is.na(meta_dt_pop$longitude)] <- 360 + meta_dt_pop$longitude[meta_dt_pop$longitude < 0 & !is.na(meta_dt_pop$longitude)]
# Get an average for each country for easier plot
meta_dt_all <- meta_dt_pop %>% 
  dplyr::group_by(COUNTRY.OF.ORIGIN) %>%
  dplyr::summarise(latitude = mean(latitude, na.rm = T),
            longitude = mean(longitude_360, na.rm = T),
            n = n())
meta_dt_all$COUNTRY.OF.ORIGIN[meta_dt_all$COUNTRY.OF.ORIGIN == "USA"] <- "Hawaii"
meta_dt_all$label_n <- paste(meta_dt_all$COUNTRY.OF.ORIGIN, ", n = ", meta_dt_all$n, sep = "")
# Make map ----------------------------------------------------------#####
## Read shape file for Common myna distribution ---------------------#####
shapefile <- readOGR(dsn=path.expand(path2shapefile))  # This path will depend on where the shapefile is stored.
shapefile_df <- fortify(shapefile)
### Alter shape file longitude to 0-360. 
# I can only do this with the distribution file and not the shoreline 
# shape file because there are some points in the polygons joining 
# different sides of the map. Using map_data() function in ggplot2 
# instead
#
shapefile_df$long.edit <- shapefile_df$long
shapefile_df$long.edit[shapefile_df$long < 0] <- 360 + shapefile_df$long[shapefile_df$long < 0]
shapefile_df$ranges <- shapefile_df$id
shapefile_df$ranges[shapefile_df$id == 0] <- "Native"
shapefile_df$ranges[shapefile_df$id == 1] <- "Invasive"

## Make base map ----------------------------------------------------#####
# Make base map of the continents and the distribution
mapWorld <- map_data('world', wrap=c(-20,340), ylim=c(-60,90))
map_dist <- geom_polygon(data = shapefile_df, aes(x = long.edit, y = lat, group = group, fill = ranges), alpha = 1) 
latbrks <- seq(-60, 90, by = 30)
lonbrks <- seq(0, 330, by = 60)
lonlabs <- lonbrks
lonlabs[lonlabs > 180] <- lonlabs[lonlabs > 180] - 360

mp <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey") +
  map_dist  +
  theme_bw() +
  theme(legend.justification = c(1,1), 
        legend.position = c(0.99,0.99),
        legend.direction = "horizontal",
        # legend.justification = "bottom", 
        # legend.position = "bottom", 
        legend.text=element_text(size=10),
        legend.title=element_text(size=12, face = 'bold'),
        legend.text.align = 0,
        legend.background = element_rect(color = "black", size = 0.01, linetype = "solid"),
        axis.text=element_text(colour = 'black', size=16),
        axis.line = element_line(colour = 'black'),
        axis.title = element_blank(),
        axis.text.x = element_text(margin = margin(3, unit = "pt")),
        plot.title = element_text(margin = margin(b = -25, unit = "pt"), size = 20),
        plot.margin = unit(c(0,20,0,5), 'pt'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "azure",
                                        colour = "azure",
                                        size = 0.5, 
                                        linetype = "solid")) +
  coord_fixed(xlim = c(-20,340), ylim = c(-60, 90), expand = F) + 
  scale_y_continuous(breaks = latbrks, 
                     labels = paste(latbrks, '°', sep = ''),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = lonbrks, 
                     labels = paste(lonlabs, '°', sep = ''),
                     expand = c(0, 0)) +
  labs(fill = NULL) +
  scale_fill_manual(values = c("#d95f02","#377eb8")) 

## Make map for ALL dataset -----------------------------------------
### Alter meta_dt_all to annotation location ------------------------
# Note that this is not necessary. However, this helps with visualising
# the map. For the final version of the figure, the annotations
# were added manually in powerpoint to be able to dictate the location
# of the labels more accurately.
meta_dt_label <- meta_dt_all
meta_dt_label$latitude[meta_dt_label$COUNTRY.OF.ORIGIN == "India"] <- 15
meta_dt_label$longitude[meta_dt_label$COUNTRY.OF.ORIGIN == "India"] <- 58

meta_dt_label$latitude[meta_dt_label$COUNTRY.OF.ORIGIN == "South Africa"] <- -30
meta_dt_label$longitude[meta_dt_label$COUNTRY.OF.ORIGIN == "South Africa"] <- 50

meta_dt_label$latitude[meta_dt_label$COUNTRY.OF.ORIGIN == "Hawaii"] <- 17
meta_dt_label$longitude[meta_dt_label$COUNTRY.OF.ORIGIN == "Hawaii"] <- 190

meta_dt_label$latitude[meta_dt_label$COUNTRY.OF.ORIGIN == "Fiji"] <- -12 
meta_dt_label$longitude[meta_dt_label$COUNTRY.OF.ORIGIN == "Fiji"] <- 185

meta_dt_label$latitude[meta_dt_label$COUNTRY.OF.ORIGIN == "New Zealand"] <- -45 
meta_dt_label$longitude[meta_dt_label$COUNTRY.OF.ORIGIN == "New Zealand"] <- 180

meta_dt_label$latitude[meta_dt_label$COUNTRY.OF.ORIGIN == "Australia"] <- -27 
meta_dt_label$longitude[meta_dt_label$COUNTRY.OF.ORIGIN == "Australia"] <- 125

### Make map
mp2a <- mp + 
  geom_point(data = meta_dt_pop, aes(x = longitude_360, y = latitude, size = n), alpha = 0.5) +
  # guides(size = "none") + 
  coord_fixed(xlim = c(20, 210), ylim = c(-50,50)) +
  theme(#legend.text=element_text(size=15),
    # legend.position = "bottom",
    # legend.justification = c(0,0),
    legend.box.just = "right") + #,
  # legend.title=element_blank()) +
  scale_fill_manual(values = c("#d95f02","#377eb8")) +
  scale_color_manual(values = c("#d95f02", "#377eb8")) +
  labs(size = "Number of samples")

mp2 <- mp2a + geom_text_repel(data = meta_dt_label, 
                              aes(x = longitude, y = latitude, label = label_n),
                              force = 0.5,
                              size = 3.5)



### Step 1
# Draw a plot with the colour legend
(p1 <- mp + scale_fill_manual(values = c("#d95f02","#377eb8")) +
    scale_color_manual(values = c("#d95f02", "#377eb8")) + 
    theme(legend.text = element_text(size = 12)))

# Extract the colour legend - leg1
leg1 <- get_legend(p1)

### Step 2
# Draw a plot with the size legend
(p2 <- ggplot(data = meta_dt_pop, aes(x = longitude_360, y = latitude, size = n)) +
    geom_point(alpha = 0.5) + 
    scale_size_continuous(breaks = c(0,5,10,20,30,40,50), 
                          labels = c(0,5,10,20,30,40,50)) +
    labs(size = "No. of indiv.") +
    theme_bw() + 
    guides(size = guide_legend(nrow = 1, byrow=T, title.position = "left")) +
    theme(legend.justification = c(1,1), 
          legend.position = c(0.99,0.99),
          legend.direction = "horizontal",
          legend.text=element_text(size=10),
          legend.title=element_text(size=12, face = 'bold'),
          legend.text.align = 0,
          legend.background = element_rect(color = "black", size = 0.01, linetype = "solid")))

# Extract the size legend - leg2
leg2 <- get_legend(p2) 

# Step 3
# Draw a plot with no legends - plot
(plot <- mp2 + guides(size="none") +
    theme(legend.position = c(0.975,0.95),
          legend.justification = c(1,1),
          axis.text = element_text(size = 8, colour = "black")))

png("results/Ourrange_datapoints_MS.png", width = 7, height = 3.5, units = "in", res = 600)
plot
dev.off()



png("results/Ourrange_datapoints_nolabelMS.png", width = 7, height = 3.5, units = "in", res = 600)
mp2a + guides(size="none") +
  theme(legend.position = c(0.975,0.95),
        legend.justification = c(1,1),
        axis.text = element_text(size = 8, colour = "black"))
dev.off()

# Make NZ map -------------------------------------------------------#####
latbrks <- seq(-41, -34, by = 2)
lonbrks <- seq(173, 179, by = 2)
meta_dt_NZ <- subset(meta_dt_pop, COUNTRY.OF.ORIGIN == "New Zealand")
meta_dt_NZ <- meta_dt_NZ[rev(order(meta_dt_NZ$latitude)),]
meta_dt_NZ$labx <- c(178.5, 178.5, 178.5, 178.5, 173, 178.5, 173, 173, 173, 178.5, 173, 173, 178.5, 178.5)
mp3a <- mp + 
  geom_point(data = meta_dt_NZ, aes(x = longitude_360, y = latitude, size = n), alpha = 0.5) +
  guides(fill = "none", colour = "none") +
  # geom_text_repel(data = meta_dt_NZ,
  #                 aes(x = longitude_360, y = latitude, label = loc_precise_time),
  #                 nudge_x = -0.2, direction = "y", hjust = "right",
  #                 size = 3.5
  # ) +
  # geom_text(data = meta_dt_NZ,
  #           aes(x = labx, y = latitude, label = loc_precise_time),
  #           hjust = 1,
  #           # max.overlaps = Inf,
  #           # box.padding = 0.5,
  #           # force = 250,
  #           size = 3.5) +
  # annotate("segment", x = c(meta_dt_NZ$labx), y = c(meta_dt_NZ$longitude_360),
  #          xend = c(meta_dt_NZ$latitude), yend = c(meta_dt_NZ$longitude_360)+1,
  #          # arrow = arrow(angle = 20, type = "closed"),
  #          colour = "black") +
  # annotate("segment", x = meta_dt_NZ$labx, y = meta_dt_NZ$latitude,
  #          xend = meta_dt_NZ$longitude_360, yend = meta_dt_NZ$latitude,
  #          # arrow = arrow(angle = 20, type = "closed"),
  #          colour = "black") +
  coord_fixed(xlim = c(172.5, 178.7), ylim = c(-41.7,-34.2)) +
  scale_y_continuous(breaks = latbrks, 
                     labels = paste(latbrks, '°', sep = ''),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = lonbrks, 
                     labels = paste(lonbrks, '°', sep = ''),
                     expand = c(0, 0)) + 
  theme(axis.text=element_text(colour = 'black', size=8))

mp3a <- mp3a + scale_size_continuous(breaks = c(0,5,10,20,30,40,50), 
                                   labels = c(0,5,10,20,30,40,50))
mp3 <- mp3a + geom_text_repel(data = meta_dt_NZ,
                              aes(x = longitude_360, y = latitude, label = pop_clean),
                              min.segment.length = 0,
                              max.overlaps = Inf,
                              box.padding = 0.5,
                              force = 250,
                              size = 3.5) 

# Make India --------------------------------------------------------#####
latbrks <- seq(15, 50, by = 15)
lonbrks <- seq(70, 115, by = 15)
meta_dt_IND <- subset(meta_dt_pop, COUNTRY.OF.ORIGIN %in% c("India", "New Zealand"))
meta_dt_INDlab <- subset(meta_dt_pop, COUNTRY.OF.ORIGIN %in% c("India"))
mp5a <- mp + 
  geom_point(data = meta_dt_IND, aes(x = longitude_360, y = latitude, size = n), alpha = 0.5) +
  guides(fill = "none", colour = "none") +
  coord_fixed(xlim = c(58, 100), ylim = c(5,50)) +
  scale_y_continuous(breaks = latbrks, 
                     labels = paste(latbrks, '°', sep = ''),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = lonbrks, 
                     labels = paste(lonbrks, '°', sep = ''),
                     expand = c(0, 0)) + 
  # theme(axis.text=element_text(colour = 'black', size=8)) +
  guides(fill = "none", colour = "none") +
  scale_size_continuous(breaks = c(0,5,10,20,30,40,50), 
                        labels = c(0,5,10,20,30,40,50)) +
  scale_y_continuous(breaks = latbrks, 
                     labels = paste(latbrks, '°', sep = '')) +
  scale_x_continuous(breaks = lonbrks, 
                     labels = paste(lonbrks, '°', sep = '')) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(colour = 'black', size=8, margin = margin(3, unit = "pt")),
        axis.line = element_line(colour = 'black'),
        axis.title = element_blank(),
        plot.title = element_text(margin = margin(b = -25, unit = "pt"), size = 20),
        plot.margin = unit(c(0,20,0,5), 'pt'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "azure",
                                        colour = "azure",
                                        size = 0.5,
                                        linetype = "solid")) 

mp5 <- mp5a + geom_text_repel(data = meta_dt_INDlab, 
                              aes(x = longitude_360, y = latitude, label = pop_clean), 
                              max.overlaps = Inf,
                              box.padding = 0.5,
                              force = 250, 
                              size = 3.5) 

prow <- plot_grid(mp3 + theme(legend.position = "none"), 
                  mp5 + theme(legend.position = "none"), ncol = 2, rel_widths = c(4,4), labels = NULL)

prow2 <- plot_grid(plot, prow, nrow = 2, rel_heights = c(1,1), labels=NULL)
prow3 <- plot_grid(prow2, leg2, nrow = 2, rel_heights = c(2, .15))

png("results/ALL_datapoints_legend.png", width = 8.3, height = 10, units = "in", res = 600)
prow3
dev.off()



#### Add legend
prow2 <- plot_grid(prow, leg2, ncol = 1, rel_heights = c(1, .15))

png("results/NZ_datapoints_nolabels.png", width = 7, height = 3.5, units = "in", res = 600)
mp3a + theme(legend.position = "none")
dev.off()

png("results/NZ_datapoints_MS.png", width = 7, height = 3.5, units = "in", res = 600)
mp3 + theme(legend.position = "none")
dev.off()
png("results/NZ_datapoints_MS_legend.png", width = 7, height = 3.5, units = "in", res = 600)
mp3 + guides(size = guide_legend(ncol=1,byrow=T, title.position = "top")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=10, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.spacing.y = unit(0.6, 'cm'))
dev.off()
png("results/India_datapoints_nolabels.png", width = 7, height = 3.5, units = "in", res = 600)
mp5a + theme(legend.position = "none")
dev.off()

png("results/India_datapoints_MS.png", width = 7, height = 3.5, units = "in", res = 600)
mp5 + theme(legend.position = "none")
dev.off()

png("results/Legend_only_MS.png", width = 7, height = 3.5, units = "in", res = 600)
plot(leg2)
dev.off()


# Australian datapoints ---------------------------------------------#####
latbrks <- seq(-40, -10, by = 5)
lonbrks <- seq(130, 160, by = 10)
meta_dt_AUS <- subset(meta_dt_pop, COUNTRY.OF.ORIGIN %in% c("Australia"))
kept_pops <- c("Melbourne", "Melbourne (ROM)", "Sydney", "Sydney (ROM)", "Gold Coast")
meta_dt_AUS$final_filt <- "No"
meta_dt_AUS$final_filt[meta_dt_AUS$loc_precise_time %in% kept_pops] <- "Yes"
mp6 <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(colour = 'black', size=8, margin = margin(3, unit = "pt")),
        axis.line = element_line(colour = 'black'),
        axis.title = element_blank(),
        plot.title = element_text(margin = margin(b = -25, unit = "pt"), size = 20),
        plot.margin = unit(c(0,20,0,5), 'pt'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "azure",
                                        colour = "azure",
                                        size = 0.5, 
                                        linetype = "solid")) +
  map_dist + 
  guides(fill = "none", colour = "none") +
  scale_fill_manual(values = c("#d95f02","#377eb8")) +
  coord_fixed(xlim = c(135, 155), ylim = c(-40,-9.5)) +
  geom_point(data = meta_dt_AUS, aes(x = longitude_360, y = latitude, size = n, col = final_filt), alpha = 0.5) + 
  scale_size_continuous(breaks = c(0,5,10,20,30,40,50),
                        labels = c(0,5,10,20,30,40,50)) +
  scale_y_continuous(breaks = latbrks, 
                     labels = paste(latbrks, '°', sep = ''),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = lonbrks, 
                     labels = paste(lonbrks, '°', sep = ''),
                     expand = c(0,0)) + 
  scale_colour_manual(values = c("black","blue"))
  
mp6 <- mp6 + guides(size = guide_legend(ncol=1,byrow=T, title.position = "top")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=10, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.spacing.y = unit(0.6, 'cm'))
mp6
png("results/AUS_datapoints_appendix.png", width = 7, height = 3.5, units = "in", res = 600)
mp6
dev.off()
