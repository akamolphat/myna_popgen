# Create_metadata_file.R --------------------------------------------
#
# This script quickly convert TableS1.2.csv to match the DArT 
# input metadata CSV file.
#
# It just renames the ID column in Table S1.2 to "id" instead
#
# Path to input TableS1.2.csv
metacsv <- "../01_download_data/TableS1.2.csv"
# Path to store the output csv file
outcsv <- "../01_download_data/DArT/metadata_all_samples.csv"

dt <- read.csv(metacsv)
colnames(dt)[1] <- "id"

write.csv(dt, outcsv, row.names = F)
