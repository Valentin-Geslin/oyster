# Load libraries ----------------------------------------------------------
library(car)
library(dplyr)
library(forcats)
library(FSA)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(readr)
library(readxl)
library(stringr)
library(tidyr)
library(tidyverse)
library(tidytable)

# IHC detection & quantification ------------------------------------------
## Import dataset ----------------------------------------------------------
### Import IHC_mask_data ---------------------------------------------------
#oyster <- "/path/to/OYSTER/folder/"
oyster <- "/home/valentin/OYSTER/"

setwd(paste(oyster, "results/IHC_tiles/", sep=""))
IHC_mask_data <- list.files(pattern = "*.csv") %>% map_df(~read_csv(., col_types = cols(.default = "c")))
IHC_mask_data$total_pixels <- as.double(IHC_mask_data$total_pixels)
IHC_mask_data$black_pixels <- as.double(IHC_mask_data$black_pixels)
IHC_mask_data$non_black_pixels <- as.double(IHC_mask_data$non_black_pixels)
IHC_mask_data <- na.omit(IHC_mask_data)

IHC_mask_data_summary <- IHC_mask_data %>%
  select(image, non_black_pixels) %>%
  group_by(image) %>%
  summarise(non_black_pixels = sum(non_black_pixels)) %>%
  ungroup()

# Remove 16045-50-16110141B02401002
IHC_mask_data_summary <- subset(IHC_mask_data_summary, image != "16045-50-16110141B02401002")
IHC_mask_data_summary <- na.omit(IHC_mask_data_summary)

### Import image_metadata ---------------------------------------------------
image_metadata <- read.table(paste(oyster, "results/image_metadata.txt", sep=""), row.names=NULL, quote="\"", comment.char="")
image_metadata <- image_metadata %>% select(-1,-2,-3,-4,-6,-8,-11,-14,-17,-20,-23,-25)
names(image_metadata) <- c("image","","width_pixels","","height_pixels","","resolution","","sizeC","","pixel_width_microns","","pixel_height_microns")
image_metadata <- image_metadata %>% select(-2,-4,-6,-8,-10,-12)

image_metadata$image <- str_replace_all(image_metadata$image, '.ndpi', '')
image_metadata$width_pixels <- str_replace_all(image_metadata$width_pixels, ',', '')
image_metadata$height_pixels <- str_replace_all(image_metadata$height_pixels, ',', '')
image_metadata$resolution <- str_replace_all(image_metadata$resolution, ',', '')
image_metadata$sizeC <- str_replace_all(image_metadata$sizeC, ',', '')
image_metadata$pixel_width_microns <- str_replace_all(image_metadata$pixel_width_microns, ',', '')

image_metadata$width_pixels <- as.numeric(image_metadata$width_pixels)
image_metadata$height_pixels <- as.numeric(image_metadata$height_pixels)
image_metadata$resolution <- as.numeric(image_metadata$resolution)
image_metadata$sizeC <- as.numeric(image_metadata$sizeC)
image_metadata$pixel_width_microns <- as.numeric(image_metadata$pixel_width_microns)
image_metadata$image_area <- image_metadata$width_pixels*image_metadata$height_pixels*image_metadata$pixel_width_microns*image_metadata$pixel_height_microns

# Remove 16045-50-16110141B02401002
image_metadata <- subset(image_metadata, image != "16045-50-16110141B02401002")

### Import tissue_area -----------------------------------------------------
tissue_area <- read_csv(paste(oyster, "results/tissue_area.csv", sep=""))
tissue_area$Image <- str_replace_all(tissue_area$Image, '.ndpi', '')
tissue_area$Image <- str_replace_all(tissue_area$Image, '16045-50-moribund-16110141B02401002', '16045-50-16110141B02401002')
colnames(tissue_area)[colnames(tissue_area) == 'Image'] <- 'image'
colnames(tissue_area)[colnames(tissue_area) == 'Area µm^2'] <- 'tissue_area'
colnames(tissue_area)[colnames(tissue_area) == 'Perimeter µm'] <- 'tissue_perimeter'

# Remove 16045-50-16110141B02401002
tissue_area <- subset(tissue_area, image != "16045-50-16110141B02401002")

### Import qPCR_data -------------------------------------------------------
qPCR_data <- read_excel(paste(oyster, "results/qPCR_results.xlsx", sep=""), 
                        sheet = "Sheet1", col_types = c("text","numeric","text","text","numeric","text","text"))

# Rename images
qPCR_data["image"][qPCR_data["number"] ==  '16045-7'] <- '64853'
qPCR_data["image"][qPCR_data["number"] ==  '16045-12'] <- '64852'
qPCR_data["image"][qPCR_data["number"] ==  '16045-26'] <- '64854'
qPCR_data["image"][qPCR_data["number"] ==  '16045-27'] <- '64856'
qPCR_data["image"][qPCR_data["number"] ==  '16045-33'] <- '64857'
qPCR_data["image"][qPCR_data["number"] ==  '16045-1'] <- '16045-1-101002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-10'] <- '16045-10-16110141B01401002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-11'] <- '16045-11-16110141B01501002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-13'] <- '16045-13-601002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-14'] <- '16045-14-16110141B01601002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-15'] <- '16045-15'
qPCR_data["image"][qPCR_data["number"] ==  '16045-16'] <- '16045-16-701002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-17'] <- '16045-17-801002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-18'] <- '16045-18-901002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-19'] <- '16045-19-1001002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-2'] <- '16045-2-16110141B00901002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-20'] <- '16045-20-1101002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-21'] <- '16045-21-1201002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-22'] <- '16045-22-16110141B01701002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-23'] <- '16045-23-16110141B01801002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-24'] <- '16045-24-16110141B01901002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-25'] <- '16045-25-1301002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-28'] <- '16045-28-16110141B02001002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-29'] <- '16045-29-16110141B02101002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-3'] <- '16045-3-16110141B01001002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-30'] <- '16045-30'
qPCR_data["image"][qPCR_data["number"] ==  '16045-31'] <- '16045-31-16110141B02201002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-32'] <- '16045-32-1701002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-4'] <- '16045-4-201002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-5'] <- '16045-5-16110141B01101002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-6'] <- '16045-6-301002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-8'] <- '16045-8-16110141B01201002'
qPCR_data["image"][qPCR_data["number"] ==  '16045-9'] <- '16045-9-16110141B01301002'
qPCR_data["image"][qPCR_data["number"] ==  '16046-1'] <- '16049-1'
qPCR_data["image"][qPCR_data["number"] ==  '16046-10'] <- "16049-10-16110141B00401002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-2'] <- "16049-2-3101002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-3'] <- "16049-3-3201002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-4'] <- '16049-4'
qPCR_data["image"][qPCR_data["number"] ==  '16046-5'] <- "16049-5-16110141B00201002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-6'] <- "16049-6-3301002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-7'] <- "16049-7-16110141B00301002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-8'] <- "16049-8-3401002"
qPCR_data["image"][qPCR_data["number"] ==  '16046-9'] <- '16049-9'

# Transform quantity_bacteria (log10)
qPCR_data$log10_quantity_bacteria <- log10(qPCR_data$quantity_bacteria)
qPCR_data$log10_quantity_bacteria[is.na(qPCR_data$log10_quantity_bacteria) | qPCR_data$log10_quantity_bacteria=="-Inf"] = 0

# Remove HL
qPCR_data <- subset(qPCR_data, tissu != "HL")

# Add column "treatment"
qPCR_data$treatment[qPCR_data$image == "16045-1-101002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-10-16110141B01401002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-11-16110141B01501002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-13-601002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-14-16110141B01601002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-15"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-17-801002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-18-901002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-19-1001002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-2-16110141B00901002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-20-1101002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-21-1201002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-22-16110141B01701002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-23-16110141B01801002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-24-16110141B01901002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-25-1301002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-28-16110141B02001002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-29-16110141B02101002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-3-16110141B01001002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-30"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-31-16110141B02201002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-32-1701002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-4-201002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-5-16110141B01101002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-50-16110141B02401002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-6-301002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-7-401003"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-16-701002"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "16045-8-16110141B01201002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16045-9-16110141B01301002"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "16049-1"] <- "control-J1"
qPCR_data$treatment[qPCR_data$image == "16049-10-16110141B00401002"] <- "control-J4"
qPCR_data$treatment[qPCR_data$image == "16049-2-3101002"] <- "control-J1"
qPCR_data$treatment[qPCR_data$image == "16049-3-3201002"] <- "control-J1"
qPCR_data$treatment[qPCR_data$image == "16049-4"] <- "control-J1"
qPCR_data$treatment[qPCR_data$image == "16049-5-16110141B00201002"] <- "control-J1"
qPCR_data$treatment[qPCR_data$image == "16049-6-3301002"] <- "control-J4"
qPCR_data$treatment[qPCR_data$image == "16049-7-16110141B00301002"] <- "control-J4"
qPCR_data$treatment[qPCR_data$image == "16049-8-3401002"] <- "control-J4"
qPCR_data$treatment[qPCR_data$image == "16049-9"] <- "control-J4"
qPCR_data$treatment[qPCR_data$image == "64851"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "64852"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "64853"] <- "infected-J1"
qPCR_data$treatment[qPCR_data$image == "64854"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "64856"] <- "infected-J4"
qPCR_data$treatment[qPCR_data$image == "64857"] <- "infected-J4"

# Add column "group" based on log10_quantity_bacteria value
# control = control
# 0 < early_infection < 3
# 3 < advanced_infection < 8
ggplot(qPCR_data, aes(x=log10_quantity_bacteria)) +
  geom_histogram() +
  theme_minimal()

qPCR_data$group[qPCR_data$image == "16045-1-101002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-10-16110141B01401002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-11-16110141B01501002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-13-601002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-14-16110141B01601002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-15"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-17-801002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-18-901002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-19-1001002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-2-16110141B00901002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-20-1101002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-21-1201002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-22-16110141B01701002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-23-16110141B01801002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-24-16110141B01901002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-25-1301002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-28-16110141B02001002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-29-16110141B02101002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-3-16110141B01001002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-30"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-31-16110141B02201002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-32-1701002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-4-201002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-5-16110141B01101002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-50-16110141B02401002"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "16045-6-301002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-7-401003"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-16-701002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-8-16110141B01201002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16045-9-16110141B01301002"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "16049-1"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-10-16110141B00401002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-2-3101002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-3-3201002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-4"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-5-16110141B00201002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-6-3301002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-7-16110141B00301002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-8-3401002"] <- "control"
qPCR_data$group[qPCR_data$image == "16049-9"] <- "control"
qPCR_data$group[qPCR_data$image == "64851"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "64852"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "64853"] <- "early_infection"
qPCR_data$group[qPCR_data$image == "64854"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "64856"] <- "advanced_infection"
qPCR_data$group[qPCR_data$image == "64857"] <- "advanced_infection"

qPCR_data$treatment <- as.factor(qPCR_data$treatment)
qPCR_data$group <- as.factor(qPCR_data$group)

# Remove 16045-50-16110141B02401002
qPCR_data <- subset(qPCR_data, image != "16045-50-16110141B02401002")

# Keep only matching samples from image_metadata
qPCR_data <- subset(qPCR_data, image %in% image_metadata$image)

# Recode image 
# I_J1
qPCR_data$code[qPCR_data$image == "16045-1-101002"] <- "I_J1_1"
qPCR_data$code[qPCR_data$image == "16045-2-16110141B00901002"] <- "I_J1_2"
qPCR_data$code[qPCR_data$image == "16045-3-16110141B01001002"] <- "I_J1_3"
qPCR_data$code[qPCR_data$image == "16045-4-201002"] <- "I_J1_4"
qPCR_data$code[qPCR_data$image == "16045-5-16110141B01101002"] <- "I_J1_5"
qPCR_data$code[qPCR_data$image == "16045-6-301002"] <- "I_J1_6"
qPCR_data$code[qPCR_data$image == "16045-8-16110141B01201002"] <- "I_J1_7"
qPCR_data$code[qPCR_data$image == "16045-9-16110141B01301002"] <- "I_J1_8"
qPCR_data$code[qPCR_data$image == "16045-10-16110141B01401002"] <- "I_J1_9"
qPCR_data$code[qPCR_data$image == "16045-11-16110141B01501002"] <- "I_J1_10"
qPCR_data$code[qPCR_data$image == "16045-13-601002"] <- "I_J1_11"
qPCR_data$code[qPCR_data$image == "16045-14-16110141B01601002"] <- "I_J1_12"
qPCR_data$code[qPCR_data$image == "64852"] <- "I_J1_13"
qPCR_data$code[qPCR_data$image == "64853"] <- "I_J1_14"
qPCR_data$code[qPCR_data$image == "64851"] <- "EI_J1_15"
qPCR_data$code[qPCR_data$image == "16045-7-401003"] <- "I_J1_16"

# I_J4
qPCR_data$code[qPCR_data$image == "16045-16-701002"] <- "I_J4_1"
qPCR_data$code[qPCR_data$image == "16045-17-801002"] <- "I_J4_2"
qPCR_data$code[qPCR_data$image == "16045-18-901002"] <- "I_J4_3"
qPCR_data$code[qPCR_data$image == "16045-19-1001002"] <- "I_J4_4"
qPCR_data$code[qPCR_data$image == "16045-20-1101002"] <- "I_J4_5"
qPCR_data$code[qPCR_data$image == "16045-21-1201002"] <- "I_J4_6"
qPCR_data$code[qPCR_data$image == "16045-22-16110141B01701002"] <- "I_J4_7"
qPCR_data$code[qPCR_data$image == "16045-23-16110141B01801002"] <- "I_J4_8"
qPCR_data$code[qPCR_data$image == "16045-24-16110141B01901002"] <- "I_J4_9"
qPCR_data$code[qPCR_data$image == "16045-25-1301002"] <- "I_J4_10"
qPCR_data$code[qPCR_data$image == "16045-28-16110141B02001002"] <- "I_J4_11"
qPCR_data$code[qPCR_data$image == "16045-29-16110141B02101002"] <- "I_J4_12"
qPCR_data$code[qPCR_data$image == "16045-31-16110141B02201002"] <- "I_J4_13"
qPCR_data$code[qPCR_data$image == "16045-32-1701002"] <- "I_J4_14"
qPCR_data$code[qPCR_data$image == "64854"] <- "I_J4_15"
qPCR_data$code[qPCR_data$image == "64856"] <- "I_J4_16"
qPCR_data$code[qPCR_data$image == "64857"] <- "I_J4_17"

# C_J1
qPCR_data$code[qPCR_data$image == "16049-2-3101002"] <- "C_J1_1"
qPCR_data$code[qPCR_data$image == "16049-3-3201002"] <- "C_J1_2"
qPCR_data$code[qPCR_data$image == "16049-5-16110141B00201002"] <- "C_J1_3"

# C_J4
qPCR_data$code[qPCR_data$image == "16049-6-3301002"] <- "C_J4_1"
qPCR_data$code[qPCR_data$image == "16049-7-16110141B00301002"] <- "C_J4_2"
qPCR_data$code[qPCR_data$image == "16049-8-3401002"] <- "C_J4_3"
qPCR_data$code[qPCR_data$image == "16049-10-16110141B00401002"] <- "C_J4_4"

qPCR_data$code <- as.factor(qPCR_data$code)

# Summarise qPCR_data
qPCR_mean <- qPCR_data %>%
  group_by(image) %>%
  summarise(mean_quantity_bacteria = mean(quantity_bacteria, na.rm = TRUE),
            sd_quantity_bacteria = sd(quantity_bacteria, na.rm = TRUE))

# log10 transform mean_quantity_bacteria
qPCR_mean$log10_mean_quantity_bacteria <- log10(qPCR_mean$mean_quantity_bacteria)
qPCR_mean$log10_mean_quantity_bacteria[qPCR_mean$log10_mean_quantity_bacteria<0] = 0

## Recode image name -------------------------------------------------------
# I_J1
qPCR_mean$code[qPCR_mean$image == "16045-1-101002"] <- "I_J1_1"
qPCR_mean$code[qPCR_mean$image == "16045-2-16110141B00901002"] <- "I_J1_2"
qPCR_mean$code[qPCR_mean$image == "16045-3-16110141B01001002"] <- "I_J1_3"
qPCR_mean$code[qPCR_mean$image == "16045-4-201002"] <- "I_J1_4"
qPCR_mean$code[qPCR_mean$image == "16045-5-16110141B01101002"] <- "I_J1_5"
qPCR_mean$code[qPCR_mean$image == "16045-6-301002"] <- "I_J1_6"
qPCR_mean$code[qPCR_mean$image == "16045-7-401003"] <- "I_J1_7"
qPCR_mean$code[qPCR_mean$image == "16045-8-16110141B01201002"] <- "I_J1_8"
qPCR_mean$code[qPCR_mean$image == "16045-9-16110141B01301002"] <- "I_J1_9"
qPCR_mean$code[qPCR_mean$image == "16045-10-16110141B01401002"] <- "I_J1_10"
qPCR_mean$code[qPCR_mean$image == "16045-11-16110141B01501002"] <- "I_J1_11"
qPCR_mean$code[qPCR_mean$image == "16045-13-601002"] <- "I_J1_12"
qPCR_mean$code[qPCR_mean$image == "16045-14-16110141B01601002"] <- "I_J1_13"
qPCR_mean$code[qPCR_mean$image == "64852"] <- "I_J1_14"
qPCR_mean$code[qPCR_mean$image == "64853"] <- "I_J1_15"
qPCR_mean$code[qPCR_mean$image == "64851"] <- "I_J1_16"

# I_J4
qPCR_mean$code[qPCR_mean$image == "16045-16-701002"] <- "I_J4_1"
qPCR_mean$code[qPCR_mean$image == "16045-17-801002"] <- "I_J4_2"
qPCR_mean$code[qPCR_mean$image == "16045-18-901002"] <- "I_J4_3"
qPCR_mean$code[qPCR_mean$image == "16045-19-1001002"] <- "I_J4_4"
qPCR_mean$code[qPCR_mean$image == "16045-20-1101002"] <- "I_J4_5"
qPCR_mean$code[qPCR_mean$image == "16045-21-1201002"] <- "I_J4_6"
qPCR_mean$code[qPCR_mean$image == "16045-22-16110141B01701002"] <- "I_J4_7"
qPCR_mean$code[qPCR_mean$image == "16045-23-16110141B01801002"] <- "I_J4_8"
qPCR_mean$code[qPCR_mean$image == "16045-24-16110141B01901002"] <- "I_J4_9"
qPCR_mean$code[qPCR_mean$image == "16045-25-1301002"] <- "I_J4_10"
qPCR_mean$code[qPCR_mean$image == "16045-28-16110141B02001002"] <- "I_J4_11"
qPCR_mean$code[qPCR_mean$image == "16045-29-16110141B02101002"] <- "I_J4_12"
qPCR_mean$code[qPCR_mean$image == "16045-31-16110141B02201002"] <- "I_J4_13"
qPCR_mean$code[qPCR_mean$image == "16045-32-1701002"] <- "I_J4_14"
qPCR_mean$code[qPCR_mean$image == "64854"] <- "I_J4_15"
qPCR_mean$code[qPCR_mean$image == "64856"] <- "I_J4_16"
qPCR_mean$code[qPCR_mean$image == "64857"] <- "I_J4_17"

# C_J1
qPCR_mean$code[qPCR_mean$image == "16049-2-3101002"] <- "C_J1_1"
qPCR_mean$code[qPCR_mean$image == "16049-3-3201002"] <- "C_J1_2"
qPCR_mean$code[qPCR_mean$image == "16049-5-16110141B00201002"] <- "C_J1_3"

# C_J4
qPCR_mean$code[qPCR_mean$image == "16049-6-3301002"] <- "C_J4_1"
qPCR_mean$code[qPCR_mean$image == "16049-7-16110141B00301002"] <- "C_J4_2"
qPCR_mean$code[qPCR_mean$image == "16049-8-3401002"] <- "C_J4_3"
qPCR_mean$code[qPCR_mean$image == "16049-10-16110141B00401002"] <- "C_J4_4"

qPCR_mean$code <- as.factor(qPCR_mean$code)
qPCR_mean$image <- as.factor(qPCR_mean$image)

qPCR_mean$group[qPCR_mean$image == "16045-1-101002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-10-16110141B01401002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-11-16110141B01501002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-13-601002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-14-16110141B01601002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-15"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-17-801002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-18-901002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-19-1001002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-2-16110141B00901002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-20-1101002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-21-1201002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-22-16110141B01701002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-23-16110141B01801002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-24-16110141B01901002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-25-1301002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-28-16110141B02001002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-29-16110141B02101002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-3-16110141B01001002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-30"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-31-16110141B02201002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-32-1701002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-4-201002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-5-16110141B01101002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-50-16110141B02401002"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "16045-6-301002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-7-401003"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-16-701002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-8-16110141B01201002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16045-9-16110141B01301002"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "16049-1"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-10-16110141B00401002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-2-3101002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-3-3201002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-4"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-5-16110141B00201002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-6-3301002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-7-16110141B00301002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-8-3401002"] <- "control"
qPCR_mean$group[qPCR_mean$image == "16049-9"] <- "control"
qPCR_mean$group[qPCR_mean$image == "64851"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "64852"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "64853"] <- "early_infection"
qPCR_mean$group[qPCR_mean$image == "64854"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "64856"] <- "advanced_infection"
qPCR_mean$group[qPCR_mean$image == "64857"] <- "advanced_infection"

qPCR_mean$group <- as.factor(qPCR_mean$group)

## Merge dataset -----------------------------------------------------------
IHC_mask_data_summary <- subset(IHC_mask_data_summary, image %in% qPCR_mean$image)
IHC_mask_data_summary$image <- as.factor(IHC_mask_data_summary$image)

image_metadata <- subset(image_metadata, image %in% qPCR_mean$image)
image_metadata$image <- as.factor(image_metadata$image)

tissue_area <- subset(tissue_area, image %in% qPCR_mean$image)
tissue_area$image <- as.factor(tissue_area$image)

image_data <- merge(IHC_mask_data_summary, tissue_area, by = "image")
image_data <- merge(image_data, image_metadata, by = "image")
image_data <- merge(image_data, qPCR_mean, by = "image")

image_data$image <- as.factor(image_data$image)

# Add column "treatment"
image_data$treatment[image_data$image == "16045-1-101002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-10-16110141B01401002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-11-16110141B01501002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-13-601002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-14-16110141B01601002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-16-701002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-17-801002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-18-901002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-19-1001002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-2-16110141B00901002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-20-1101002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-21-1201002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-22-16110141B01701002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-23-16110141B01801002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-24-16110141B01901002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-25-1301002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-28-16110141B02001002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-29-16110141B02101002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-3-16110141B01001002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-31-16110141B02201002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-32-1701002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-4-201002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-5-16110141B01101002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-50-16110141B02401002"] <- "infected-J4"
image_data$treatment[image_data$image == "16045-6-301002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-7-401003"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-8-16110141B01201002"] <- "infected-J1"
image_data$treatment[image_data$image == "16045-9-16110141B01301002"] <- "infected-J1"
image_data$treatment[image_data$image == "16049-10-16110141B00401002"] <- "control-J4"
image_data$treatment[image_data$image == "16049-2-3101002"] <- "control-J1"
image_data$treatment[image_data$image == "16049-3-3201002"] <- "control-J1"
image_data$treatment[image_data$image == "16049-5-16110141B00201002"] <- "control-J1"
image_data$treatment[image_data$image == "16049-6-3301002"] <- "control-J4"
image_data$treatment[image_data$image == "16049-7-16110141B00301002"] <- "control-J4"
image_data$treatment[image_data$image == "16049-8-3401002"] <- "control-J4"
image_data$treatment[image_data$image == "64851"] <- "infected-J1"
image_data$treatment[image_data$image == "64852"] <- "infected-J1"
image_data$treatment[image_data$image == "64853"] <- "infected-J1"
image_data$treatment[image_data$image == "64854"] <- "infected-J4"
image_data$treatment[image_data$image == "64856"] <- "infected-J4"
image_data$treatment[image_data$image == "64857"] <- "infected-J4"

# Add column "group"
image_data$group[image_data$image == "16045-1-101002"] <- "early_infection"
image_data$group[image_data$image == "16045-10-16110141B01401002"] <- "early_infection"
image_data$group[image_data$image == "16045-11-16110141B01501002"] <- "early_infection"
image_data$group[image_data$image == "16045-13-601002"] <- "early_infection"
image_data$group[image_data$image == "16045-14-16110141B01601002"] <- "early_infection"
image_data$group[image_data$image == "16045-16-701002"] <- "early_infection"
image_data$group[image_data$image == "16045-17-801002"] <- "early_infection"
image_data$group[image_data$image == "16045-18-901002"] <- "early_infection"
image_data$group[image_data$image == "16045-19-1001002"] <- "early_infection"
image_data$group[image_data$image == "16045-2-16110141B00901002"] <- "early_infection"
image_data$group[image_data$image == "16045-20-1101002"] <- "early_infection"
image_data$group[image_data$image == "16045-21-1201002"] <- "early_infection"
image_data$group[image_data$image == "16045-22-16110141B01701002"] <- "early_infection"
image_data$group[image_data$image == "16045-23-16110141B01801002"] <- "early_infection"
image_data$group[image_data$image == "16045-24-16110141B01901002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-25-1301002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-28-16110141B02001002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-29-16110141B02101002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-3-16110141B01001002"] <- "early_infection"
image_data$group[image_data$image == "16045-31-16110141B02201002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-32-1701002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-4-201002"] <- "early_infection"
image_data$group[image_data$image == "16045-5-16110141B01101002"] <- "early_infection"
image_data$group[image_data$image == "16045-50-16110141B02401002"] <- "advanced_infection"
image_data$group[image_data$image == "16045-6-301002"] <- "early_infection"
image_data$group[image_data$image == "16045-7-401003"] <- "early_infection"
image_data$group[image_data$image == "16045-8-16110141B01201002"] <- "early_infection"
image_data$group[image_data$image == "16045-9-16110141B01301002"] <- "early_infection"
image_data$group[image_data$image == "16049-10-16110141B00401002"] <- "control"
image_data$group[image_data$image == "16049-2-3101002"] <- "control"
image_data$group[image_data$image == "16049-3-3201002"] <- "control"
image_data$group[image_data$image == "16049-5-16110141B00201002"] <- "control"
image_data$group[image_data$image == "16049-6-3301002"] <- "control"
image_data$group[image_data$image == "16049-7-16110141B00301002"] <- "control"
image_data$group[image_data$image == "16049-8-3401002"] <- "control"
image_data$group[image_data$image == "64851"] <- "early_infection"
image_data$group[image_data$image == "64852"] <- "early_infection"
image_data$group[image_data$image == "64853"] <- "early_infection"
image_data$group[image_data$image == "64854"] <- "advanced_infection"
image_data$group[image_data$image == "64856"] <- "advanced_infection"
image_data$group[image_data$image == "64857"] <- "advanced_infection"

image_data$treatment <- as.factor(image_data$treatment)
image_data$group <- as.factor(image_data$group)

# Calculate prop_IHC_tissue
image_data$image <- as.factor(image_data$image)
image_data$IHC_area <- image_data$non_black_pixels*image_data$pixel_width_microns*image_data$pixel_height_microns

# Correct IHC area based on "control"
max_control_IHC_area <- max(image_data$IHC_area[image_data$group == "control"])
image_data$IHC_area_corrected <- image_data$IHC_area-max_control_IHC_area
image_data$IHC_area_corrected[image_data$IHC_area_corrected < 0] <- 0

image_data$prop_tissue_image <- (image_data$tissue_area*100)/image_data$image_area
image_data$prop_IHC_tissue <- (image_data$IHC_area*100)/image_data$tissue_area
image_data$prop_IHC_tissue_corrected <- (image_data$IHC_area_corrected*100)/image_data$tissue_area

#rm(IHC_mask_data, IHC_mask_data_summary, image_metadata, tissue_area, qPCR_mean)

## Data visualization & analysis ------------------------------------------------------
output_path <- paste(oyster, "output/", sep="")

### qPCR_data --------------------------------------------------------------
ggplot(qPCR_data, aes(x=log10_quantity_bacteria, fill=group)) +
  geom_histogram(binwidth = 1, color='black', position="dodge") +
  theme_minimal()
#ggsave("distribution_quantity_bacteria_image_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

f1 <- ggplot(qPCR_mean, aes(x=log10_mean_quantity_bacteria, fill=group)) +
  geom_histogram(binwidth = 1, color='black', position="dodge") +
  labs(y = "Count", x ="Mean quantity of bacteria cells per 25ng of total DNA (log10)") +
  theme_minimal() +
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.title.x = element_text( size=15),
        axis.title.y = element_text(size=15))
f1 + scale_fill_discrete(name = "Infection group", labels = c("Advanced infection","Control","Early infection"))
ggsave("distribution_quantity_bacteria_image_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

(p1 <- ggplot(qPCR_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                             "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                             "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                             "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                             y=log10_quantity_bacteria, fill=group)) +
  geom_boxplot() +
  labs(x = "Images", y ="Quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme_minimal() +
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.title.x = element_text( size=15),
        axis.title.y = element_text(size=15)))
p1 + scale_fill_discrete(name = "Infection group", labels = c("Advanced infection","Control","Early infection"))
ggsave("quantity_bacteria_image_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(qPCR_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                             "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                             "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                             "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                      y=log10_quantity_bacteria, fill=treatment)) +
  geom_boxplot() +
  labs(x = "Images", y ="Quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("quantity_bacteria_image_treatment.png", path = output_path, width = 29.7, height = 21, units = "cm")

qPCR_data_summary <- qPCR_data %>%
  group_by(image) %>%
  summarize(mean_qPCR = mean(log10_quantity_bacteria, na.rm = TRUE),
            median_qPCR = median(log10_quantity_bacteria, na.rm = TRUE),
            sd_qPCR = sd(log10_quantity_bacteria, na.rm = TRUE))

write_csv(qPCR_data_summary, paste(oyster, "results/qPCR_data_summary.csv", sep=""))

summary(qPCR_data$log10_quantity_bacteria)
summary(qPCR_data$log10_quantity_bacteria[qPCR_data$group == "advanced_infection"])
summary(qPCR_data$log10_quantity_bacteria[qPCR_data$group == "early_infection"])
summary(qPCR_data$log10_quantity_bacteria[qPCR_data$group == "control"])

kruskal.test(log10_quantity_bacteria ~ group, data = qPCR_data)
dunnTest(log10_quantity_bacteria ~ group, data = qPCR_data)

### prop_IHC_tissue_corrected --------------------------------------------------------
(p2 <- ggplot(image_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                                     "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                                     "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                                     "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                              y=prop_IHC_tissue_corrected, fill=group)) +
   geom_bar(stat="identity", position= "dodge") +
   labs(x = "Images", y ="Proportion of tissue infected by IHC stained bacteria (%)") +
   theme_minimal() +
   theme(legend.title = element_text(size = 15), 
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle=45, vjust=1, hjust=1),
         axis.title.x = element_text( size=15),
         axis.title.y = element_text(size=15)))
p2 + scale_fill_discrete(name = "Infection group", labels = c("Advanced infection","Control","Early infection"))
ggsave("prop_IHC_area_tissue_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=reorder(code, prop_IHC_tissue_corrected), y=prop_IHC_tissue_corrected, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Proportion of tissue infected by IHC stained bacteria (%)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("prop_IHC_area_tissue_group_reordered.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                              "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                              "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                              "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                       y=prop_IHC_tissue_corrected, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Proportion of tissue infected by IHC stained bacteria (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("prop_IHC_area_tissue_treatment.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=reorder(code, prop_IHC_tissue_corrected), y=prop_IHC_tissue_corrected, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Proportion of tissue infected by IHC stained bacteria (%)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# Merge p1 & p2
ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE) 
ggsave("quantity_bacteria_prop_IHC_area_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

# summary prop_IHC_tissue_corrected
summary(image_data$prop_IHC_tissue_corrected)
summary(image_data$prop_IHC_tissue_corrected[image_data$group == "control"])
summary(image_data$prop_IHC_tissue_corrected[image_data$group == "early_infection"])
summary(image_data$prop_IHC_tissue_corrected[image_data$group == "advanced_infection"])

# test for normality
shapiro.test(image_data$prop_IHC_tissue_corrected)
qqnorm(image_data$prop_IHC_tissue_corrected)
qqline(image_data$prop_IHC_tissue_corrected)
hist(image_data$prop_IHC_tissue_corrected)

# wilcox.test is used to determine whether the median of the sample is equal to a known standard value (mu)
wilcox.test(x = image_data$prop_IHC_tissue_corrected[image_data$group== "control"], mu = 0, alternative = "greater")
wilcox.test(x = image_data$prop_IHC_tissue_corrected[image_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = image_data$prop_IHC_tissue_corrected[image_data$group== "advanced_infection"], mu = 0, alternative = "greater")

kruskal.test(prop_IHC_tissue_corrected ~ group, data = image_data) 
dunnTest(prop_IHC_tissue_corrected ~ group, data = image_data)

### IHC_area_corrected --------------------------------------------------------
ggplot(image_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                              "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                              "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                              "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                       y=IHC_area_corrected, fill=group)) +
  geom_bar(stat="identity", position= "dodge") +
  labs(x = "Images", y ="Area of tissue infected by IHC stained bacteria (?m^2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("IHC_area_corrected_tissue_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=reorder(code, IHC_area_corrected), y=IHC_area, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Area of tissue infected by IHC stained bacteria (?m^2)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("IHC_area_corrected_tissue_group_reordered.png", path = output_path, width = 29.7, height = 21, units = "cm")

# summary(IHC_area_corrected)
summary(image_data$IHC_area_corrected)
summary(image_data$IHC_area_corrected[image_data$group == "advanced_infection"])
summary(image_data$IHC_area_corrected[image_data$group == "early_infection"])
summary(image_data$IHC_area_corrected[image_data$group == "control"])

# test for normality
shapiro.test(image_data$IHC_area_corrected)
qqnorm(image_data$IHC_area_corrected)
qqline(image_data$IHC_area_corrected)
hist(image_data$IHC_area_corrected)

# wilcox.test is used to determine whether the median of the sample is equal to a known standard value (mu)
wilcox.test(x = image_data$IHC_area_corrected[image_data$group== "advanced_infection"], mu = 0, alternative = "greater")
wilcox.test(x = image_data$IHC_area_corrected[image_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = image_data$IHC_area_corrected[image_data$group== "control"], mu = 0, alternative = "greater")

kruskal.test(IHC_area_corrected ~ group, data = image_data) 
dunnTest(IHC_area_corrected ~ group, data = image_data)

### log10_mean_quantity_bacteria -------------------------------------------
ggplot(image_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                              "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                              "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                              "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                       y=log10_mean_quantity_bacteria, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Mean quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("mean_quantity_bacteria_image_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=reorder(code, log10_mean_quantity_bacteria), y=log10_mean_quantity_bacteria, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Mean quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("mean_quantity_bacteria_image_group_reordered.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                              "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                              "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                              "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                       y=log10_mean_quantity_bacteria, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Mean quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("mean_quantity_bacteria_image_treatment.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=reorder(code, log10_mean_quantity_bacteria), y=log10_mean_quantity_bacteria, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Mean quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#summary(log10_mean_quantity_bacteria)
summary(image_data$log10_mean_quantity_bacteria)
summary(image_data$log10_mean_quantity_bacteria[image_data$group == "advanced_infection"])
summary(image_data$log10_mean_quantity_bacteria[image_data$group == "early_infection"])
summary(image_data$log10_mean_quantity_bacteria[image_data$group == "control"])

# test for normality
shapiro.test(image_data$log10_mean_quantity_bacteria)
qqnorm(image_data$log10_mean_quantity_bacteria)
qqline(image_data$log10_mean_quantity_bacteria)
hist(image_data$log10_mean_quantity_bacteria)

# wilcox.test is used to determine whether the median of the sample is equal to a known standard value (mu)
wilcox.test(x = image_data$log10_mean_quantity_bacteria[image_data$group== "advanced_infection"], mu = 0, alternative = "greater")
wilcox.test(x = image_data$log10_mean_quantity_bacteria[image_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = image_data$log10_mean_quantity_bacteria[image_data$group== "control"], mu = 0, alternative = "greater")

kruskal.test(log10_mean_quantity_bacteria ~ group, data = image_data)
dunnTest(log10_mean_quantity_bacteria ~ group, data = image_data)
prop_IHC_tissue_corrected

### log10_mean_quantity_bacteria ~ IHC_area_corrected -----------------------
ggplot(image_data, aes(x=IHC_area_corrected, y=log10_mean_quantity_bacteria)) +
  geom_point(aes(colour = group)) +
  labs(x="IHC positive tissue area (?m^2))", y="Mean quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme_minimal() 

# Correlation test
cor.test(image_data$IHC_area_corrected, image_data$log10_mean_quantity_bacteria, method = "kendall")


### log10_mean_quantity_bacteria ~ prop_IHC_tissue_corrected --------------------------
(g1 <- ggplot(image_data, aes(x=prop_IHC_tissue_corrected, y=log10_mean_quantity_bacteria)) +
  geom_point(aes(colour = group)) +
  labs(x="Proportion of tissue infected by IHC stained bacteria (%)", y="Mean quantity of bacteria cells per 25ng of total DNA (log10)") +
  theme_minimal() +
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.title.x = element_text( size=15),
        axis.title.y = element_text(size=15)))
g1 + scale_colour_discrete(name = "Infection group", labels = c("Advanced infection","Control","Early infection"))
ggsave("log10_mean_quantity_bacteria_prop_IHC_tissue_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(image_data, aes(x=prop_IHC_tissue_corrected, y=log10_mean_quantity_bacteria)) +
  geom_point(aes(colour = treatment)) + 
  labs(x="Proportion of tissue infected by IHC stained bacteria (%)", y="Mean quantity of bacteria cells per 25ng of total DNA for each sample (log10)") +
  theme_minimal() #+ geom_smooth(method = "loess", formula = y ~ log10(x))
ggsave("log10_mean_quantity_bacteria_prop_IHC_tissue__treatment.png", path = output_path, width = 29.7, height = 21, units = "cm")

# Correlation test
cor.test(image_data$prop_IHC_tissue_corrected, image_data$log10_mean_quantity_bacteria, method = "kendall")

# Logarithmic regression model (model <- lm(y ~ log(x)))
model1 <- lm(image_data$log10_mean_quantity_bacteria ~ log(image_data$prop_IHC_tissue_corrected))
summary(model1)
coef(model1)

res1 <- resid(model1)
hist(res1) # check normality of the residuals 
plot(fitted(model1), res1) # check variance of the residuals
abline(0,0)

# Extrapolate quantity of bacteria from IHC_area with logarithmic regression model
image_data$extrapol_qt_bacteria_1 <- predict(model1, newdata=image_data)
image_data$extrapol_qt_bacteria_1[image_data$extrapol_qt_bacteria_1 < 0] <- 0

ggplot(image_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                              "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                              "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                              "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                       y=extrapol_qt_bacteria_1, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Extrapolated quantity of bacteria cells per animal (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

ggplot(image_data, aes(x=prop_IHC_tissue, y=extrapol_qt_bacteria_1)) +
  geom_point(aes(colour = group)) +
  labs(x="Proportion of tissue infected by IHC stained bacteria (%)", y="Extrapolated mean quantity of bacteria cells per animal (log10)") +
  theme_minimal() #+ geom_smooth(method = "lm", formula = y ~ log10(x))

ggplot(data=image_data) + 
  geom_point(aes(x=prop_IHC_tissue,y=log10_mean_quantity_bacteria), colour = "black") + # original data
  geom_point(aes(x=prop_IHC_tissue,y=extrapol_qt_bacteria_1), colour = "blue")

# Colocalization ----------------------------------------------------------
## Import dataset ----------------------------------------------------------
setwd(paste(oyster, "results/colocalization", sep=""))
folder_path <- paste(oyster, "results/colocalization", sep="")
file_list <- list.files(folder_path, full.names = TRUE)
data_list <- list()

for (file in file_list) {
  data <- read.delim(file, stringsAsFactors = FALSE)
  colnames(data)[1] <- 'Coloc'
  data$Coloc <- str_replace(data$Coloc, "Image B:", "Image B=")
  data[c('Coef', 'Value')] <- str_split_fixed(data$Coloc, '=', 2)
  data <- data %>% select(-Coloc)
  data <- data[-c(2), ]
  data$Value <- str_replace(data$Value, ".png", "")
  data$Coef <- str_replace(data$Coef, "Image B", "image")
  data <- data %>% pivot_wider(names_from = "Coef", values_from = "Value")
  data_list[[file]] <- data
}

coloc <- bind_rows(data_list)
coloc$image <- str_replace(coloc$image, " ", "")
coloc <- na.omit(coloc)
coloc$image <- as.factor(coloc$image)
coloc$r <- as.numeric(coloc$r)
coloc <- subset(coloc, image %in% image_data$image)

## Merge data set ----------------------------------------------------------
coloc_data <- merge(coloc, image_data)
#rm(coloc, data, data_list, file, file_list, folder_path)

## Data visualisation ------------------------------------------------------
output_path <- paste(oyster, "output/colocalization/", sep="")

(p3 <- ggplot(coloc_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                                     "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                                     "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                                     "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                              y=r, color=group)) +
  geom_point(size=3) +
  labs(x = "Images", y ="Colocalization index (Pearson's coefficient value)") +
  theme_minimal() +
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.title.x = element_text( size=15),
        axis.title.y = element_text(size=15)))
p3 + scale_color_discrete(name = "Infection group", labels = c("Advanced infection","Control","Early infection"))
ggsave("pearson_coef_r_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(coloc_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                              "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                              "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                              "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                       y=r, color=treatment)) +
  geom_point(size=3) +
  labs(x = "Images", y ="Colocalization index (Pearson's coefficient value)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme_minimal()
ggsave("pearson_coef_r_treatment.png", path = output_path, width = 29.7, height = 21, units = "cm")

'''
### WE ARE NOT TESTING RELATIONSHIP BTW Pearsons Coefficient & prop_IHC_tissue_corrected
(p4 <- ggplot(coloc_data, aes(y=r, x=prop_IHC_tissue_corrected)) +
  geom_point(aes(colour = group)) +
  labs(y="Pearson Coefficient", x="Proportion of tissue infected by IHC stained bacteria (%)") +
  theme_minimal())
ggsave("pearson_coef_r_prop_IHC_area_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

# Merge p3 & p4
ggarrange(p3, p4, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("pearson_coef_prop_IHC_area_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

### WE ARE NOT TESTING RELATIONSHIP BTW Pearsons Coefficient & log10_mean_quantity_bacteria
(p5 <- ggplot(coloc_data, aes(y=r, x=log10_mean_quantity_bacteria)) +
    geom_point(aes(colour = group)) +
    labs(y="Pearson Coefficient", x="log10_mean_quantity_bacteria") +
    theme_minimal())
ggsave("pearson_coef_r_log10_mean_quantity_bacteria_group.png", path = output_path, width = 29.7, height = 21, units = "cm")
'''

## Data analysis -----------------------------------------------------------
summary(coloc_data$r)
summary(coloc_data$r[coloc_data$group == "control"])
summary(coloc_data$r[coloc_data$group == "early_infection"])
summary(coloc_data$r[coloc_data$group == "advanced_infection"])

coloc_data %>%
  group_by(group) %>%
  summarise(mean = mean(r),
            median = median(r),
            sd = sd(r))

shapiro.test(coloc_data$r)
hist(coloc_data$r)
qqnorm(coloc_data$r)
qqline(coloc_data$r)

wilcox.test(x = coloc_data$r[coloc_data$group== "advanced_infection"], mu = 0, alternative = "greater")
wilcox.test(x = coloc_data$r[coloc_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = coloc_data$r[coloc_data$group== "control"], mu = 0, alternative = "greater")

kruskal.test(r ~ group, data = coloc_data)
dunnTest(r ~ group, data = coloc_data, method="bonferroni")

'''
### WE ARE NOT TESTING RELATIONSHIP BTW Pearsons Coefficient & prop_IHC_tissue_corrected
cor.test(coloc_data$r, coloc_data$prop_IHC_tissue_corrected, method = "kendall")

# Logarithmic regression model (model <- lm(y ~ log(x)))
model2 <- lm(coloc_data$r ~ log(coloc_data$prop_IHC_tissue_corrected))
summary(model2)

res <- resid(model2)
hist(res) # check normality of the residuals 
plot(fitted(model2), res) # check variance of the residuals
abline(0,0)
'''

# Dispersion ---------------------------------------------------
## Import dataset ----------------------------------------------------------
setwd(paste(oyster, "results/dispersion/", sep=""))
folder_path <- paste(oyster, "results/dispersion", sep="")
file_list <- list.files(folder_path, full.names = TRUE)
data_list <- list()

for (file in file_list) {
  data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  data$Value <- parse_number(data$V1)
  data[c('Factor', 'Value_1')] <- str_split_fixed(data$V1, '=', 2)
  data <- data[-c(10,12,16,17,18,19,20),]
  data$test <- str_extract(data$Value_1, "[a-z]+")
  data$image <- file
  data_list[[file]] <- data
}

dispersion <- bind_rows(data_list)
dispersion <- dispersion[,-c(1, 2)]
dispersion$image <- str_replace(dispersion$image, folder_path, "")
dispersion$image <- str_replace(dispersion$image, "/", "")
dispersion$image <- str_replace(dispersion$image, ".txt", "")
dispersion$Factor <- str_replace(dispersion$Factor, "Theoretical random nearest neighbor distance ", "theoretical_random_NN_distance")
dispersion$Factor <- str_replace(dispersion$Factor, "Measured average nearest neighbor distance ", "measured_mean_NN_distance")
dispersion$Factor <- str_replace(dispersion$Factor, "Measured median nearest neighbor distance ", "measured_median_NN_distance")
dispersion$Factor <- str_replace(dispersion$Factor, "confidence interval ", "confidence_interval_percentage")
dispersion$Factor <- str_replace(dispersion$Factor, "Sample size n ", "sample_size_n")
dispersion$test <- str_replace(dispersion$test, "elch", "Welch_test")
dispersion$test <- str_replace(dispersion$test, "tudent", "Student_test")
dispersion$test <- str_replace(dispersion$test, "nfinity", "Student_test")
dispersion$Value <- str_replace(dispersion$Value, "(Student's t-test)", "")
dispersion$Value <- str_replace(dispersion$Value, "(Welch's t-test)", "")
dispersion$Value <- gsub("[()]", "", dispersion$Value) 
dispersion$image <- str_replace(dispersion$image, " ", "")
dispersion$Value <- str_replace(dispersion$Value, "%", "")
dispersion <- subset(dispersion, !is.na(Value))
dispersion <- dispersion[,-c(2)]

rm(data, data_list, file, file_list, folder_path)

dispersion_summary <- dispersion %>%
  pivot_wider(
    names_from = Factor, 
    values_from = Value)

colnames(dispersion_summary)[colnames(dispersion_summary) == "Variance "] <- 'variance'
colnames(dispersion_summary)[colnames(dispersion_summary) == "StdDev "] <- 'StdDev'
colnames(dispersion_summary)[colnames(dispersion_summary) == "Variance (mean) "] <- 'variance_mean_NN_distance'
colnames(dispersion_summary)[colnames(dispersion_summary) == "StdDev (mean) "] <- 'StdDev_mean_NN_distance'
colnames(dispersion_summary)[colnames(dispersion_summary) == "Variance (median) "] <- 'variance_median_NN_distance'
colnames(dispersion_summary)[colnames(dispersion_summary) == "StdDev (median) "] <- 'StdDev_median_NN_distance'
colnames(dispersion_summary)[colnames(dispersion_summary) == "t "] <- 't'
colnames(dispersion_summary)[colnames(dispersion_summary) == "confidence_interval_%"] <- 'confidence_interval_percentage'
colnames(dispersion_summary)[colnames(dispersion_summary) == "critical t-value "] <- 'critical_t_value'

dispersion_summary <- dispersion_summary %>% group_by(image) %>% summarise(across(everything(), ~ na.omit(unique(.))))

dispersion_summary$theoretical_random_NN_distance <- as.numeric(dispersion_summary$theoretical_random_NN_distance)
dispersion_summary$variance <- as.numeric(dispersion_summary$variance)
dispersion_summary$StdDev <- as.numeric(dispersion_summary$StdDev)
dispersion_summary$measured_mean_NN_distance <- as.numeric(dispersion_summary$measured_mean_NN_distance)
dispersion_summary$variance_mean_NN_distance <- as.numeric(dispersion_summary$variance_mean_NN_distance)
dispersion_summary$StdDev_mean_NN_distance <- as.numeric(dispersion_summary$StdDev_mean_NN_distance)
dispersion_summary$measured_median_NN_distance <- as.numeric(dispersion_summary$measured_median_NN_distance)
dispersion_summary$variance_median_NN_distance <- as.numeric(dispersion_summary$variance_median_NN_distance)
dispersion_summary$StdDev_median_NN_distance <- as.numeric(dispersion_summary$StdDev_median_NN_distance)
dispersion_summary$sample_size_n <- as.numeric(dispersion_summary$sample_size_n)
dispersion_summary$t <- as.numeric(dispersion_summary$t)
dispersion_summary$t[dispersion_summary$t=="Inf"] = NA
dispersion_summary$critical_t_value <- as.numeric(dispersion_summary$critical_t_value)
dispersion_summary$confidence_interval_percentage <- as.numeric(dispersion_summary$confidence_interval_percentage)

# Keep only matching samples from image_data (05_Data_analysis.R)
dispersion_summary <- subset(dispersion_summary, image %in% image_data$image)

## Merge dataset -----------------------------------------------------------
dispersion_data <- merge(dispersion_summary, image_data)
#rm(dispersion, dispersion_summary)

dispersion_data$image <- as.factor(dispersion_data$image)
dispersion_data$code <- as.character(dispersion_data$code)

## Data_visualization & analysis ------------------------------------------------------
output_path <- paste(oyster, "output/dispersion/", sep="")

### StdDev_measured_mean_NN_distance ----------------------------------------
ggplot(dispersion_data, aes(x=StdDev_mean_NN_distance, fill=group)) +
  geom_histogram(binwidth = 500, color='black', position="dodge") +
  labs(x = "StdDev Mean nearest neighbor distance (\u00b5m)", y ="Count") +
  theme_minimal()

# "Surface area concentration (",mu, m^2, "/", m^3,")"
summary(dispersion_data$StdDev_mean_NN_distance)
summary(dispersion_data$StdDev_mean_NN_distance[dispersion_data$group == "control"])
summary(dispersion_data$StdDev_mean_NN_distance[dispersion_data$group == "early_infection"])
summary(dispersion_data$StdDev_mean_NN_distance[dispersion_data$group == "advanced_infection"])

shapiro.test(dispersion_data$StdDev_mean_NN_distance) # p-value < 0.05 <- not normally distributed
hist(dispersion_data$StdDev_mean_NN_distance)
qqnorm(dispersion_data$StdDev_mean_NN_distance)
qqline(dispersion_data$StdDev_mean_NN_distance)

wilcox.test(x = dispersion_data$StdDev_mean_NN_distance[dispersion_data$group== "advanced_infection"], mu = 0, alternative = "greater")
wilcox.test(x = dispersion_data$StdDev_mean_NN_distance[dispersion_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = dispersion_data$StdDev_mean_NN_distance[dispersion_data$group== "control"], mu = 0, alternative = "greater")

kruskal.test(StdDev_mean_NN_distance ~ group, data = dispersion_data)
dunnTest(StdDev_mean_NN_distance ~ group, data = dispersion_data, method="bonferroni")

cor.test(dispersion_data$StdDev_mean_NN_distance, dispersion_data$prop_IHC_tissue_corrected, method = "kendall")
cor.test(dispersion_data$StdDev_mean_NN_distance, dispersion_data$log10_mean_quantity_bacteria, method = "kendall")

### !
# no logarithmic model between StdDev_mean_NN_distance & log10_mean_quantity_bacteria
### ! 

### measured_mean_NN_distance -----------------------------------------------
(g1 <- ggplot(dispersion_data, aes(x=measured_mean_NN_distance, fill=group)) +
   geom_histogram(binwidth = 500, color='black', position="dodge") +
   labs(x = "Mean nearest neighbor distance (\u00b5m)", y ="Count") +
   theme_minimal())
ggsave("distribution_measured_mean_NN_distance_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

(g2 <- ggplot(dispersion_data , aes(x=reorder(code, measured_mean_NN_distance), y=measured_mean_NN_distance, fill=group)) + 
  geom_errorbar(aes(ymin=measured_mean_NN_distance - StdDev_mean_NN_distance, ymax=measured_mean_NN_distance + StdDev_mean_NN_distance), width=.3) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Images", y ="Mean nearest neighbor distance (\u00b5m)") +
  scale_y_continuous(limits=c(0, 10000)) +
  theme_minimal() +
  theme(legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          axis.title.x = element_text( size=15),
          axis.title.y = element_text(size=15)))
g2 + scale_fill_discrete(name = "Infection group", labels = c("Advanced infection","Control","Early infection"))
ggsave("measured_mean_NN_distance_group_reordered.png", path = output_path, width = 29.7, height = 21, units = "cm")

# Merge g1 & g2
ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("distribution_measured_mean_NN_distance_group_reordered.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(dispersion_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                                   "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                                   "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                                   "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                            y=measured_mean_NN_distance, colour = group)) +
  geom_point(size=3) +
  labs(x = "Images", y ="Mean nearest neighbor distance (\u00b5m)")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

summary(dispersion_data$measured_mean_NN_distance)
summary(dispersion_data$measured_mean_NN_distance[dispersion_data$group == "control"])
summary(dispersion_data$measured_mean_NN_distance[dispersion_data$group == "early_infection"])
summary(dispersion_data$measured_mean_NN_distance[dispersion_data$group == "advanced_infection"])

shapiro.test(dispersion_data$measured_mean_NN_distance) # p-value < 0.05 <- not normally distributed
hist(dispersion_data$measured_mean_NN_distance)
qqnorm(dispersion_data$measured_mean_NN_distance)
qqline(dispersion_data$measured_mean_NN_distance)

wilcox.test(x = dispersion_data$measured_mean_NN_distance[dispersion_data$group== "advanced_infection"], mu = 0, alternative = "greater")
wilcox.test(x = dispersion_data$measured_mean_NN_distance[dispersion_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = dispersion_data$measured_mean_NN_distance[dispersion_data$group== "control"], mu = 0, alternative = "greater")

kruskal.test(measured_mean_NN_distance ~ group, data = dispersion_data)
dunnTest(measured_mean_NN_distance ~ group, data = dispersion_data, method="bonferroni")

# Correlation test
cor.test(dispersion_data$measured_mean_NN_distance, dispersion_data$log10_mean_quantity_bacteria, method = "kendall")

### !
# no logarithmic model between measured_mean_NN_distance & log10_mean_quantity_bacteria
### ! 

### Dispersion index (t-test) -----------------------------------------------
(p5 <- ggplot(dispersion_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                                   "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                                   "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                                   "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                            y=t, colour = group)) +
  geom_point(size=3) +
  labs(x = "Images", y ="t-test")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)))
ggsave("particle_dispersion_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

ggplot(dispersion_data, aes(x=factor(code, level=c("C_J1_1","C_J1_2","C_J1_3",
                                                   "C_J4_1","C_J4_2","C_J4_3","C_J4_4",
                                                   "I_J1_1","I_J1_2","I_J1_3","I_J1_4","I_J1_5","I_J1_6","I_J1_7","I_J1_8","I_J1_9","I_J1_10","I_J1_11","I_J1_12","I_J1_13","I_J1_14","I_J1_15","I_J1_16",
                                                   "I_J4_1","I_J4_2","I_J4_3","I_J4_4","I_J4_5","I_J4_6","I_J4_7","I_J4_8","I_J4_9","I_J4_10","I_J4_11","I_J4_12","I_J4_13","I_J4_14","I_J4_15","I_J4_16","I_J4_17")),
                            y=t, colour = treatment)) +
  geom_point(size=3) +
  labs(x = "Images", y ="t-test")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.text.y = element_text(angle=45, vjust=1, hjust=1) +
  theme_minimal())
ggsave("particle_dispersion_treatment.png", path = output_path, width = 29.7, height = 21, units = "cm")

summary(dispersion_data$t)
summary(dispersion_data$t[dispersion_data$group == "control"])
summary(dispersion_data$t[dispersion_data$group == "early_infection"])
summary(dispersion_data$t[dispersion_data$group == "advanced_infection"])

shapiro.test(dispersion_data$t)
hist(dispersion_data$t)
qqnorm(dispersion_data$t)
qqline(dispersion_data$t)

wilcox.test(x = dispersion_data$t[dispersion_data$group== "advanced_infection"], mu = 0, alternative = "greater")
wilcox.test(x = dispersion_data$t[dispersion_data$group== "early_infection"], mu = 0, alternative = "greater")
wilcox.test(x = dispersion_data$t[dispersion_data$group== "control"], mu = 0, alternative = "greater")

kruskal.test(t ~ group, data = dispersion_data)
dunnTest(t ~ group, data = dispersion_data, method="bonferroni")

### Dispersion index (t-test) ~ prop_IHC_tissue_corrected -------------------------------------------
(p6 <- ggplot(dispersion_data, aes(x=t, y=prop_IHC_tissue_corrected))+
  geom_point(aes(colour = group)) +
  labs(x="t value", y="Proportion of IHC stained bacteria area (%)") +
  theme_minimal())
ggsave("t_prop_IHC_area_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

# Merge p5 & p6
ggarrange(p5, p6, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("t_particle_dispersion_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

cor.test(dispersion_data$t, dispersion_data$prop_IHC_tissue_corrected, method = "kendall")

# Logarithmic regression model (model <- lm(y ~ log(x)))
model3 <- lm(dispersion_data$prop_IHC_tissue_corrected ~ log(dispersion_data$t))
summary(model3)
coef(model3)

res3 <- resid(model3)
hist(res3) # check normality of the residuals 
plot(fitted(model3), res3) # check variance of the residuals
abline(0,0)

### Dispersion index (t-test) ~ log10_mean_quantity_bacteria -------------------------------------------
(p7 <- ggplot(dispersion_data, aes(x=t, y=log10_mean_quantity_bacteria))+
    geom_point(aes(colour = group)) +
    labs(x="t value", y="Mean quantity bacteria (log10)") +
    theme_minimal())
ggsave("t_mean_quantity_bacteria_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

# Merge p5 & p7
ggarrange(p5, p7, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("t_mean_quantity_bacteria_group.png", path = output_path, width = 29.7, height = 21, units = "cm")

cor.test(dispersion_data$t, dispersion_data$log10_mean_quantity_bacteria, method = "kendall")

# Logarithmic regression model (model <- lm(y ~ log(x)))
model4 <- lm(dispersion_data$log10_mean_quantity_bacteria ~ log(dispersion_data$t))
summary(model4)
coef(model4)

res4 <- resid(model4)
hist(res4) # check normality of the residuals 
plot(fitted(model4), res4) # check variance of the residuals
abline(0,0)
