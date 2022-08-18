####### Writing began in August 2022 by Thales
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script imports species occurrence data from GBIF and cleans it

# Library
library(rgbif)
library(phytools)
library(ape)
library(dplyr)
library(raster)
library(CoordinateCleaner)

# Importing tree and getting the terminal names (i.e., extant species)
tree <- read.nexus("./output/ultra_rerooted_tree.nex")

# Dropping the outgroup species
# The ancestral area reconstruction will be done only for the species of Calamoideae subfamily
tree_core <- ape::drop.tip(tree, tip = c("Nypa_fruticans",
                                         "Trachycarpus_martianus",
                                         "Pseudophoenix_vinifera"))
tree_core
plot(tree_core)

# Getting species names
spp_names <- tree_core$tip.label

# Download data from GBIF for all 36 species
occs_spp <- list()
for (i in 1:length(spp_names)) {
  species <- spp_names[i]
  occs_spp[[i]] <- occ_search(scientificName = species, # Specify that you are calling the species name
                     limit = 100000,           # Specify the limit of occurences
                     basisOfRecord = "PRESERVED_SPECIMEN")
  message(paste0("Importing GBIF date for species ", spp_names[i]))
}

# Importing species subtribe
subtribes <- read.csv("./data/raw/species_subtribe.csv", sep = ";")

for (i in 1:length(occs_spp)) {
  occs_spp[[i]]$data$scientificName_update <- spp_names[i]
  occs_spp[[i]]$data$subtribe <- subtribes$subtribo[i]
}

# Getting the coordinates (and removing the NAs) ------------------------------
coor_spp <- NULL

for (i in 1:length(occs_spp)) {
  species <- occs_spp[[i]]
  if ("decimalLatitude" %in% colnames(species$data)) {
    sp_subtribe <-  species$data$subtribe
    sp_name <- species$data$scientificName_update
    sp_lat <- species$data$decimalLatitude
    sp_lon <- species$data$decimalLongitude
    sp_country <- species$data$country
    sp_country[sp_country == "unknown or invalid"] <- NA
    sp_occ <- data.frame(sp_name, sp_subtribe, sp_lat, sp_lon)

    coor_spp <- rbind(coor_spp, sp_occ)
  }
}

coor_spp <- data.frame(na.omit(coor_spp))
unique(coor_spp$sp_name)
row.names(coor_spp) <- 1:length(coor_spp$sp_name)

write.csv(coor_spp, file = "./data/processed/coor_spp.csv")

# Saving the coordinates

## Getting the countries names (and removing the NAs) ------------------------------
spp_countries <- NULL

for (i in 1:length(occs_spp)) {
  species <- occs_spp[[i]]
  if ("country" %in% colnames(species$data)){
    sp_subtribe <- species$data$subtribe
    sp_name <- species$data$scientificName_update
    sp_country <- species$data$country
    sp_country[sp_country == "unknown or invalid"] <- NA
    sp_occ <- data.frame(sp_name, sp_subtribe, sp_country)
    spp_countries <- rbind(spp_countries, sp_occ)
  }
}


spp_countries <- na.omit(spp_countries)
unique(spp_countries$sp_name)
row.names(spp_countries) <- 1:length(spp_countries$sp_name)

write.csv(spp_countries, file = "./data/processed/countries_spp.csv")

## Cleaning data -----------------------------------------
coor_spp_clean <- clean_coordinates(coor_spp, lon = "sp_lon",
                                    lat = "sp_lat", species = "sp_name",
                                    value = "clean")

write.csv(coor_spp_clean, file = "./data/processed/coor_spp_clean.csv")
