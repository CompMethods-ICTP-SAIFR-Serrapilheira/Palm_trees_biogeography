####### Writing began in August 2022 by Thales
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script uses the spatial elements of occurence and biogeographic regions created
####### by the scritp 04_spatial_objects and determines the biogeographical area occupied by
####### each species

# Library
library(raster)
library(dplyr)
library(ape)
library(tmap)


###### Ploting species occurrences and visualizing their distribution ----------------------

# Importing clean species data
coor_spp_clean <- read.csv("data/processed/coor_spp_clean.csv")

# Importing spatial objects
spatial_coord_spp <- raster::shapefile("output/GIS/spatial_coord_spp.shp")
spatial_bio_regs <- raster::shapefile("output/GIS/bio_regs.shp")

# Getting species subtribe information
n_tribes <- unique(coor_spp_clean$sp_subtribe)
tribes <- unique(coor_spp_clean$sp_subtribe)

# Importing the biogeographic regions shapefiles
MyPalette1 <- c("darkorange", "red", "yellow", "darkgreen")
MyPalette2 <- c("#00E5FF", "#FF8F00", "#FF3D00", "#000000",
                "#76FF03", "#0F00C4", "#6200EA", "#FFEE07")

# Writing the png map with the biogeographic regions and subtribes species occurences
pdf("./figs/t_subfamilies_maps.pdf", width = 7, height = 12)
p <- tm_shape(shapefile("data/GIS/World_Boundaries/World_Countries__Generalized_.shp")) +
  tm_borders(col = "grey") +
  tm_shape(spatial_bio_regs) +
  tm_fill('bio_reg',
          palette = MyPalette1, alpha = 0.2) +
  tm_shape(spatial_coord_spp) +
  tm_symbols(size = 0.2, shape = 21, col = "subtribes",
             palette = MyPalette2) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom",
            legend.text.size = 1.5, legend.title.color = "white",
            legend.stack = "horizontal")
p
dev.off()



###### Getting each species biogeographic regions ---------------------------------------

# Overlying occurrence points and biogeographic regions shapefile to extract each species
# biogeographic areas
bio_regs <- sp::over(spatial_coord_spp, spatial_bio_regs)

# Joining species biographic vector with the occurence vector
spp_coord_bio_reg <- cbind(coor_spp_clean, bio_regs)
spp_coord_bio_reg <- na.omit(spp_coord_bio_reg) # removing the NAs

# Getting species names
spp_names <- unique(spp_coord_bio_reg$sp_name)

# Creating a presence/absence matrix relating species and biogeographic regions
spp_bio_reg <- NULL
br <- NULL
for (i in 1:length(spp_names)) {
  br  <- spp_coord_bio_reg %>%
    filter(sp_name == spp_names[i]) %>%
    select(bio_reg) %>%
    unique() %>%
    as.vector() %>%
    unlist()
  if ("Afrotropic" %in% br) {
    br_vec <- 1
  } else {br_vec <- 0}

  if ("Australasia" %in% br) {
    br_vec[2] <- 1
  } else {br_vec[2] <- 0}

  if ("Indo_Malay" %in% br) {
    br_vec[3] <- 1
  } else {br_vec[3] <- 0}

  if ("Neotropic" %in% br) {
    br_vec[4] <- 1
  } else {br_vec[4] <- 0}
  spp_bio_reg <- rbind(spp_bio_reg, br_vec) %>%
                    as.data.frame()
}
spp_bio_reg <- cbind(spp_names, spp_bio_reg)

# Naming the column appropriately
bio_reg_names <- unique(bio_regs) %>%
                  na.omit() %>%
                  as.vector() %>%
                  unlist() %>%
                  sort

for (i in 1:length(bio_reg_names)) {
  names(spp_bio_reg)[i+1] <- bio_reg_names[i]
}

# Checking names
head(spp_bio_reg)


###### Assigning biogeographic areas to species that had no georeferenced occurrence records
# Some species had no occurrence records with latitude/longitude
# Let us check which of the initial species in the tree are not in the final biogeographic file

# Getting the name of all species in the Calamoideae subfamily tree used in this project
tree <- read.nexus("./output/ultra_rerooted_tree.nex")
tree_core <- ape::drop.tip(tree, tip = c("Nypa_fruticans",
                                         "Trachycarpus_martianus",
                                         "Pseudophoenix_vinifera"))
spp_names <- tree_core$tip.label

# Checking which species had no records and were excluded during the previous steps
coorless_spp <- spp_names[!spp_names %in% spp_bio_reg$spp_names]

# One of the options now is to check in which countries these species occur and than
# assign a biographic area with this information

# Importing country data
countries_spp <- read.csv("data/processed/countries_spp.csv")

countries_spp %>%
  filter(sp_name == coorless_spp[1])

countries_spp %>%
  filter(sp_name == coorless_spp[2])

countries_spp %>%
  filter(sp_name == coorless_spp[3])


# Checking in the map in which biogeographic regios the countries where the species occur are
World <- shapefile("data/GIS/World_Boundaries/World_Countries__Generalized_.shp")
World$isN <- ifelse(World$COUNTRY=="Malaysia", "darkorange", "white")
World$isN[World$COUNTRY == "Indonesia"] <- "purple"
World$isN[World$COUNTRY == "Philippines"] <- "blue"

tmap_mode("view")
tm_shape(spatial_bio_regs) +
  tm_fill('bio_reg',
          palette = MyPalette1) +
  tm_shape(World) +
  tm_fill("isN", alpha = 0.7)

# Using the tmap in the interactive mode we can see that:
# Malaysia and the Philippines are within the Indo-Malay region
# While Indonesia territory spans both Indo-Malay and Australasia regions

# Attributing the biogeographic regions to the three coordinate less species
# using the information about their countries

# Calamus calospathus and Salacca ramosiana
countries_spp %>%
  filter(sp_name == coorless_spp[1])
sp1 <- data.frame("Calamus_calospathus",0, 0, 1, 0)
for (i in 1:length(bio_reg_names)) {
  names(sp1)[i+1] <- bio_reg_names[i]
}
names(sp1)[1] <- "spp_names"

countries_spp %>%
  filter(sp_name == coorless_spp[3])
sp2 <- data.frame("Salacca_ramosiana",0, 0, 1, 0)
for (i in 1:length(bio_reg_names)) {
  names(sp2)[i+1] <- bio_reg_names[i]
}
names(sp2)[1] <- "spp_names"


# Pigafetta_elata
countries_spp %>%
  filter(sp_name == coorless_spp[2])
sp3 <- data.frame("Pigafetta_elata",0, 1, 1, 0)
for (i in 1:length(bio_reg_names)) {
  names(sp3)[i+1] <- bio_reg_names[i]
}
names(sp3)[1] <- "spp_names"


##### Now we can join all the biogeographic information for all species
spp_bio_reg <- rbind(spp_bio_reg, sp1, sp2, sp3)

# Saving the file with biogeography areas of each species
write.table(spp_bio_reg, "./output/spp_bio_reg.txt")
