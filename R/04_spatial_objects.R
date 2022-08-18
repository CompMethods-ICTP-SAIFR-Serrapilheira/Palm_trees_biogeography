####### Writing began in August 2022 by Thales
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script treats the shapefile and occurrence data, outing them in the proper
####### former to be used in the spp_biogeo_obtention

####### The shapefiles are not in the GitHub repository for the sake of memory saving
####### Nevertheless they are all just fragments of a larger shapefile presented by
####### Olson et al. 2001 Terrestrial Ecorregions of the World: A New Map of Life on Earth
####### and the instructions to access it are in their article

# Library
library(raster)
library(dplyr)

# Importing biogeographical regions shapefiles and species occurrence points
Neotropic <- shapefile("data/GIS/Neotropic/Neotropic.shp")
Indo_Malay <- shapefile("data/GIS/Indo-Malay/Indo-Malay.shp")
Australasia <- shapefile("data/GIS/Australasia/Australasia_corrected.shp")
Afrotropic <- shapefile("data/GIS/Afrotropic/Afrotropic.shp")

coor_spp_clean <- read.csv("./data/processed/coor_spp_clean.csv")

# Merging all 4 biogeographic regions polygons --------------------------------------------
t <- list(Neotropic@polygons[[1]], Afrotropic@polygons[[1]],
          Indo_Malay@polygons[[1]], Australasia@polygons[[1]])
ID <- c("A", "B", "C", "D")

for (i in 1:length(ID)) {
  t[[i]]@ID <- ID[i]
}

t_p  <- SpatialPolygons(list(t[[1]], t[[2]], t[[3]], t[[4]]),
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

df <- data.frame(bio_reg = c("Neotropic", "Afrotropic",
                             "Indo_Malay", "Australasia"), row.names = c("A", "B", "C", "D"))

bio_regs <- SpatialPolygonsDataFrame(t_p, data = as.data.frame(df))



# Creating a spatial object for the occurence points (including subtribe information) ------------------------------------------------
spatial_coor_spp <- coor_spp_clean %>%
  dplyr::select(sp_lon, sp_lat) %>%
  sp::SpatialPoints(, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

df1 <- data.frame(subtribes = coor_spp_clean$sp_subtribe,
                  row.names = c(1:length(coor_spp_clean$sp_lat)))

spatial_coor_spp <- SpatialPointsDataFrame(spatial_coor_spp, data = as.data.frame(df1))


# writing spatial objects of occurrences points and biogeographic regions polygons
if (!dir.exists("./output/GIS/")) {
  dir.create("./output/GIS/")
}

shapefile(spatial_coor_spp, "./output/GIS/spatial_coord_spp.shp")
shapefile(bio_regs, "./output/GIS/bio_regs.shp")
