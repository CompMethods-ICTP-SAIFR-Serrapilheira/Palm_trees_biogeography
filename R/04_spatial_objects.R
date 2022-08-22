####### Writing began in August 2022 by Thales de Lima
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script creates spatial objects for the occurrence records of Calamoideae
####### subfamily species and for the biogrographic regions inhabited by them.
####### These spatial objects are used in the 05_spp_biogeo_obtention.R script

####### The shapefile used here with the biogeographic regions of the world
####### is from the article Olson et al.2001 Terrestrial Ecorregions of the World:
####### A New Map of Life on Earth. It is not in the GitHub repository for the sake
####### of memory saving (it takes up 62.3MB of memory).
####### Nevertheless it can be downloaded from the following website:
####### https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
####### The downloaded file will be a zip file. You must unzip it and
####### copy the folder "official" to inside the folder data/GIS/WWF_terr_ecos/

####### The part that creates the biogeographic spatial object (bio_regs.shp) is
####### all commented, and I suggest skipping this part and just downloading the
####### processed shapefile from the GitHub repository. I say this because the
####### original shapefile from WWF is way to heavy and might take I large time to
####### import it to R environment

# Library
library(raster)
library(dplyr)
library(maptools)

### Creating shapefile for biogeographic regions inhabited by subfamily Calamoideae ----------------

# Reading the WWF shapefile
#WWF <- shapefile("data/GIS/WWF_terr_ecos/official/wwf_terr_ecos.shp")

# Extracting and merging the polygons of the Neotropic, Afrotropic, Indo-Malay and
# Australasia biogeographic regions
#NT_poly <- raster::subset(WWF, WWF$REALM == "NT")
#NT_poly_m <- maptools::unionSpatialPolygons(NT_poly, NT_poly$REALM == "NT")


#AT_poly <- raster::subset(WWF, WWF$REALM == "AT")
#AT_poly_m <- maptools::unionSpatialPolygons(AT_poly, AT_poly$REALM == "AT")


#IM_poly <- raster::subset(WWF, WWF$REALM == "IM")
#IM_poly_m <- maptools::unionSpatialPolygons(IM_poly, IM_poly$REALM == "IM")


#AA_poly <- raster::subset(WWF, WWF$REALM == "AA")
#AA_poly_m <- maptools::unionSpatialPolygons(AA_poly, AA_poly$REALM == "AA")



# Joining all 4 biogeographic regions polygons into one spatial object
#t <- list(NT_poly_m@polygons[[1]], AT_poly_m@polygons[[1]],
#          IM_poly_m@polygons[[1]], AA_poly_m@polygons[[1]])
#ID <- c("A", "B", "C", "D")

#for (i in 1:length(ID)) {
#  t[[i]]@ID <- ID[i]
#}

#t_p  <- SpatialPolygons(list(t[[1]], t[[2]], t[[3]], t[[4]]),
#                        proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

#df <- data.frame(bio_reg = c("Neotropic", "Afrotropic",
#                             "Indo_Malay", "Australasia"), row.names = c("A", "B", "C", "D"))

#bio_regs <- SpatialPolygonsDataFrame(t_p, data = as.data.frame(df))


### Creating a spatial object for the occurence points (including subtribe information) ------------------------------------------------
coor_spp_clean <- read.csv("./data/processed/coor_spp_clean.csv")

spatial_coor_spp <- coor_spp_clean %>%
  dplyr::select(sp_lon, sp_lat) %>%
  sp::SpatialPoints(, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

df1 <- data.frame(subtribes = coor_spp_clean$sp_subtribe,
                  row.names = c(1:length(coor_spp_clean$sp_lat)))

spatial_coor_spp <- SpatialPointsDataFrame(spatial_coor_spp, data = as.data.frame(df1))


### Writing spatial objects of occurrences points and biogeographic regions polygons -------------
if (!dir.exists("./output/GIS/")) {
  dir.create("./output/GIS/")
}

#shapefile(bio_regs, "./output/GIS/bio_regs.shp") # if you have runned the commented shapefile
# codelines in the first session, uncomment this line to save the generated shapefile
# As said, this shapefile is already in the GitHub repository.
shapefile(spatial_coor_spp, "./output/GIS/spatial_coord_spp.shp")
