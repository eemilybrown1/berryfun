library("tidyverse")
library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library("here")

setwd(here("data"))

##Loading and prepping data

#Get bioclimatic variable data from worldclim 
# bioclim_data <- getData(name = "worldclim",
#                        var = "bio",
#                        res = 2.5, #Using this coarse resolution for now because our range spans multiple SRTM tiles
#                        path = here("data"))
# **Edit: shouldn't need this anymore because Terrell has uploaded the climate data from ClimateNA.
# **Don't recommend running this code because it downloads a very large file. Just keeping it here in case we need it later.


#Loading GBIF csvs
obs_aa <- read.csv(file = "amelanchieralnifolia.csv")
obs_rl <- read.csv(file = "rubuslasiococcus.csv")
obs_rn <- read.csv(file = "rubusnivalis.csv")
obs_vp <- read.csv(file = "vacciniumparvifolium.csv")

#Dropping values with no Lat/Lon data and limiting to PNW
obs_aa <- obs_aa[!is.na(obs_aa$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))
obs_rl <- obs_rl[!is.na(obs_rl$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))
obs_rn <- obs_rn[!is.na(obs_rn$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))
obs_vp <- obs_vp[!is.na(obs_vp$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))



##Mapping A. alnifolia

#Determining geographic extent
max_lat_aa <- ceiling(max(obs_aa$decimalLatitude))
min_lat_aa <- floor(min(obs_aa$decimalLatitude))
max_lon_aa <- ceiling(max(obs_aa$decimalLongitude))
min_lon_aa <- floor(min(obs_aa$decimalLongitude))
geographic_extent_aa <- extent(x = c(min_lon_aa, max_lon_aa, min_lat_aa, max_lat_aa))

#Load data for base map
data(wrld_simpl)

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_aa, max_lon_aa),
     ylim = c(min_lat_aa, max_lat_aa),
     axes = TRUE, 
     col = "grey95")

#Add the points for individual observation
points(x = obs_aa$decimalLongitude, 
       y = obs_aa$decimalLatitude, 
       col = "#E7298A", 
       pch = 20, 
       cex = 0.75)
box()
title(main = "A. alnifolia observations")



##Mapping R. lasiococcus

#Determining geographic extent
max_lat_rl <- ceiling(max(obs_rl$decimalLatitude))
min_lat_rl <- floor(min(obs_rl$decimalLatitude))
max_lon_rl <- ceiling(max(obs_rl$decimalLongitude))
min_lon_rl <- floor(min(obs_rl$decimalLongitude))
geographic_extent_rl <- extent(x = c(min_lon_rl, max_lon_rl, min_lat_rl, max_lat_rl))

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_rl, max_lon_rl),
     ylim = c(min_lat_rl, max_lat_rl),
     axes = TRUE, 
     col = "grey95")

#Add the points for individual observation
points(x = obs_rl$decimalLongitude, 
       y = obs_rl$decimalLatitude, 
       col = "#1B9E77", 
       pch = 20, 
       cex = 0.75)
box()
title(main = "R. lasiococcus observations")



##Mapping R. nivalis

#Determining geographic extent
max_lat_rn <- ceiling(max(obs_rn$decimalLatitude))
min_lat_rn <- floor(min(obs_rn$decimalLatitude))
max_lon_rn <- ceiling(max(obs_rn$decimalLongitude))
min_lon_rn <- floor(min(obs_rn$decimalLongitude))
geographic_extent_rn <- extent(x = c(min_lon_rn, max_lon_rn, min_lat_rn, max_lat_rn))

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_rn, max_lon_rn),
     ylim = c(min_lat_rn, max_lat_rn),
     axes = TRUE, 
     col = "grey95")

#Add the points for individual observation
points(x = obs_rn$decimalLongitude, 
       y = obs_rn$decimalLatitude, 
       col = "#D95F02", 
       pch = 20, 
       cex = 0.75)
box()
title(main = "R. nivalis observations")



##Mapping V. parvifolium

#Determining geographic extent
max_lat_vp <- ceiling(max(obs_vp$decimalLatitude))
min_lat_vp <- floor(min(obs_vp$decimalLatitude))
max_lon_vp <- ceiling(max(obs_vp$decimalLongitude))
min_lon_vp <- floor(min(obs_vp$decimalLongitude))
geographic_extent_vp <- extent(x = c(min_lon_vp, max_lon_vp, min_lat_vp, max_lat_vp))

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_vp, max_lon_vp),
     ylim = c(min_lat_vp, max_lat_vp),
     axes = TRUE, 
     col = "grey95")

#Add the points for individual observation
points(x = obs_vp$decimalLongitude, 
       y = obs_vp$decimalLatitude, 
       col = "#7570B3", 
       pch = 20, 
       cex = 0.75)
box()
title(main = "V. parvifolium observations")



