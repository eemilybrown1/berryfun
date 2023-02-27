library("tidyverse")
library("sp")
library("dplyr")
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library("here")
library("sf")

setwd(here("data"))

##Loading and prepping data

#Get bioclimatic variable data from worldclim 
bioclim_data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5, #Using this coarse resolution for now because our range spans multiple SRTM tiles
                        path = here("data"))
# DON'T PUSH THIS DOWNLOAD, TOO BIG FOR GITHUB

#Loading GBIF csvs
gbif_aa <- read.csv(file = "amelanchieralnifolia.csv")
gbif_rl <- read.csv(file = "rubuslasiococcus.csv")
gbif_rn <- read.csv(file = "rubusnivalis.csv")
gbif_vp <- read.csv(file = "vacciniumparvifolium.csv")

#Dropping values with no Lat/Lon data and limiting to PNW
obs_aa <- gbif_aa[!is.na(gbif_aa$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))
obs_rl <- gbif_rl[!is.na(gbif_rl$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))
obs_rn <- gbif_rn[!is.na(gbif_rn$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))
obs_vp <- gbif_vp[!is.na(gbif_vp$decimalLatitude), ] %>%
  filter(stateProvince %in% c("British Columbia", "Washington", "Oregon"))

#Making lat/lon only datasets
latlon_aa <- obs_aa %>% dplyr::select(decimalLongitude, decimalLatitude)

latlon_rl <- obs_rl %>% dplyr::select(decimalLongitude, decimalLatitude)

latlon_rn <- obs_rn %>% dplyr::select(decimalLongitude, decimalLatitude)

latlon_vp <- obs_vp %>% dplyr::select(decimalLongitude, decimalLatitude)

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


#Clip climateNA data to geographic extent
bioclim_data_aa <- crop(x = bioclim_data, y = geographic_extent_aa)

bc_model_aa <- bioclim(x = bioclim_data_aa, p=latlon_aa)

#Predict presence from model
predict_presence_aa <- dismo::predict(object = bc_model_aa,
                                      x = bioclim_data_aa,
                                      ext = geographic_extent_aa)

#plotting predicted presence
plot(wrld_simpl,
     xlim = c(min_lon_aa, max_lon_aa),
     ylim = c(min_lat_aa, max_lat_aa),
     axes = TRUE,
     col = "grey95")

plot(predict_presence_aa, add = TRUE)

plot(wrld_simpl, add = TRUE, border = "grey5")

points(x = obs_aa$decimalLongitude,
       y = obs_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

box()

#Adding in pseudo-absence points to improve model 

#Use bioclim data files for sampling resolution
bil_files_aa <- list.files(path = "wc2-5",
                           pattern = "*.bil$",
                           full.names = TRUE)

#only need one file, so just choose the first one in the list
mask_aa <- raster(bil_files_aa[1])

#Setting the seed for random number generator
set.seed(20230226)

# Randomly sample points (same number as our observed points)
background_aa <- randomPoints(mask = mask_aa,     # Provides resolution of sampling points
                           n = nrow(latlon_aa),      # Number of random points
                           ext = geographic_extent_aa, # Spatially restricts sampling
                           extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_aa)

#Can visualize pseudo-absence points:
plot(wrld_simpl, 
     xlim = c(min_lon_aa, max_lon_aa),
     ylim = c(min_lat_aa, max_lat_aa),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

# Add the background points
points(background_aa, col = "grey30", pch = 1, cex = 0.75)

# Add the observations
points(x = obs_aa$decimalLongitude,
       y = obs_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

box()

#Need to build model using training data, and reserve som etesting data. 
#Reserve 20% of the data for testing. 
#Using the kfold function in the dismo package evenly assigns each observation to a random group.

testing_group_aa <- 1

group_presence_aa <- kfold(x=latlon_aa, k=5)

table(group_presence_aa)

# Separate observations into training and testing groups
presence_train_aa <- latlon_aa[group_presence_aa != testing_group_aa, ]
presence_test_aa <- latlon_aa[group_presence_aa == testing_group_aa, ]

# Repeat the process for pseudo-absence points
group_background_aa <- kfold(x = background_aa, k = 5)
background_train_aa <- background_aa[group_background_aa != testing_group_aa, ]
background_test_aa <- background_aa[group_background_aa == testing_group_aa, ]

##Training and testing the model

# Build a model using training data
bc_model_aa <- bioclim(x= bioclim_data_aa, p = presence_train_aa)

# Predict presence from model (same as previously, but with the update model)
predict_presence_aa <- dismo::predict(object = bc_model_aa, 
                                   x = bioclim_data_aa, 
                                   ext = geographic_extent_aa)

# Use testing data for model evaluation
bc_eval_aa <- evaluate(p = presence_test_aa,   # The presence testing data
                    a = background_test_aa, # The absence testing data
                    model = bc_model_aa,    # The model we are evaluating
                    x = bioclim_data_aa)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc_threshold_aa <- threshold(x = bc_eval_aa, stat = "spec_sens")

#We want to use the threshold to paint a map with the predicted range:
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_aa, max_lon_aa),
     ylim = c(min_lat_aa, max_lat_aa),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict_presence_aa > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# And add those observations
points(x = obs_aa$decimalLongitude, 
       y = obs_aa$decimalLatitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#The map drawn here presents a categorical classification of whether a particular point on the landscape will be suitable or not to the species of interest


#Next: forecasting distributions.








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


