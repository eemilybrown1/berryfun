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

##PNW extent
pnw <- raster(paste0(here("data/pnw_raster_10000.tif")))
plot(pnw)

geographic_extent_pnw <- extent(pnw)

max_lat_pnw <- geographic_extent_pnw[4]
min_lat_pnw <- geographic_extent_pnw[3]
max_lon_pnw <- geographic_extent_pnw[2]
min_lon_pnw <- geographic_extent_pnw[1]

#Clip climateNA data to geographic extent of pnw
bioclim_data_pnw <- crop(x = bioclim_data, y = geographic_extent_pnw)

#Creating a mask from a wc file bil file
bil_files <- list.files(path = "wc2-5",
                        pattern = "*.bil$",
                        full.names = TRUE)

#only need one file, so just choose the first one in the list
mask <- raster(bil_files[1])

#Load data for base map
data(wrld_simpl)

#Setting the seed for random number generator
set.seed(20230227)


#####Mapping A. alnifolia

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
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


#Modelling aa using presence points in pnw extent
bc_model_aa <- bioclim(x = bioclim_data_pnw, p=latlon_aa)

#Predict presence from model
predict_presence_aa <- dismo::predict(object = bc_model_aa,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw) 

#plotting predicted presence with just presence points
plot(wrld_simpl,
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
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

# Randomly sample points (same number as our observed points)
background_aa <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(latlon_aa),      # Number of random points
                           ext = geographic_extent_pnw, # Spatially restricts sampling
                           extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_aa)

#Can visualize pseudo-absence points:
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
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
bc_model_aa <- bioclim(x= bioclim_data_pnw, p = presence_train_aa)

# Predict presence from model (same as previously, but with the update model)
predict_presence_aa <- dismo::predict(object = bc_model_aa, 
                                   x = bioclim_data_pnw, 
                                   ext = geographic_extent_pnw)

# Use testing data for model evaluation
bc_eval_aa <- evaluate(p = presence_test_aa,   # The presence testing data
                    a = background_test_aa, # The absence testing data
                    model = bc_model_aa,    # The model we are evaluating
                    x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc_threshold_aa <- threshold(x = bc_eval_aa, stat = "spec_sens")

#We want to use the threshold to paint a map with the predicted range:
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Predicted Landscape Suitability to A. alnifolia")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict_presence_aa > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

#And add those observations
points(x = obs_aa$decimalLongitude, 
       y = obs_aa$decimalLatitude, 
       col = "black",
       pch = "+", 
       cex = 0.4)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#The map drawn here presents a categorical classification of whether a particular point on the landscape will be suitable or not to the species of interest


#Next: forecasting distributions.
#Info on worldclim forecasting:
#CMIP5 uses 5th IPCC report- 6th not available in R yet
#RCP 26 is lowest emissions scenario, RCP 85 is highest emissions scenario)
#year should be 50 or 70
#'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO".



#Climate data 2021-2040
modelBC_rcp26_year50 <- getData(name = "CMIP5",
                                var = "bio",
                                res = 2.5,
                                model = "BC", #model BCC_CSM1.1
                                rcp = 26,
                                year = 50,
                                path = here("data"))

modelBC_rcp26_year70 <- getData(name = "CMIP5",
                                var = "bio",
                                res = 2.5,
                                model = "BC", #model BCC_CSM1.1
                                rcp = 26,
                                year = 70,
                                path = here("data"))

modelBC_rcp85_year50 <- getData(name = "CMIP5",
                                var = "bio",
                                res = 2.5,
                                model = "BC", #model BCC_CSM1.1
                                rcp = 85,
                                year = 50,
                                path = here("data"))

modelBC_rcp85_year70 <- getData(name = "CMIP5",
                                var = "bio",
                                res = 2.5,
                                model = "BC", #model BCC_CSM1.1
                                rcp = 85,
                                year = 70,
                                path = here("data"))

# DON'T PUSH THIS DOWNLOAD, TOO BIG FOR GITHUB

#RCP 26 year 50

modelBC_rcp26_year50 <- brick(modelBC_rcp26_year50)

names(modelBC_rcp26_year50) <- names(bioclim_data_pnw)

forecast_presence_aa_26_50 <- dismo::predict(object = bc_model_aa,
                                             x = modelBC_rcp26_year50,
                                             ext = geographic_extent_pnw)

plot(wrld_simpl,
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE,
     col = "grey95")

plot(forecast_presence_aa_26_50, add = TRUE)

plot(wrld_simpl, add=TRUE, border = "grey5")

points(x = latlon_aa$decimalLongitude,
       y = latlon_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

box()


#Plotting based on presence/absence
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_26_50 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()


#RCP 26 year 70

modelBC_rcp26_year70 <- brick(modelBC_rcp26_year70)

names(modelBC_rcp26_year70) <- names(bioclim_data_pnw)

forecast_presence_aa_26_70 <- dismo::predict(object = bc_model_aa,
                                             x = modelBC_rcp26_year70,
                                             ext = geographic_extent_pnw)

plot(wrld_simpl,
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE,
     col = "grey95")

plot(forecast_presence_aa_26_70, add = TRUE)

plot(wrld_simpl, add=TRUE, border = "grey5")

points(x = latlon_aa$decimalLongitude,
       y = latlon_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

box()


#Plotting based on presence/absence
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_26_70 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()


#RCP 85 year 50

modelBC_rcp85_year50 <- brick(modelBC_rcp85_year50)

names(modelBC_rcp85_year50) <- names(bioclim_data_pnw)

forecast_presence_aa_85_50 <- dismo::predict(object = bc_model_aa,
                                             x = modelBC_rcp85_year50,
                                             ext = geographic_extent_pnw)

plot(wrld_simpl,
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE,
     col = "grey95")

plot(forecast_presence_aa_85_50, add = TRUE)

plot(wrld_simpl, add=TRUE, border = "grey5")

points(x = latlon_aa$decimalLongitude,
       y = latlon_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

box()


#Plotting based on presence/absence
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_85_50 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()


#RCP 85 year 70

modelBC_rcp85_year70 <- brick(modelBC_rcp85_year70)

names(modelBC_rcp85_year70) <- names(bioclim_data_pnw)

forecast_presence_aa_85_70 <- dismo::predict(object = bc_model_aa,
                                             x = modelBC_rcp85_year70,
                                             ext = geographic_extent_pnw)

plot(wrld_simpl,
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE,
     col = "grey95")

plot(forecast_presence_aa_85_70, add = TRUE)

plot(wrld_simpl, add=TRUE, border = "grey5")

points(x = latlon_aa$decimalLongitude,
       y = latlon_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

box()


#Plotting based on presence/absence
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_85_70 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()



#Plotting current and future

par(mfrow=c(2,3))

#Panel 1
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Current Predicted Landscape Suitability")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict_presence_aa > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#Panel 2
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Forecast: RCP 26 Year 2050")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_26_50 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#Panel 3
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Forecast: RCP 26 Year 2070")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_26_70 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#Panel 4
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "GBIF species observations")

#Plotting point observations
points(x = latlon_aa$decimalLongitude,
       y = latlon_aa$decimalLatitude,
       col = "#E7298A",
       pch = 20,
       cex = 0.75)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()


#Panel 5
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Forecast: RCP 85 Year 2050")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_85_50 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#Panel 3
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Forecast: RCP 85 Year 2070")

# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence_aa_85_70 > bc_threshold_aa, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#E7298A"))

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()














##Mapping R. lasiococcus

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
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


#Clip climateNA data to geographic extent
bc_model_rl <- bioclim(x = bioclim_data_pnw, p=latlon_rl)

#Predict presence from model
predict_presence_rl <- dismo::predict(object = bc_model_rl,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw)

#plotting predicted presence with just presence points
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95")

plot(predict_presence_rl, add = TRUE)

plot(wrld_simpl, add = TRUE, border = "grey5")

points(x = obs_rl$decimalLongitude,
       y = obs_rl$decimalLatitude,
       col = "#1B9E77",
       pch = 20,
       cex = 0.75)

box()

#Adding in pseudo-absence points to improve model 

# Randomly sample points (same number as our observed points)
background_rl <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_rl),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_rl)

#Can visualize pseudo-absence points:
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

# Add the background points
points(background_rl, col = "grey30", pch = 1, cex = 0.75)

# Add the observations
points(x = obs_rl$decimalLongitude,
       y = obs_rl$decimalLatitude,
       col = "#1B9E77",
       pch = 20,
       cex = 0.75)

box()

#Need to build model using training data, and reserve som etesting data. 
#Reserve 20% of the data for testing. 
#Using the kfold function in the dismo package evenly assigns each observation to a random group.

testing_group_rl <- 1

group_presence_rl <- kfold(x=latlon_rl, k=5)

table(group_presence_rl)

# Separate observations into training and testing groups
presence_train_rl <- latlon_rl[group_presence_rl != testing_group_rl, ]
presence_test_rl <- latlon_rl[group_presence_rl == testing_group_rl, ]

# Repeat the process for pseudo-absence points
group_background_rl <- kfold(x = background_rl, k = 5)
background_train_rl <- background_rl[group_background_rl != testing_group_rl, ]
background_test_rl <- background_rl[group_background_rl == testing_group_rl, ]

##Training and testing the model

# Build a model using training data
bc_model_rl <- bioclim(x= bioclim_data_pnw, p = presence_train_rl)

# Predict presence from model (same as previously, but with the update model)
predict_presence_rl <- dismo::predict(object = bc_model_rl, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

# Use testing data for model evaluation
bc_eval_rl <- evaluate(p = presence_test_rl,   # The presence testing data
                       a = background_test_rl, # The absence testing data
                       model = bc_model_rl,    # The model we are evaluating
                       x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc_threshold_rl <- threshold(x = bc_eval_rl, stat = "spec_sens")

#We want to use the threshold to paint a map with the predicted range:
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Predicted Landscape Suitability to R. lasiococcus")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict_presence_rl > bc_threshold_rl, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#1B9E77"))

# And add those observations
points(x = obs_rl$decimalLongitude, 
       y = obs_rl$decimalLatitude, 
       col = "black",
       pch = "+", 
       cex = 0.3)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()



















##Mapping R. nivalis

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
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

#Model with just presence data
bc_model_rn <- bioclim(x = bioclim_data_pnw, p=latlon_rn)

#Predict presence from model
predict_presence_rn <- dismo::predict(object = bc_model_rn,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw)

#plotting predicted presence with just presence points
plot(wrld_simpl,
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE,
     col = "grey95")

plot(predict_presence_rn, add = TRUE)

plot(wrld_simpl, add = TRUE, border = "grey5")

points(x = obs_rn$decimalLongitude,
       y = obs_rn$decimalLatitude,
       col = "#D95F02",
       pch = 20,
       cex = 0.75)

box()

#Adding in pseudo-absence points to improve model 

# Randomly sample points (same number as our observed points)
background_rn <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_rn),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_rn)

#Can visualize pseudo-absence points:
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

# Add the background points
points(background_rn, col = "grey30", pch = 1, cex = 0.75)

# Add the observations
points(x = obs_rn$decimalLongitude,
       y = obs_rn$decimalLatitude,
       col = "#D95F02",
       pch = 20,
       cex = 0.75)

box()

#Need to build model using training data, and reserve som etesting data. 
#Reserve 20% of the data for testing. 
#Using the kfold function in the dismo package evenly assigns each observation to a random group.

testing_group_rn <- 1

group_presence_rn <- kfold(x=latlon_rn, k=5)

table(group_presence_rn)

# Separate observations into training and testing groups
presence_train_rn <- latlon_rn[group_presence_rn != testing_group_rn, ]
presence_test_rn <- latlon_rn[group_presence_rn == testing_group_rn, ]

# Repeat the process for pseudo-absence points
group_background_rn <- kfold(x = background_rn, k = 5)
background_train_rn <- background_rn[group_background_rn != testing_group_rn, ]
background_test_rn <- background_rn[group_background_rn == testing_group_rn, ]

##Training and testing the model

# Build a model using training data
bc_model_rn <- bioclim(x= bioclim_data_pnw, p = presence_train_rn)

# Predict presence from model (same as previously, but with the update model)
predict_presence_rn <- dismo::predict(object = bc_model_rn, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

# Use testing data for model evaluation
bc_eval_rn <- evaluate(p = presence_test_rn,   # The presence testing data
                       a = background_test_rn, # The absence testing data
                       model = bc_model_rn,    # The model we are evaluating
                       x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc_threshold_rn <- threshold(x = bc_eval_rn, stat = "spec_sens")

#We want to use the threshold to paint a map with the predicted range:
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Predicted Landscape Suitability to R. nivalis")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict_presence_rn > bc_threshold_rn, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#D95F02"))

# And add those observations
points(x = obs_rn$decimalLongitude, 
       y = obs_rn$decimalLatitude, 
       col = "black",
       pch = "+", 
       cex = 0.4)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

















##Mapping V. parvifolium

#Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
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

#Clip climateNA data to geographic extent
bc_model_vp <- bioclim(x = bioclim_data_pnw, p=latlon_vp)

#Predict presence from model
predict_presence_vp <- dismo::predict(object = bc_model_vp,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw)

#plotting predicted presence with just presence points
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95")

plot(predict_presence_rn, add = TRUE)

plot(wrld_simpl, add = TRUE, border = "grey5")

points(x = obs_vp$decimalLongitude,
       y = obs_vp$decimalLatitude,
       col = "#7570B3",
       pch = 20,
       cex = 0.75)

box()

#Adding in pseudo-absence points to improve model 

# Randomly sample points (same number as our observed points)
background_vp <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_vp),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_vp)

#Can visualize pseudo-absence points:
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

# Add the background points
points(background_vp, col = "grey30", pch = 1, cex = 0.75)

# Add the observations
points(x = obs_vp$decimalLongitude,
       y = obs_vp$decimalLatitude,
       col = "#7570B3",
       pch = 20,
       cex = 0.75)

box()

#Need to build model using training data, and reserve som etesting data. 
#Reserve 20% of the data for testing. 
#Using the kfold function in the dismo package evenly assigns each observation to a random group.

testing_group_vp <- 1

group_presence_vp <- kfold(x=latlon_vp, k=5)

table(group_presence_vp)

# Separate observations into training and testing groups
presence_train_vp <- latlon_vp[group_presence_vp != testing_group_vp, ]
presence_test_vp <- latlon_vp[group_presence_vp == testing_group_vp, ]

# Repeat the process for pseudo-absence points
group_background_vp <- kfold(x = background_vp, k = 5)
background_train_vp <- background_vp[group_background_vp != testing_group_vp, ]
background_test_vp <- background_vp[group_background_vp == testing_group_vp, ]

##Training and testing the model

# Build a model using training data
bc_model_vp <- bioclim(x= bioclim_data_pnw, p = presence_train_vp)

# Predict presence from model (same as previously, but with the update model)
predict_presence_vp <- dismo::predict(object = bc_model_vp, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

# Use testing data for model evaluation
bc_eval_vp <- evaluate(p = presence_test_vp,   # The presence testing data
                       a = background_test_vp, # The absence testing data
                       model = bc_model_vp,    # The model we are evaluating
                       x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc_threshold_vp <- threshold(x = bc_eval_vp, stat = "spec_sens")

#We want to use the threshold to paint a map with the predicted range:
# Plot base map
plot(wrld_simpl, 
     xlim = c(min_lon_pnw, max_lon_pnw),
     ylim = c(min_lat_pnw, max_lat_pnw),
     axes = TRUE, 
     col = "grey95",
     main = "Predicted Landscape Suitability to V. parvifolium")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict_presence_vp > bc_threshold_vp, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "#7570B3"))

# And add those observations
points(x = obs_vp$decimalLongitude, 
       y = obs_vp$decimalLatitude, 
       col = "black",
       pch = "+", 
       cex = 0.3)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()
