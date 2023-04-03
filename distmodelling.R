#Top----
library("tidyverse")
library("sp")
library("dplyr")
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library("here")
library("sf")
library("geodata")
library(ggplot2)
library(gridExtra)
library(ggspatial)
library(terra)
library(tidyterra)
library(ggeasy)
library(ggtext)


setwd(here("data"))

#Loading and prepping data----
#Get bioclimatic variable data from worldclim 
bioclim_data <- worldclim_global(var = "bio",
                        res = 2.5, #Using this coarse resolution for now because our range spans multiple SRTM tiles
                        version = '2.1', 
                        path = here("data"))
bioclim_data <- brick(bioclim_data)

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

pnw <- raster(paste0("pnw_raster_10000.tif"))
pnw <- raster(paste0(here("data/pnw_raster_10000.tif")))

geographic_extent_pnw <- extent(pnw)

max_lat_pnw <- geographic_extent_pnw[4]
min_lat_pnw <- geographic_extent_pnw[3]
max_lon_pnw <- geographic_extent_pnw[2]
min_lon_pnw <- geographic_extent_pnw[1]

#also do this using the spat vector
bound <- c('Oregon', 'Washington', 'British Columbia')
pnw.bound <- gadm(country=c('USA', 'CAN'), level =1, path = here('data')) 
pnw.bound <- pnw.bound[pnw.bound$NAME_1 %in% bound,]
pnw.bound.ext <- ext(pnw.bound)
crs(pnw.bound)

#PNW Basemap----
basemap <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey', colour = 'black') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow() 

basemap

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

#CMIP6
#SSP 245 is intermediate emissions scenario, SSP 585 is highest emissions scenario
#using model "MRI-ESM2-0" based on evaluation from this paper: https://link.springer.com/article/10.1007/s00382-022-06410-1

#Forecasted climate data----
forecast_ssp245_year30 <- cmip6_world(model = "MRI-ESM2-0",
                                      ssp = "245",
                                      time = "2021-2040",
                                      var = "bioc",
                                      res = 2.5,
                                      path = here("data"))
forecast_ssp245_year50 <- cmip6_world(model = "MRI-ESM2-0",
                                      ssp = "245",
                                      time = "2041-2060",
                                      var = "bioc",
                                      res = 2.5,
                                      path = here("data"))

forecast_ssp245_year70 <- cmip6_world(model = "MRI-ESM2-0",
                                      ssp = "245",
                                      time = "2061-2080",
                                      var = "bioc",
                                      res = 2.5,
                                      path = here("data"))
forecast_ssp585_year30 <- cmip6_world(model = "MRI-ESM2-0",
                                      ssp = "585",
                                      time = "2021-2040",
                                      var = "bioc",
                                      res = 2.5,
                                      path = here("data"))

forecast_ssp585_year50 <- cmip6_world(model = "MRI-ESM2-0",
                                      ssp = "585",
                                      time = "2041-2060",
                                      var = "bioc",
                                      res = 2.5,
                                      path = here("data"))

forecast_ssp585_year70 <- cmip6_world(model = "MRI-ESM2-0",
                                      ssp = "585",
                                      time = "2061-2080",
                                      var = "bioc",
                                      res = 2.5,
                                      path = here("data"))

# DON'T PUSH THIS DOWNLOAD, TOO BIG FOR GITHUB

#Getting all the climate versions in the right format 
forecast_ssp245_year30 <- brick(forecast_ssp245_year30)
names(forecast_ssp245_year30) <- names(bioclim_data_pnw)

forecast_ssp245_year50 <- brick(forecast_ssp245_year50)
names(forecast_ssp245_year50) <- names(bioclim_data_pnw)

forecast_ssp245_year70 <- brick(forecast_ssp245_year70)
names(forecast_ssp245_year70) <- names(bioclim_data_pnw)

forecast_ssp585_year30 <- brick(forecast_ssp585_year30)
names(forecast_ssp585_year30) <- names(bioclim_data_pnw)

forecast_ssp585_year50 <- brick(forecast_ssp585_year50)
names(forecast_ssp585_year50) <- names(bioclim_data_pnw)

forecast_ssp585_year70 <- brick(forecast_ssp585_year70)
names(forecast_ssp585_year70) <- names(bioclim_data_pnw)

#A alnioflia
#A alnioflia modeling----
#Modelling aa using presence points in pnw extent

model_aa <- bioclim(x = bioclim_data_pnw, p=latlon_aa)

#Predict presence from model
predict_presence_aa <- dismo::predict(object = model_aa,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw) 

#Adding in pseudo-absence points to improve model 

# Randomly sample points (same number as our observed points)
background_aa <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_aa),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_aa)



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


# Predict presence from model (same as previously, but with the update model)

model_aa <- bioclim(x= bioclim_data_pnw, p = presence_train_aa)

# Predict presence from model (same as previously, but with the update model)
predict_presence_aa <- dismo::predict(object = model_aa, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

# Use testing data for model evaluation
eval_aa <- evaluate(p = presence_test_aa,   # The presence testing data
                    a = background_test_aa, # The absence testing data
                    model = model_aa,    # The model we are evaluating
                    x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
threshold_aa <- threshold(x = eval_aa, stat = "spec_sens")

#A alnifolia forecasting suitability----
forecast_presence_aa_245_30 <- dismo::predict(object = model_aa,
                                              x = forecast_ssp245_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_aa_245_50 <- dismo::predict(object = model_aa,
                                              x = forecast_ssp245_year50,
                                              ext = geographic_extent_pnw)
forecast_presence_aa_245_70 <- dismo::predict(object = model_aa,
                                              x = forecast_ssp245_year70,
                                              ext = geographic_extent_pnw)
forecast_presence_aa_585_30 <- dismo::predict(object = model_aa,
                                              x = forecast_ssp585_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_aa_585_50 <- dismo::predict(object = model_aa,
                                              x = forecast_ssp585_year50,
                                              ext = geographic_extent_pnw)
forecast_presence_aa_585_70 <- dismo::predict(object = model_aa,
                                              x = forecast_ssp585_year70,
                                              ext = geographic_extent_pnw)
#A. alnifolia ggplot----
#data prep for bellow
#historical
aa_currentpa <- predict_presence_aa>threshold_aa
aa_currentpa_rast <- terra::rast(aa_currentpa)
aa_currentpa_vect <- terra::as.polygons(aa_currentpa_rast)
aa_currentpa_crop<- crop(x=aa_currentpa_vect, y=pnw.bound)
aa_currentpa_mask <- mask(x=aa_currentpa_crop, mask=pnw.bound)
aa_currentpa_sub <- terra::subset(aa_currentpa_mask, aa_currentpa_mask$layer==1)
#ssp245 2030
aa_pa_24530 <- forecast_presence_aa_245_30>threshold_aa
aa_pa_24530_rast <- terra::rast(aa_pa_24530)
aa_pa_24530_vect <- terra::as.polygons(aa_pa_24530_rast)
aa_pa_24530_crop<- crop(x=aa_pa_24530_vect, y=pnw.bound)
aa_pa_24530_mask <- mask(x=aa_pa_24530_crop, mask=pnw.bound)
aa_pa_24530_sub <- terra::subset(aa_pa_24530_mask, aa_pa_24530_mask$layer==1)
#ssp245 2050
aa_pa_24550 <- forecast_presence_aa_245_50>threshold_aa
aa_pa_24550_rast <- terra::rast(aa_pa_24550)
aa_pa_24550_vect <- terra::as.polygons(aa_pa_24550_rast)
aa_pa_24550_crop<- crop(x=aa_pa_24550_vect, y=pnw.bound)
aa_pa_24550_mask <- mask(x=aa_pa_24550_crop, mask=pnw.bound)
aa_pa_24550_sub <- terra::subset(aa_pa_24550_mask, aa_pa_24550_mask$layer==1)
#ssp245 2070
aa_pa_24570 <- forecast_presence_aa_245_70>threshold_aa
aa_pa_24570_rast <- terra::rast(aa_pa_24570)
aa_pa_24570_vect <- terra::as.polygons(aa_pa_24570_rast)
aa_pa_24570_crop<- crop(x=aa_pa_24570_vect, y=pnw.bound)
aa_pa_24570_mask <- mask(x=aa_pa_24570_crop, mask=pnw.bound)
aa_pa_24570_sub <- terra::subset(aa_pa_24570_mask, aa_pa_24570_mask$layer==1)
#ssp585 2030
aa_pa_58530 <- forecast_presence_aa_585_30>threshold_aa
aa_pa_58530_rast <- terra::rast(aa_pa_58530)
aa_pa_58530_vect <- terra::as.polygons(aa_pa_58530_rast)
aa_pa_58530_crop<- crop(x=aa_pa_58530_vect, y=pnw.bound)
aa_pa_58530_mask <- mask(x=aa_pa_58530_crop, mask=pnw.bound)
aa_pa_58530_sub <- terra::subset(aa_pa_58530_mask, aa_pa_58530_mask$layer==1)
#ssp585 2050
aa_pa_58550 <- forecast_presence_aa_585_50>threshold_aa
aa_pa_58550_rast <- terra::rast(aa_pa_58550)
aa_pa_58550_vect <- terra::as.polygons(aa_pa_58550_rast)
aa_pa_58550_crop<- crop(x=aa_pa_58550_vect, y=pnw.bound)
aa_pa_58550_mask <- mask(x=aa_pa_58550_crop, mask=pnw.bound)
aa_pa_58550_sub <- terra::subset(aa_pa_58550_mask, aa_pa_58550_mask$layer==1)
#ssp585 2070
aa_pa_58570 <- forecast_presence_aa_585_70>threshold_aa
aa_pa_58570_rast <- terra::rast(aa_pa_58570)
aa_pa_58570_vect <- terra::as.polygons(aa_pa_58570_rast)
aa_pa_58570_crop<- crop(x=aa_pa_58570_vect, y=pnw.bound)
aa_pa_58570_mask <- mask(x=aa_pa_58570_crop, mask=pnw.bound)
aa_pa_58570_sub <- terra::subset(aa_pa_58570_mask, aa_pa_58570_mask$layer==1)

aa_observations <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey', colour = 'black') +
  geom_point(data=latlon_aa, aes(x=decimalLongitude, y=decimalLatitude ), colour = '#3333FF' ) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia')), paste('GBIF Research Grade Observations'))))+
  labs(x=NULL, y=NULL) +
  easy_center_title()

aa_observations


#Historical habitat suitability 
aa_historical <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = "Historical", values = c('Historical' = '#3333FF'), labels = 'Historical', name = 'Habitat Suitability') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia')), paste('Historical Habitat Suitability'))))+
  easy_center_title()

aa_historical

#ssp245 2030 suitability 
aa_ssp245_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_24530_sub, aes(fill= '2030'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030'), values = c('2030'= '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030')) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('2030 SSP245 Habitat Suitability'))))+
  easy_center_title()

aa_ssp245_2030

#ssp245 2050 suitability 
aa_ssp245_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_24550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2050'), values = c('2050'= '#ED83D1', 'Historical' = '#3333FF'), labels = c('Historical', '2050')) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('2050 SSP245 Habitat Suitability'))))+
  easy_center_title()

aa_ssp245_2050

#ssp245 2070 suitability 
aa_ssp245_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_24570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2070'), values = c('2070'= '#E7298A', 'Historical' = '#3333FF'), labels = c('Historical', '2070')) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('2070 SSP245 Habitat Suitability'))))+
  easy_center_title()

aa_ssp245_2070

#ssp245 arrangment
aa_ssp245_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_24530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=aa_pa_24550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=aa_pa_24570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('SSP245 Habitat Suitability'))))+
  easy_center_title()

aa_ssp245_legend

#ssp585 2030 suitability 
aa_ssp585_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_58530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030'), values = c('2030'= '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030')) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('2030 SSP585 Habitat Suitability'))))+
  easy_center_title()

aa_ssp585_2030

#ssp245 2050 suitability 
aa_ssp585_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_58550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2050'), values = c('2050'= '#ED83D1', 'Historical' = '#3333FF'), labels = c('Historical', '2050')) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('2050 SSP585 Habitat Suitability'))))+
  easy_center_title()

aa_ssp585_2050

#ssp585 2070 suitability 
aa_ssp585_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_58570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2070'), values = c('2070'= '#E7298A', 'Historical' = '#3333FF'), labels = c('Historical', '2070')) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('2070 SSP585 Habitat Suitability'))))+
  easy_center_title()

aa_ssp585_2070

#ssp585 arrangment
aa_ssp585_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=aa_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=aa_pa_58530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=aa_pa_58550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=aa_pa_58570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Amelanchier alnifolia'),' ','Forecasted'), paste('SSP585 Habitat Suitability'))))+
  easy_center_title()

aa_ssp585_legend

#function to take grob of legend to plot with individial slice plots 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

aa_legend <- get_legend(aa_ssp245_legend)

#grid arrangement of habitat suitability slices 
#ssp245
grid.arrange(aa_ssp245_2030, aa_ssp245_2050, aa_ssp245_2070, aa_legend,
             nrow = 1, ncol = 4)


#A. alnifolia area----
#calculate area  (km2) of each raster layer
aa_area_historical <- terra::expanse(x=aa_currentpa_sub, unit = 'km') %>% 
                        sum() 
aa_area_24530 <- terra::expanse(x=aa_pa_24530_sub, unit = 'km') %>% 
                        sum() 
aa_area_24550 <- terra::expanse(x=aa_pa_24550_sub, unit = 'km') %>% 
                        sum() 
aa_area_24570 <- terra::expanse(x=aa_pa_24570_sub, unit = 'km') %>% 
                        sum()
aa_area_58530 <- terra::expanse(x=aa_pa_58530_sub, unit = 'km') %>% 
                        sum()
aa_area_58550 <- terra::expanse(x=aa_pa_58550_sub, unit = 'km') %>% 
                        sum()
aa_area_58570 <- terra::expanse(x=aa_pa_58570_sub, unit = 'km') %>% 
                        sum()


##density plots for lat/lon distribution

aa_currentpa <- predict_presence_aa>threshold_aa 
aa_currentpapoints <- rasterToPoints(aa_currentpa, function(x)x==1)

aa_pa24530 <- forecast_presence_aa_245_30 > threshold_aa 
aa_pa24530points <- rasterToPoints(aa_pa24530, function(x)x==1)

aa_pa24550 <- forecast_presence_aa_245_50 > threshold_aa 
aa_pa24550points <- rasterToPoints(aa_pa24550, function(x)x==1)

aa_pa24570 <- forecast_presence_aa_245_70 > threshold_aa 
aa_pa24570points <- rasterToPoints(aa_pa24570, function(x)x==1)

aa_pa58530 <- forecast_presence_aa_585_30 > threshold_aa 
aa_pa58530points <- rasterToPoints(aa_pa58530, function(x)x==1)

aa_pa58550 <- forecast_presence_aa_585_50 > threshold_aa 
aa_pa58550points <- rasterToPoints(aa_pa58550, function(x)x==1)

aa_pa58570 <- forecast_presence_aa_585_70 > threshold_aa 
aa_pa58570points <- rasterToPoints(aa_pa58570, function(x)x==1)

#A. alnifolia density plots----
#Lattitude density ssp245
ggplot() +
  geom_density(aes(y = aa_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = aa_pa24530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = aa_pa24550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = aa_pa24570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP245") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp245
  ggplot() +
    geom_density(aes(x = aa_currentpapoints[,1], fill = "Historical"), alpha = 1) +
    #geom_density(aes(x = aa_pa24530points[,1], fill = "2030"), alpha = 0.2) +
    #geom_density(aes(x = aa_pa24550points[,1], fill = "2050"), alpha = 0.2) +
    geom_density(aes(x = aa_pa24570points[,1], fill = "2070"), alpha = .8) +
    #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP245") +
    xlab("Longitude") +
    ylab('Density') +
    theme_bw() +
    scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
    scale_x_continuous(breaks = c(-135, -125, -115)) +
    theme(legend.background = element_rect(colour = 'black'),
          legend.position = 'top',
          legend.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 12),
          title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   


#Lattitude density ssp585
  ggplot() +
    geom_density(aes(y = aa_currentpapoints[,2], fill = "Historical"), alpha = 1) +
    #geom_density(aes(y = aa_pa58530points[,2], fill = "2030"), alpha = 0.2) +
    #geom_density(aes(y = aa_pa58550points[,2], fill = "2050"), alpha = 0.2) +
    geom_density(aes(y = aa_pa58570points[,2], fill = "2070"), alpha = .8) +
    #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP585") +
    ylab("Latitude") +
    xlab('Density') +
    theme_bw() +
    scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
    theme(legend.background = element_rect(colour = 'black'),
          legend.position = 'top',
          legend.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 12),
          title = element_text(colour = 'black', size = 16)) +
    scale_y_continuous(breaks = c(45, 50, 55, 60)) +
    ylim(c(41, 60)) +
  easy_center_title()
  
#Longitude density ssp585
ggplot() +
    geom_density(aes(x = aa_currentpapoints[,1], fill = "Historical"), alpha = 1) +
    #geom_density(aes(x = aa_pa58530points[,1], fill = "2030"), alpha = 0.2) +
    #geom_density(aes(x = aa_pa58550points[,1], fill = "2050"), alpha = 0.2) +
    geom_density(aes(x = aa_pa58570points[,1], fill = "2070"), alpha = .8) +
    #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP585") +
    xlab("Longitude") +
    ylab('Density') +
    theme_bw() +
    scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
    theme(legend.background = element_rect(colour = 'black'),
          legend.position = 'top',
          legend.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 12),
          title = element_text(colour = 'black', size = 16)) +
    easy_center_title()   
  


#R. lasiococcus modeling----

model_rl <- bioclim(x = bioclim_data_pnw, p=latlon_rl)

#Predict presence from model
predict_presence_rl <- dismo::predict(object = model_rl,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw)



#Adding in pseudo-absence points to improve model 
# Randomly sample points (same number as our observed points)
background_rl <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_rl),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_rl)

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
model_rl <- bioclim(x= bioclim_data_pnw, p = presence_train_rl)

# Predict presence from model (same as previously, but with the update model)
predict_presence_rl <- dismo::predict(object = model_rl, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

eval_rl <- evaluate(p = presence_test_rl,   # The presence testing data
                       a = background_test_rl, # The absence testing data
                       model = model_rl,    # The model we are evaluating
                       x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
threshold_rl <- threshold(x = eval_rl, stat = "spec_sens")

#R. lasiococcus forecasting suitability----
forecast_presence_rl_245_30 <- dismo::predict(object = model_rl,
                                              x = forecast_ssp245_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_rl_245_50 <- dismo::predict(object = model_rl,
                                             x = forecast_ssp245_year50,
                                             ext = geographic_extent_pnw)
forecast_presence_rl_245_70 <- dismo::predict(object = model_rl,
                                             x = forecast_ssp245_year70,
                                             ext = geographic_extent_pnw)
forecast_presence_rl_585_30 <- dismo::predict(object = model_rl,
                                              x = forecast_ssp585_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_rl_585_50 <- dismo::predict(object = model_rl,
                                             x = forecast_ssp585_year50,
                                             ext = geographic_extent_pnw)
forecast_presence_rl_585_70 <- dismo::predict(object = model_rl,
                                             x = forecast_ssp585_year70,
                                             ext = geographic_extent_pnw)

#R. lasiococcus ggplot----
#data prep for bellow
#historical
rl_currentpa <- predict_presence_rl>threshold_rl
rl_currentpa_rast <- terra::rast(rl_currentpa)
rl_currentpa_vect <- terra::as.polygons(rl_currentpa_rast)
rl_currentpa_crop<- crop(x=rl_currentpa_vect, y=pnw.bound)
rl_currentpa_mask <- mask(x=rl_currentpa_crop, mask=pnw.bound)
rl_currentpa_sub <- terra::subset(rl_currentpa_mask, rl_currentpa_mask$layer==1)
#ssp245 2030
rl_pa_24530 <- forecast_presence_rl_245_30>threshold_rl
rl_pa_24530_rast <- terra::rast(rl_pa_24530)
rl_pa_24530_vect <- terra::as.polygons(rl_pa_24530_rast)
rl_pa_24530_crop<- crop(x=rl_pa_24530_vect, y=pnw.bound)
rl_pa_24530_mask <- mask(x=rl_pa_24530_crop, mask=pnw.bound)
rl_pa_24530_sub <- terra::subset(rl_pa_24530_mask, rl_pa_24530_mask$layer==1)
#ssp245 2050
rl_pa_24550 <- forecast_presence_rl_245_50>threshold_rl
rl_pa_24550_rast <- terra::rast(rl_pa_24550)
rl_pa_24550_vect <- terra::as.polygons(rl_pa_24550_rast)
rl_pa_24550_crop<- crop(x=rl_pa_24550_vect, y=pnw.bound)
rl_pa_24550_mask <- mask(x=rl_pa_24550_crop, mask=pnw.bound)
rl_pa_24550_sub <- terra::subset(rl_pa_24550_mask, rl_pa_24550_mask$layer==1)
#ssp245 2070
rl_pa_24570 <- forecast_presence_rl_245_70>threshold_rl
rl_pa_24570_rast <- terra::rast(rl_pa_24570)
rl_pa_24570_vect <- terra::as.polygons(rl_pa_24570_rast)
rl_pa_24570_crop<- crop(x=rl_pa_24570_vect, y=pnw.bound)
rl_pa_24570_mask <- mask(x=rl_pa_24570_crop, mask=pnw.bound)
rl_pa_24570_sub <- terra::subset(rl_pa_24570_mask, rl_pa_24570_mask$layer==1)
#ssp585 2030
rl_pa_58530 <- forecast_presence_rl_585_30>threshold_rl
rl_pa_58530_rast <- terra::rast(rl_pa_58530)
rl_pa_58530_vect <- terra::as.polygons(rl_pa_58530_rast)
rl_pa_58530_crop<- crop(x=rl_pa_58530_vect, y=pnw.bound)
rl_pa_58530_mask <- mask(x=rl_pa_58530_crop, mask=pnw.bound)
rl_pa_58530_sub <- terra::subset(rl_pa_58530_mask, rl_pa_58530_mask$layer==1)
#ssp585 2050
rl_pa_58550 <- forecast_presence_rl_585_50>threshold_rl
rl_pa_58550_rast <- terra::rast(rl_pa_58550)
rl_pa_58550_vect <- terra::as.polygons(rl_pa_58550_rast)
rl_pa_58550_crop<- crop(x=rl_pa_58550_vect, y=pnw.bound)
rl_pa_58550_mask <- mask(x=rl_pa_58550_crop, mask=pnw.bound)
rl_pa_58550_sub <- terra::subset(rl_pa_58550_mask, rl_pa_58550_mask$layer==1)
#ssp585 2070
rl_pa_58570 <- forecast_presence_rl_585_70>threshold_rl
rl_pa_58570_rast <- terra::rast(rl_pa_58570)
rl_pa_58570_vect <- terra::as.polygons(rl_pa_58570_rast)
rl_pa_58570_crop<- crop(x=rl_pa_58570_vect, y=pnw.bound)
rl_pa_58570_mask <- mask(x=rl_pa_58570_crop, mask=pnw.bound)
rl_pa_58570_sub <- terra::subset(rl_pa_58570_mask, rl_pa_58570_mask$layer==1)

#Historical habitat suitability 
rl_historical <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = "Historical", values = c('Historical' = '#FF0000'), labels = 'Historical', name = 'Habitat Suitability') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle('Historical Habitat\nSuitability (1970-2000)') +
  easy_center_title()

rl_historical

#ssp245 2030 suitability 
rl_ssp245_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rl_pa_24530_sub, fill = '#FFD5FF', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030'), values = c('2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('2030 SSP245 Habitat Suitability'))))+
  easy_center_title()

rl_ssp245_2030

#ssp245 2050 suitability 
rl_ssp245_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rl_pa_24550_sub, fill = '#ED83D1', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical",'2050'), values = c( '2050' = '#189E77', 'Historical' = '#FF0000'), labels = c('Historical', '2050'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('2050 SSP245 Habitat Suitability'))))+
  easy_center_title()

rl_ssp245_2050

#ssp245 2070 suitability 
rl_ssp245_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rl_pa_24570_sub, fill = '#E7298A', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2070'), values = c('2070' = '#0E543F','Historical' = '#FF0000'), labels = c('Historical',  '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('2070 SSP245 Habitat Suitability'))))+
  easy_center_title()

rl_ssp245_2070

#ssp245 arrangment
rl_ssp245_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=rl_pa_24530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=rl_pa_24550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=rl_pa_24570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#0E543F',  '2050' = '#189E77','2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('SSP245 Habitat Suitability'))))+
  easy_center_title()

rl_ssp245_legend

#ssp585 2030 suitability 
rl_ssp585_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rl_pa_58530_sub, fill = '#FFD5FF', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030'), values = c('2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('2030 SSP585 Habitat Suitability'))))+
  easy_center_title()

rl_ssp585_2030

#ssp245 2050 suitability 
rl_ssp585_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rl_pa_58550_sub, fill = '#ED83D1', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical",  '2050'), values = c('2050' = '#189E77', 'Historical' = '#FF0000'), labels = c('Historical', '2050' ), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('2050 SSP585 Habitat Suitability'))))+
  easy_center_title()

rl_ssp585_2050

#ssp585 2070 suitability 
rl_ssp585_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rl_pa_58570_sub, fill = '#E7298A', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  #scale_fill_manual(breaks = c("Historical", 'SSP585'), values = c('SSP585'= '#E7298A', 'Historical' = '#3333FF'), labels = c('Historical', 'SSP245'), name = 'Habitat Suitability') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2070 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

rl_ssp585_2070

#ssp585 arrangment
rl_ssp585_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rl_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=rl_pa_58530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=rl_pa_58550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=rl_pa_58570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#0E543F',  '2050' = '#189E77','2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus lasiococcus'),' ','Forecasted'), paste('SSP585 Habitat Suitability'))))+
  easy_center_title()

rl_ssp585_legend

#function to take grob of legend to plot with individial slice plots 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

rl_legend <- get_legend(rl_ssp245_legend)

#grid arrangement of habitat suitability slices 
#ssp245
grid.arrange(rl_ssp245_2030, rl_ssp245_2050, rl_ssp245_2070, rl_legend,
             nrow = 1, ncol = 4)


#R. lasiococcus area----
#calculate area  (km2) of each raster layer
rl_area_historical <- terra::expanse(x=rl_currentpa_sub, unit = 'km') %>% 
  sum() 
rl_area_24530 <- terra::expanse(x=rl_pa_24530_sub, unit = 'km') %>% 
  sum() 
rl_area_24550 <- terra::expanse(x=rl_pa_24550_sub, unit = 'km') %>% 
  sum() 
rl_area_24570 <- terra::expanse(x=rl_pa_24570_sub, unit = 'km') %>% 
  sum()
rl_area_58530 <- terra::expanse(x=rl_pa_58530_sub, unit = 'km') %>% 
  sum()
rl_area_58550 <- terra::expanse(x=rl_pa_58550_sub, unit = 'km') %>% 
  sum()
rl_area_58570 <- terra::expanse(x=rl_pa_58570_sub, unit = 'km') %>% 
  sum()


##density plots for lat/lon distribution

rl_currentpa <- predict_presence_rl>threshold_rl 
rl_currentpapoints <- rasterToPoints(rl_currentpa, function(x)x==1)

rl_pa24530 <- forecast_presence_rl_245_30 > threshold_rl 
rl_pa24530points <- rasterToPoints(rl_pa24530, function(x)x==1)

rl_pa24550 <- forecast_presence_rl_245_50 > threshold_rl 
rl_pa24550points <- rasterToPoints(rl_pa24550, function(x)x==1)

rl_pa24570 <- forecast_presence_rl_245_70 > threshold_rl 
rl_pa24570points <- rasterToPoints(rl_pa24570, function(x)x==1)

rl_pa58530 <- forecast_presence_rl_585_30 > threshold_rl 
rl_pa58530points <- rasterToPoints(rl_pa58530, function(x)x==1)

rl_pa58550 <- forecast_presence_rl_585_50 > threshold_rl 
rl_pa58550points <- rasterToPoints(rl_pa58550, function(x)x==1)

rl_pa58570 <- forecast_presence_rl_585_70 > threshold_rl 
rl_pa58570points <- rasterToPoints(rl_pa58570, function(x)x==1)

#R. lasiococcus density plots----
#Lattitude density ssp245
ggplot() +
  geom_density(aes(y = rl_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = rl_pa24530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = rl_pa24550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = rl_pa24570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP245") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#0E543F',  '2050' = '#189E77','2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp245
ggplot() +
  geom_density(aes(x = rl_currentpapoints[,1], fill = "Historical"), alpha = 1) +
  #geom_density(aes(x = rl_pa24530points[,1], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(x = rl_pa24550points[,1], fill = "2050"), alpha = 0.2) +
  geom_density(aes(x = rl_pa24570points[,1], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP245") +
  xlab("Longitude") +
  ylab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#0E543F',  '2050' = '#189E77','2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   


#Lattitude density ssp585
ggplot() +
  geom_density(aes(y = rl_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = rl_pa58530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = rl_pa58550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = rl_pa58570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP585") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#0E543F',  '2050' = '#189E77','2030' = '#98E4B8', 'Historical' = '#FF0000'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp585
ggplot() +
  geom_density(aes(x = rl_currentpapoints[,1], fill = "Historical"), alpha = 1) +
  #geom_density(aes(x = rl_pa58530points[,1], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(x = rl_pa58550points[,1], fill = "2050"), alpha = 0.2) +
  geom_density(aes(x = rl_pa58570points[,1], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP585") +
  xlab("Longitude") +
  ylab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#E7298A',  '2050' = '#ED83D1','2030' = '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   



#R. nivalis modeling----

model_rn <- bioclim(x = bioclim_data_pnw, p=latlon_rn)

#Predict presence from model
predict_presence_rn <- dismo::predict(object = model_rn,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw)

#Adding in pseudo-absence points to improve model 

# Randomly sample points (same number as our observed points)
background_rn <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_rn),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_rn)

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

model_rn <- bioclim(x= bioclim_data_pnw, p = presence_train_rn)

# Predict presence from model (same as previously, but with the update model)
predict_presence_rn <- dismo::predict(object = model_rn, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

# Use testing data for model evaluation

eval_rn <- evaluate(p = presence_test_rn,   # The presence testing data
                       a = background_test_rn, # The absence testing data
                       model = model_rn,    # The model we are evaluating
                       x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
threshold_rn <- threshold(x = eval_rn, stat = "spec_sens")

#R. nivalis forecasting suitability----
forecast_presence_rn_245_30 <- dismo::predict(object = model_rn,
                                              x = forecast_ssp245_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_rn_245_50 <- dismo::predict(object = model_rn,
                                             x = forecast_ssp245_year50,
                                             ext = geographic_extent_pnw)
forecast_presence_rn_245_70 <- dismo::predict(object = model_rn,
                                             x = forecast_ssp245_year70,
                                             ext = geographic_extent_pnw)
forecast_presence_rn_585_30 <- dismo::predict(object = model_rn,
                                              x = forecast_ssp585_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_rn_585_50 <- dismo::predict(object = model_rn,
                                             x = forecast_ssp585_year50,
                                             ext = geographic_extent_pnw)
forecast_presence_rn_585_70 <- dismo::predict(object = model_rn,
                                             x = forecast_ssp585_year70,
                                             ext = geographic_extent_pnw)

#R. nivalis ggplot----
#data prep for bellow
#historical
rn_currentpa <- predict_presence_rn>threshold_rn
rn_currentpa_rast <- terra::rast(rn_currentpa)
rn_currentpa_vect <- terra::as.polygons(rn_currentpa_rast)
rn_currentpa_crop<- crop(x=rn_currentpa_vect, y=pnw.bound)
rn_currentpa_mask <- mask(x=rn_currentpa_crop, mask=pnw.bound)
rn_currentpa_sub <- terra::subset(rn_currentpa_mask, rn_currentpa_mask$layer==1)
#ssp245 2030
rn_pa_24530 <- forecast_presence_rn_245_30>threshold_rn
rn_pa_24530_rast <- terra::rast(rn_pa_24530)
rn_pa_24530_vect <- terra::as.polygons(rn_pa_24530_rast)
rn_pa_24530_crop<- crop(x=rn_pa_24530_vect, y=pnw.bound)
rn_pa_24530_mask <- mask(x=rn_pa_24530_crop, mask=pnw.bound)
rn_pa_24530_sub <- terra::subset(rn_pa_24530_mask, rn_pa_24530_mask$layer==1)
#ssp245 2050
rn_pa_24550 <- forecast_presence_rn_245_50>threshold_rn
rn_pa_24550_rast <- terra::rast(rn_pa_24550)
rn_pa_24550_vect <- terra::as.polygons(rn_pa_24550_rast)
rn_pa_24550_crop<- crop(x=rn_pa_24550_vect, y=pnw.bound)
rn_pa_24550_mask <- mask(x=rn_pa_24550_crop, mask=pnw.bound)
rn_pa_24550_sub <- terra::subset(rn_pa_24550_mask, rn_pa_24550_mask$layer==1)
#ssp245 2070
rn_pa_24570 <- forecast_presence_rn_245_70>threshold_rn
rn_pa_24570_rast <- terra::rast(rn_pa_24570)
rn_pa_24570_vect <- terra::as.polygons(rn_pa_24570_rast)
rn_pa_24570_crop<- crop(x=rn_pa_24570_vect, y=pnw.bound)
rn_pa_24570_mask <- mask(x=rn_pa_24570_crop, mask=pnw.bound)
rn_pa_24570_sub <- terra::subset(rn_pa_24570_mask, rn_pa_24570_mask$layer==1)
#ssp585 2030
rn_pa_58530 <- forecast_presence_rn_585_30>threshold_rn
rn_pa_58530_rast <- terra::rast(rn_pa_58530)
rn_pa_58530_vect <- terra::as.polygons(rn_pa_58530_rast)
rn_pa_58530_crop<- crop(x=rn_pa_58530_vect, y=pnw.bound)
rn_pa_58530_mask <- mask(x=rn_pa_58530_crop, mask=pnw.bound)
rn_pa_58530_sub <- terra::subset(rn_pa_58530_mask, rn_pa_58530_mask$layer==1)
#ssp585 2050
rn_pa_58550 <- forecast_presence_rn_585_50>threshold_rn
rn_pa_58550_rast <- terra::rast(rn_pa_58550)
rn_pa_58550_vect <- terra::as.polygons(rn_pa_58550_rast)
rn_pa_58550_crop<- crop(x=rn_pa_58550_vect, y=pnw.bound)
rn_pa_58550_mask <- mask(x=rn_pa_58550_crop, mask=pnw.bound)
rn_pa_58550_sub <- terra::subset(rn_pa_58550_mask, rn_pa_58550_mask$layer==1)
#ssp585 2070
rn_pa_58570 <- forecast_presence_rn_585_70>threshold_rn
rn_pa_58570_rast <- terra::rast(rn_pa_58570)
rn_pa_58570_vect <- terra::as.polygons(rn_pa_58570_rast)
rn_pa_58570_crop<- crop(x=rn_pa_58570_vect, y=pnw.bound)
rn_pa_58570_mask <- mask(x=rn_pa_58570_crop, mask=pnw.bound)
rn_pa_58570_sub <- terra::subset(rn_pa_58570_mask, rn_pa_58570_mask$layer==1)

#NOTE NEED TO CHANGE FILL COLOURS OF MOST PLOTS BELLOW 
#Historical habitat suitability 
rn_historical <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = "Historical", values = c('Historical' = '#3333FF'), labels = 'Historical', name = 'Habitat Suitability') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle('Historical Habitat\nSuitability (1970-2000)') +
  easy_center_title()

rn_historical

#ssp245 2030 suitability 
rn_ssp245_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rn_pa_24530_sub, fill = '#FFD5FF', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2030 Forecasted Habitat\nSuitability - SSP245') +
#easy_center_title()

rn_ssp245_2030

#ssp245 2050 suitability 
rn_ssp245_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rn_pa_24550_sub, fill = '#ED83D1', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2050 Forecasted Habitat\nSuitability - SSP245') +
#easy_center_title()

rn_ssp245_2050

#ssp245 2070 suitability 
rn_ssp245_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rn_pa_24570_sub, fill = '#E7298A', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2070 Forecasted Habitat\nSuitability - SSP245') +
#easy_center_title()

rn_ssp245_2070

#ssp245 arrangment
rn_ssp245_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=rn_pa_24530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=rn_pa_24550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=rn_pa_24570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus nivalis'),' ','Forecasted'), paste('SSP245 Habitat Suitability'))))+
  easy_center_title()

rn_ssp245_legend

#ssp585 2030 suitability 
rn_ssp585_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rn_pa_58530_sub, fill = '#FFD5FF', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2030 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

rn_ssp585_2030

#ssp245 2050 suitability 
rn_ssp585_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rn_pa_58550_sub, fill = '#ED83D1', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2050 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

rn_ssp585_2050

#ssp585 2070 suitability 
rn_ssp585_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=rn_pa_58570_sub, fill = '#E7298A', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2070 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

rn_ssp585_2070

#ssp585 arrangment
rn_ssp585_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=rn_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=rn_pa_58530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=rn_pa_58550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=rn_pa_58570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Rubus nivalis'),' ','Forecasted'), paste('SSP585 Habitat Suitability'))))+
  easy_center_title()

rn_ssp585_legend

#function to take grob of legend to plot with individial slice plots 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

rn_legend <- get_legend(rn_ssp245_legend)

#grid arrangement of habitat suitability slices 
#ssp245
grid.arrange(rn_ssp245_2030, rn_ssp245_2050, rn_ssp245_2070, rn_legend,
             nrow = 1, ncol = 4)


#R. nivalis area----
#calculate area  (km2) of each raster layer
rn_area_historical <- terra::expanse(x=rn_currentpa_sub, unit = 'km') %>% 
  sum() 
rn_area_24530 <- terra::expanse(x=rn_pa_24530_sub, unit = 'km') %>% 
  sum() 
rn_area_24550 <- terra::expanse(x=rn_pa_24550_sub, unit = 'km') %>% 
  sum() 
rn_area_24570 <- terra::expanse(x=rn_pa_24570_sub, unit = 'km') %>% 
  sum()
rn_area_58530 <- terra::expanse(x=rn_pa_58530_sub, unit = 'km') %>% 
  sum()
rn_area_58550 <- terra::expanse(x=rn_pa_58550_sub, unit = 'km') %>% 
  sum()
rn_area_58570 <- terra::expanse(x=rn_pa_58570_sub, unit = 'km') %>% 
  sum()


##density plots for lat/lon distribution

rn_currentpa <- predict_presence_rn>threshold_rn 
rn_currentpapoints <- rasterToPoints(rn_currentpa, function(x)x==1)

rn_pa24530 <- forecast_presence_rn_245_30 > threshold_rn 
rn_pa24530points <- rasterToPoints(rn_pa24530, function(x)x==1)

rn_pa24550 <- forecast_presence_rn_245_50 > threshold_rn 
rn_pa24550points <- rasterToPoints(rn_pa24550, function(x)x==1)

rn_pa24570 <- forecast_presence_rn_245_70 > threshold_rn 
rn_pa24570points <- rasterToPoints(rn_pa24570, function(x)x==1)

rn_pa58530 <- forecast_presence_rn_585_30 > threshold_rn 
rn_pa58530points <- rasterToPoints(rn_pa58530, function(x)x==1)

rn_pa58550 <- forecast_presence_rn_585_50 > threshold_rn 
rn_pa58550points <- rasterToPoints(rn_pa58550, function(x)x==1)

rn_pa58570 <- forecast_presence_rn_585_70 > threshold_rn 
rn_pa58570points <- rasterToPoints(rn_pa58570, function(x)x==1)

#R. nivalis density plots----
#Lattitude density ssp245
ggplot() +
  geom_density(aes(y = rn_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = rn_pa24530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = rn_pa24550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = rn_pa24570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP245") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp245
ggplot() +
  geom_density(aes(x = rn_currentpapoints[,1], fill = "Historical"), alpha = 1) +
  #geom_density(aes(x = rn_pa24530points[,1], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(x = rn_pa24550points[,1], fill = "2050"), alpha = 0.2) +
  geom_density(aes(x = rn_pa24570points[,1], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP245") +
  xlab("Longitude") +
  ylab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   


#Lattitude density ssp585
ggplot() +
  geom_density(aes(y = rn_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = rn_pa58530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = rn_pa58550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = rn_pa58570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP585") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp585
ggplot() +
  geom_density(aes(x = rn_currentpapoints[,1], fill = "Historical"), alpha = 1) +
  #geom_density(aes(x = rn_pa58530points[,1], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(x = rn_pa58550points[,1], fill = "2050"), alpha = 0.2) +
  geom_density(aes(x = rn_pa58570points[,1], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP585") +
  xlab("Longitude") +
  ylab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#B34D00',  '2050' = '#DB6F1D','2030' = '#F0A461', 'Historical' = '#FFFF00'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   




#V. parvifolium mapping----
model_vp <- bioclim(x = bioclim_data_pnw, p=latlon_vp)

#Predict presence from model
predict_presence_vp <- dismo::predict(object = model_vp,
                                      x = bioclim_data_pnw,
                                      ext = geographic_extent_pnw)


#Adding in pseudo-absence points to improve model 
# Randomly sample points (same number as our observed points)
background_vp <- randomPoints(mask = mask,     # Provides resolution of sampling points
                              n = nrow(latlon_vp),      # Number of random points
                              ext = geographic_extent_pnw, # Spatially restricts sampling
                              extf = 1.25)  # Expands sampling a little bit

#take a look at background object
head(background_vp)

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
model_vp <- bioclim(x= bioclim_data_pnw, p = presence_train_vp)

# Predict presence from model (same as previously, but with the update model)
predict_presence_vp <- dismo::predict(object = model_vp, 
                                      x = bioclim_data_pnw, 
                                      ext = geographic_extent_pnw)

# Use testing data for model evaluation
eval_vp <- evaluate(p = presence_test_vp,   # The presence testing data
                       a = background_test_vp, # The absence testing data
                       model = model_vp,    # The model we are evaluating
                       x = bioclim_data_pnw)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
threshold_vp <- threshold(x = eval_vp, stat = "spec_sens")


#V parvifoliumn forecasting suitability----
forecast_presence_vp_245_30 <- dismo::predict(object = model_vp,
                                              x = forecast_ssp245_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_vp_245_50 <- dismo::predict(object = model_vp,
                                             x = forecast_ssp245_year50,
                                             ext = geographic_extent_pnw)
forecast_presence_vp_245_70 <- dismo::predict(object = model_vp,
                                             x = forecast_ssp245_year70,
                                             ext = geographic_extent_pnw)
forecast_presence_vp_585_30 <- dismo::predict(object = model_vp,
                                              x = forecast_ssp585_year30,
                                              ext = geographic_extent_pnw)
forecast_presence_vp_585_50 <- dismo::predict(object = model_vp,
                                             x = forecast_ssp585_year50,
                                             ext = geographic_extent_pnw)
forecast_presence_vp_585_70 <- dismo::predict(object = model_vp,
                                             x = forecast_ssp585_year70,
                                             ext = geographic_extent_pnw)

#V parvifolium ggplot----
#data prep for bellow
#historical
vp_currentpa <- predict_presence_vp>threshold_vp
vp_currentpa_rast <- terra::rast(vp_currentpa)
vp_currentpa_vect <- terra::as.polygons(vp_currentpa_rast)
vp_currentpa_crop<- crop(x=vp_currentpa_vect, y=pnw.bound)
vp_currentpa_mask <- mask(x=vp_currentpa_crop, mask=pnw.bound)
vp_currentpa_sub <- terra::subset(vp_currentpa_mask, vp_currentpa_mask$layer==1)
#ssp245 2030
vp_pa_24530 <- forecast_presence_vp_245_30>threshold_vp
vp_pa_24530_rast <- terra::rast(vp_pa_24530)
vp_pa_24530_vect <- terra::as.polygons(vp_pa_24530_rast)
vp_pa_24530_crop<- crop(x=vp_pa_24530_vect, y=pnw.bound)
vp_pa_24530_mask <- mask(x=vp_pa_24530_crop, mask=pnw.bound)
vp_pa_24530_sub <- terra::subset(vp_pa_24530_mask, vp_pa_24530_mask$layer==1)
#ssp245 2050
vp_pa_24550 <- forecast_presence_vp_245_50>threshold_vp
vp_pa_24550_rast <- terra::rast(vp_pa_24550)
vp_pa_24550_vect <- terra::as.polygons(vp_pa_24550_rast)
vp_pa_24550_crop<- crop(x=vp_pa_24550_vect, y=pnw.bound)
vp_pa_24550_mask <- mask(x=vp_pa_24550_crop, mask=pnw.bound)
vp_pa_24550_sub <- terra::subset(vp_pa_24550_mask, vp_pa_24550_mask$layer==1)
#ssp245 2070
vp_pa_24570 <- forecast_presence_vp_245_70>threshold_vp
vp_pa_24570_rast <- terra::rast(vp_pa_24570)
vp_pa_24570_vect <- terra::as.polygons(vp_pa_24570_rast)
vp_pa_24570_crop<- crop(x=vp_pa_24570_vect, y=pnw.bound)
vp_pa_24570_mask <- mask(x=vp_pa_24570_crop, mask=pnw.bound)
vp_pa_24570_sub <- terra::subset(vp_pa_24570_mask, vp_pa_24570_mask$layer==1)
#ssp585 2030
vp_pa_58530 <- forecast_presence_vp_585_30>threshold_vp
vp_pa_58530_rast <- terra::rast(vp_pa_58530)
vp_pa_58530_vect <- terra::as.polygons(vp_pa_58530_rast)
vp_pa_58530_crop<- crop(x=vp_pa_58530_vect, y=pnw.bound)
vp_pa_58530_mask <- mask(x=vp_pa_58530_crop, mask=pnw.bound)
vp_pa_58530_sub <- terra::subset(vp_pa_58530_mask, vp_pa_58530_mask$layer==1)
#ssp585 2050
vp_pa_58550 <- forecast_presence_vp_585_50>threshold_vp
vp_pa_58550_rast <- terra::rast(vp_pa_58550)
vp_pa_58550_vect <- terra::as.polygons(vp_pa_58550_rast)
vp_pa_58550_crop<- crop(x=vp_pa_58550_vect, y=pnw.bound)
vp_pa_58550_mask <- mask(x=vp_pa_58550_crop, mask=pnw.bound)
vp_pa_58550_sub <- terra::subset(vp_pa_58550_mask, vp_pa_58550_mask$layer==1)
#ssp585 2070
vp_pa_58570 <- forecast_presence_vp_585_70>threshold_vp
vp_pa_58570_rast <- terra::rast(vp_pa_58570)
vp_pa_58570_vect <- terra::as.polygons(vp_pa_58570_rast)
vp_pa_58570_crop<- crop(x=vp_pa_58570_vect, y=pnw.bound)
vp_pa_58570_mask <- mask(x=vp_pa_58570_crop, mask=pnw.bound)
vp_pa_58570_sub <- terra::subset(vp_pa_58570_mask, vp_pa_58570_mask$layer==1)

#NOTE NEED TO CHANGE FILL COLOURS OF MOST PLOTS BELLOW 
#AND OBJECT AND RASTER NAMES TO vp OBJECTS
#Historical habitat suitability 
vp_historical <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = "Historical", values = c('Historical' = '#3333FF'), labels = 'Historical', name = 'Habitat Suitability') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle('Historical Habitat\nSuitability (1970-2000)') +
  easy_center_title()

vp_historical

#ssp245 2030 suitability 
vp_ssp245_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=vp_pa_24530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030'), values = c('2030'= '#FFD5FF', 'Historical' = '#3333FF'), labels = c('Historical', '2030'), name = 'Habitat Suitability') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2030 Forecasted Habitat\nSuitability - SSP245') +
#easy_center_title()

vp_ssp245_2030

#ssp245 2050 suitability 
vp_ssp245_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=vp_pa_24550_sub, fill = '#ED83D1', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2050 Forecasted Habitat\nSuitability - SSP245') +
#easy_center_title()

vp_ssp245_2050

#ssp245 2070 suitability 
vp_ssp245_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=vp_pa_24570_sub, fill = '#E7298A', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2070 Forecasted Habitat\nSuitability - SSP245') +
#easy_center_title()

vp_ssp245_2070

#ssp245 arrangment
vp_ssp245_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=vp_pa_24530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=vp_pa_24550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=vp_pa_24570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Vaccinium parvifolium'),' ','Forecasted'), paste('SSP245 Habitat Suitability'))))+
  easy_center_title()

vp_ssp245_legend

#ssp585 2030 suitability 
vp_ssp585_2030 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=vp_pa_58530_sub, fill = '#FFD5FF', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2030 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

vp_ssp585_2030

#ssp245 2050 suitability 
vp_ssp585_2050 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=vp_pa_58550_sub, fill = '#ED83D1', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2050 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

vp_ssp585_2050

#ssp585 2070 suitability 
vp_ssp585_2070 <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, fill = '#3333FF', colour = 'transparent') +
  geom_sf(data=vp_pa_58570_sub, fill = '#E7298A', colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = c(.2, .3),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) #+
#ggtitle('2070 Forecasted Habitat\nSuitability - SSP585') +
#easy_center_title()

vp_ssp585_2070

#ssp585 arrangment
vp_ssp585_legend <- ggplot() + geom_sf(data=pnw.bound, fill = 'grey') +
  geom_sf(data=vp_currentpa_sub, aes(fill = 'Historical'), colour = 'transparent') +
  geom_sf(data=vp_pa_58530_sub, aes(fill = '2030'), colour = 'transparent') +
  geom_sf(data=vp_pa_58550_sub, aes(fill = '2050'), colour = 'transparent') +
  geom_sf(data=vp_pa_58570_sub, aes(fill = '2070'), colour = 'transparent') +
  geom_sf(data=pnw.bound, colour = 'black', fill = 'transparent') +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  annotation_north_arrow(style = north_arrow_orienteering, pad_y = unit(0.7, 'cm')) +
  annotation_scale(width_hint=.3) +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size =12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  ggtitle(expression(atop(paste(italic('Vaccinium parvifolium'),' ','Forecasted'), paste('SSP585 Habitat Suitability'))))+
  easy_center_title()

vp_ssp585_legend

#function to take grob of legend to plot with individial slice plots 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

vp_legend <- get_legend(vp_ssp245_legend)

#grid arrangement of habitat suitability slices 
#ssp245
grid.arrange(vp_ssp245_2030, vp_ssp245_2050, vp_ssp245_2070, vp_legend,
             nrow = 1, ncol = 4)

#V parvifolium  area----
#calculate area  (km2) of each raster layer
vp_area_historical <- terra::expanse(x=vp_currentpa_sub, unit = 'km') %>% 
  sum() 
vp_area_24530 <- terra::expanse(x=vp_pa_24530_sub, unit = 'km') %>% 
  sum() 
vp_area_24550 <- terra::expanse(x=vp_pa_24550_sub, unit = 'km') %>% 
  sum() 
vp_area_24570 <- terra::expanse(x=vp_pa_24570_sub, unit = 'km') %>% 
  sum()
vp_area_58530 <- terra::expanse(x=vp_pa_58530_sub, unit = 'km') %>% 
  sum()
vp_area_58550 <- terra::expanse(x=vp_pa_58550_sub, unit = 'km') %>% 
  sum()
vp_area_58570 <- terra::expanse(x=vp_pa_58570_sub, unit = 'km') %>% 
  sum()


##density plots for lat/lon distribution

vp_currentpa <- predict_presence_vp>threshold_vp 
vp_currentpapoints <- rasterToPoints(vp_currentpa, function(x)x==1)

vp_pa24530 <- forecast_presence_vp_245_30 > threshold_vp 
vp_pa24530points <- rasterToPoints(vp_pa24530, function(x)x==1)

vp_pa24550 <- forecast_presence_vp_245_50 > threshold_vp 
vp_pa24550points <- rasterToPoints(vp_pa24550, function(x)x==1)

vp_pa24570 <- forecast_presence_vp_245_70 > threshold_vp 
vp_pa24570points <- rasterToPoints(vp_pa24570, function(x)x==1)

vp_pa58530 <- forecast_presence_vp_585_30 > threshold_vp 
vp_pa58530points <- rasterToPoints(vp_pa58530, function(x)x==1)

vp_pa58550 <- forecast_presence_vp_585_50 > threshold_vp 
vp_pa58550points <- rasterToPoints(vp_pa58550, function(x)x==1)

vp_pa58570 <- forecast_presence_vp_585_70 > threshold_vp 
vp_pa58570points <- rasterToPoints(vp_pa58570, function(x)x==1)

#V parvifolium  density plots----
#Lattitude density ssp245
ggplot() +
  geom_density(aes(y = vp_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = vp_pa24530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = vp_pa24550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = vp_pa24570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP245") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp245
ggplot() +
  geom_density(aes(x = vp_currentpapoints[,1], fill = "Historical"), alpha = 1) +
  #geom_density(aes(x = vp_pa24530points[,1], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(x = vp_pa24550points[,1], fill = "2050"), alpha = 0.2) +
  geom_density(aes(x = vp_pa24570points[,1], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP245") +
  xlab("Longitude") +
  ylab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  scale_x_continuous(breaks = c(-135, -125, -115)) +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   


#Lattitude density ssp585
ggplot() +
  geom_density(aes(y = vp_currentpapoints[,2], fill = "Historical"), alpha = 1) +
  #geom_density(aes(y = vp_pa58530points[,2], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(y = vp_pa58550points[,2], fill = "2050"), alpha = 0.2) +
  geom_density(aes(y = vp_pa58570points[,2], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Latitudinal Distribution of A. alnifolia under SSP585") +
  ylab("Latitude") +
  xlab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_continuous(breaks = c(45, 50, 55, 60)) +
  ylim(c(41, 60)) +
  easy_center_title()

#Longitude density ssp585
ggplot() +
  geom_density(aes(x = vp_currentpapoints[,1], fill = "Historical"), alpha = 1) +
  #geom_density(aes(x = vp_pa58530points[,1], fill = "2030"), alpha = 0.2) +
  #geom_density(aes(x = vp_pa58550points[,1], fill = "2050"), alpha = 0.2) +
  geom_density(aes(x = vp_pa58570points[,1], fill = "2070"), alpha = .8) +
  #ggtitle("Change in Longitudinal Distribution of A. alnifolia under SSP585") +
  xlab("Longitude") +
  ylab('Density') +
  theme_bw() +
  scale_fill_manual(breaks = c("Historical", '2030', '2050', '2070'), values = c('2070' = '#310354',  '2050' = '#7570B3','2030' = '#A5A4ED', 'Historical' = '#FF6666'), labels = c('Historical', '2030', '2050', '2070'), name = '') +
  theme(legend.background = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  easy_center_title()   

#Bottom----