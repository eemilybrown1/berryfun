########
#Terrell Roulston
#terrell.roulston@ubc.ca
#Started Mar  09, 2023
#ZOOL/RES509
#Berry Bunch Project
#MaxEnt SDM
########
#Referenced from 'ENMeval and Maxent for Species Distribution Modeling'
#https://rsh249.github.io/bioinformatics/ENMeval.html
library(raster)
library(ggplot2)
library(spocc)
library(dplyr)
library(sp)
library(tidyverse)
library(here)
library(geodata)
library(ENMeval)
library(spThin)
library(maptools)
library(viridis)
library(tidyterra)
library(dismo)
library(ggpubr)
library(qpdf)
#setwd to data to load occurrence data
setwd(here('data'))


#species codes
#'aa' = Amelanchier alnifolia
#'rl' = Rubus lasiococcus
#'rn' = Rubus nivalis
#'vp' = Vaccinium parvifolium

#call GBIF csvs
aa.occ <- read.csv(file = 'amelanchieralnifolia.csv') %>% filter(stateProvince %in% c('British Columbia', 'Washington', 'Oregon')) %>% #filter data found within PNW
  dplyr::select(decimalLongitude, decimalLatitude, species) %>% #select only lat and long and species (sp needed for thinning function)
  filter(!is.na(decimalLatitude)) #drop occurrences without lat/long

rl.occ <- read.csv(file = 'rubuslasiococcus.csv') %>% filter(stateProvince %in% c('British Columbia', 'Washington', 'Oregon')) %>% #filter data found within PNW
  dplyr::select(decimalLongitude, decimalLatitude, species) %>% #select only lat and long and species (sp needed for thinning function)
  filter(!is.na(decimalLatitude)) #drop occurrences without lat/long

rn.occ <- read.csv(file = 'rubusnivalis.csv') %>% filter(stateProvince %in% c('British Columbia', 'Washington', 'Oregon')) %>% #filter data found within PNW
  dplyr::select(decimalLongitude, decimalLatitude, species) %>% #select only lat and long and species (sp needed for thinning function)
  filter(!is.na(decimalLatitude)) #drop occurrences without lat/long

vp.occ <- read.csv(file = 'vacciniumparvifolium.csv') %>% filter(stateProvince %in% c('British Columbia', 'Washington', 'Oregon')) %>% #filter data found within PNW
  dplyr::select(decimalLongitude, decimalLatitude, species) %>% #select only lat and long and species (sp needed for thinning function)
  filter(!is.na(decimalLatitude)) #drop occurrences without lat/long

#extract world clim data
#WARNING: DO NOT COMMIT CLIMATE DATA - TOO LARGE FOR GITHUB
wc <- worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = here('data')) #download wc - data 1970-2000 'baseline'

#CMIP6 - SSP 245
cmip6.245.2040 <- cmip6_world(model='GISS-E2-1-H', ssp='245',time='2021-2040', var='bioc', res = 2.5, path = here('data')) #CMG is NASA's GISS-E2-1-H - cited as performing best
cmip6.245.2060 <- cmip6_world(model='GISS-E2-1-H', ssp='245',time='2041-2060', var='bioc', res = 2.5, path = here('data')) #CMG is NASA's GISS-E2-1-H - cited as performing best
cmip6.245.2080 <- cmip6_world(model='GISS-E2-1-H', ssp='245',time='2061-2080', var='bioc', res = 2.5, path = here('data')) #CMG is NASA's GISS-E2-1-H - cited as performing best

#CMIP6 - SSP 585
cmip6.585.2040 <- cmip6_world(model='GISS-E2-1-H', ssp='585',time='2021-2040', var='bioc', res = 2.5, path = here('data')) #CMG is NASA's GISS-E2-1-H - cited as performing best
cmip6.585.2060 <- cmip6_world(model='GISS-E2-1-H', ssp='585',time='2041-2060', var='bioc', res = 2.5, path = here('data')) #CMG is NASA's GISS-E2-1-H - cited as performing best
cmip6.585.2080 <- cmip6_world(model='GISS-E2-1-H', ssp='585',time='2061-2080', var='bioc', res = 2.5, path = here('data')) #CMG is NASA's GISS-E2-1-H - cited as performing best


pnw <- raster(paste0(here('data/pnw_raster_10000.tif'))) #import raster tiff for geo ref
geographic_extent_pnw <- extent(pnw) #determine extent for world clim cropping
max_lat_pnw <- geographic_extent_pnw[4]
min_lat_pnw <- geographic_extent_pnw[3]
max_lon_pnw <- geographic_extent_pnw[2]
min_lon_pnw <- geographic_extent_pnw[1]

wc.pnw <- crop(x = wc, y = geographic_extent_pnw) #crop wc data
cmip6.245.2040.pnw <- crop(x = cmip6.245.2040, y = geographic_extent_pnw)
cmip6.245.2060.pnw <- crop(x = cmip6.245.2060, y = geographic_extent_pnw)
cmip6.245.2080.pnw <- crop(x = cmip6.245.2080, y = geographic_extent_pnw)
cmip6.585.2040.pnw <- crop(x = cmip6.585.2040, y = geographic_extent_pnw)
cmip6.585.2060.pnw <- crop(x = cmip6.585.2060, y = geographic_extent_pnw)
cmip6.585.2080.pnw <- crop(x = cmip6.585.2080, y = geographic_extent_pnw)
#plot(wc.pnw) #check extent is correct

#thinning occurrence data to avoid sampling biases
#thinning to keep only one record in each cell of the climate data
data(wrld_simpl)

coordinates(aa.occ) <- ~decimalLongitude+decimalLatitude #defining coordinates
crs(aa.occ) <- crs(wrld_simpl) #align crs of world map with occurences

#plot(wc.pnw[[1]], legend=T) #plot wc variable 1 (mean annual temp)
#plot(aa.occ, pch = 20, add = T) #plot occurrences

set.seed(1337)
aa.occ.thin <- as.data.frame(gridSample(aa.occ, wc.pnw[[1]], n=1)) #keep one sample per cell

#visualize points retained points to ensure thinning worked correctly
wrld_crop <- crop(x=wrld_simpl, y=geographic_extent_pnw)

wc.pnw.crop <- crop(wc.pnw, extent(c(-123.6, -123.4), #you will need to zoom in to see points that are kept/dropped
                                   ylim = c(48.4, 48.6)))

#plot(wc.pnw.crop[[1]], xlim =c(-123.6, -123.4), ylim = c(48.4, 48.6))
#points(aa.occ, pch = 20, add = T)
#points(aa.occ.thin, col = 'blue', cex = 2) #blue is points that are kept
#thinning runs successfully

bound <- c('Oregon', 'Washington', 'British Columbia')
pnw.bound <- gadm(country=c('USA', 'CAN'), level =1, path = here('data')) 
pnw.bound <-pnw.bound[pnw.bound$NAME_1 %in% bound,]

#GGPLOT 
#first need to convert raster to df for ggplot
wc.pnw.df.all <- terra::as.data.frame(wc.pnw, xy = T)
cmip6.245.2040.pnw.df <- terra::as.data.frame(cmip6.245.2040.pnw, xy = T)
cmip6.245.2060.pnw.df <- terra::as.data.frame(cmip6.245.2060.pnw, xy = T)
cmip6.245.2080.pnw.df <- terra::as.data.frame(cmip6.245.2080.pnw, xy = T)
cmip6.585.2040.pnw.df <- terra::as.data.frame(cmip6.585.2040.pnw, xy = T)
cmip6.585.2060.pnw.df <- terra::as.data.frame(cmip6.585.2060.pnw, xy = T)
cmip6.585.2080.pnw.df <- terra::as.data.frame(cmip6.585.2080.pnw, xy = T)
#make theme for plots
theme_maps <- theme(
  panel.background = element_rect(fill='transparent', colour = 'black', size = 1.5),
  plot.background = element_rect(fill='transparent', color = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_rect(fill='transparent'),
  legend.box.background = element_rect(fill='transparent'),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = c(0.15, 0.3),
  axis.text = element_text(size = 4, colour = 'black'),
  axis.ticks = element_line(colour = 'black'),
  plot.margin = margin(t=10, r=-100000 ,b=10,l=-100000, unit = 'pt'),
  legend.title = element_text(size=4),
  legend.text = element_text(size=4),
  legend.key.size = unit(.2, 'cm'),
  plot.title = element_text(size = 4)
)
#PLOTTING
#plotting 1970-2000 climatic variables
#####
bio.1 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000') 

bio.2 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.3 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.4 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.5 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.6 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.7 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.8 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')
#can see there is a strong effect of mountain ranges - should this be dropped?

bio.9 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.10 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.11 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.12 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.13 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.14 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.15 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.16 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.17 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.18 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')

bio.19 <- ggplot() + geom_raster(data = wc.pnw.df.all, aes(x=x, y=y, fill= wc2.1_2.5m_bio_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('1970-2000')
#arrange plots together
temp_plots <- ggarrange(bio.1, bio.2, bio.3, bio.4, bio.5, bio.6, bio.7, bio.8, bio.9, bio.10, bio.11, 
                       ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots <- ggarrange(bio.12, bio.13, bio.14, bio.15, bio.16, bio.17, bio.18, bio.19, 
                        ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

#save plots
lapply(names(temp_plots), function(x)ggsave(filename=paste(x, 'temp_plots', '.pdf', sep =''), plot = temp_plots[[x]], path = here('figures/world clim/wc_1970_2000')))
lapply(names(perc_plots), function(x)ggsave(filename=paste(x, 'perc_plots', '.pdf', sep =''), plot = perc_plots[[x]], path = here('figures/world clim/wc_1970_2000')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/wc_1970_2000/1temp_plots.pdf'), here('figures/world clim/wc_1970_2000/2temp_plots.pdf'), here('figures/world clim/wc_1970_2000/3temp_plots.pdf'), 
                      here('figures/world clim/wc_1970_2000/1perc_plots.pdf'), here('figures/world clim/wc_1970_2000/2perc_plots.pdf')),
            output = here('figures/world clim/wc_1970_2000/1970-2000_bioc.pdf'))
#####
#CMIP6 SSP245 - 2021-2040
#SSP_245_2040
#####
ssp245.2040.bio.1 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.2 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.3 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.4 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.5 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.6 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.7 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.8 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')
#can see there is a strong effect of mountain ranges - should this be dropped?

ssp245.2040.bio.9 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.10 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.11 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.12 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.13 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.14 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.15 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.16 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.17 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.18 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040')

ssp245.2040.bio.19 <- ggplot() + geom_raster(data = cmip6.245.2040.pnw.df, aes(x=x, y=y, fill= wc2_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2040') 

#arrange plots together
temp_plots_245_2040 <- ggarrange(ssp245.2040.bio.1, ssp245.2040.bio.2, ssp245.2040.bio.3, ssp245.2040.bio.4, ssp245.2040.bio.5, ssp245.2040.bio.6, ssp245.2040.bio.7, ssp245.2040.bio.8, ssp245.2040.bio.9, ssp245.2040.bio.10, ssp245.2040.bio.11, 
                        ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots_245_2040 <- ggarrange(ssp245.2040.bio.12, ssp245.2040.bio.13, ssp245.2040.bio.14, ssp245.2040.bio.15, ssp245.2040.bio.16, ssp245.2040.bio.17, ssp245.2040.bio.18, ssp245.2040.bio.19, 
                        ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
#save plots
lapply(names(temp_plots_245_2040), function(x)ggsave(filename=paste(x, 'temp_plots_245_2040', '.pdf', sep =''), plot = temp_plots_245_2040[[x]], path = here('figures/world clim/ssp245_2040')))
lapply(names(perc_plots_245_2040), function(x)ggsave(filename=paste(x, 'perc_plots_245_2040', '.pdf', sep =''), plot = perc_plots_245_2040[[x]], path = here('figures/world clim/ssp245_2040')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/ssp245_2040/1temp_plots_245_2040.pdf'), here('figures/world clim/ssp245_2040/2temp_plots_245_2040.pdf'), here('figures/world clim/ssp245_2040/3temp_plots_245_2040.pdf'), 
                      here('figures/world clim/ssp245_2040/1perc_plots_245_2040.pdf'), here('figures/world clim/ssp245_2040/2perc_plots_245_2040.pdf')),
            output = here('figures/world clim/ssp245_2040/wc_245_2040.pdf'))
#####
#CMIP6 SSP245 - 2041-2060
#SSP_245_2060
#####
ssp245.2060.bio.1 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.2 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.3 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.4 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.5 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.6 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.7 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.8 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')
#can see there is a strong effect of mountain ranges - should this be dropped?

ssp245.2060.bio.9 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.10 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.11 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.12 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2021-2060')

ssp245.2060.bio.13 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.14 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.15 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.16 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.17 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.18 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060')

ssp245.2060.bio.19 <- ggplot() + geom_raster(data = cmip6.245.2060.pnw.df, aes(x=x, y=y, fill= wc2_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2041-2060') 

#arrange plots together
temp_plots_245_2060 <- ggarrange(ssp245.2060.bio.1, ssp245.2060.bio.2, ssp245.2060.bio.3, ssp245.2060.bio.4, ssp245.2060.bio.5, ssp245.2060.bio.6, ssp245.2060.bio.7, ssp245.2060.bio.8, ssp245.2060.bio.9, ssp245.2060.bio.10, ssp245.2060.bio.11, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots_245_2060 <- ggarrange(ssp245.2060.bio.12, ssp245.2060.bio.13, ssp245.2060.bio.14, ssp245.2060.bio.15, ssp245.2060.bio.16, ssp245.2060.bio.17, ssp245.2060.bio.18, ssp245.2060.bio.19, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
#save plots
lapply(names(temp_plots_245_2060), function(x)ggsave(filename=paste(x, 'temp_plots_245_2060', '.pdf', sep =''), plot = temp_plots_245_2060[[x]], path = here('figures/world clim/ssp245_2060')))
lapply(names(perc_plots_245_2060), function(x)ggsave(filename=paste(x, 'perc_plots_245_2060', '.pdf', sep =''), plot = perc_plots_245_2060[[x]], path = here('figures/world clim/ssp245_2060')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/ssp245_2060/1temp_plots_245_2060.pdf'), here('figures/world clim/ssp245_2060/2temp_plots_245_2060.pdf'), here('figures/world clim/ssp245_2060/3temp_plots_245_2060.pdf'), 
                      here('figures/world clim/ssp245_2060/1perc_plots_245_2060.pdf'), here('figures/world clim/ssp245_2060/2perc_plots_245_2060.pdf')),
            output = here('figures/world clim/ssp245_2060/wc_245_2060.pdf'))
#####
#CMIP6 SSP245 - 2061-2080
#SSP_245_2080
#####
ssp245.2080.bio.1 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.2 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.3 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.4 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.5 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.6 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.7 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.8 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')
#can see there is a strong effect of mountain ranges - should this be dropped?

ssp245.2080.bio.9 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.10 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.11 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.12 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.13 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.14 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.15 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.16 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.17 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.18 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080')

ssp245.2080.bio.19 <- ggplot() + geom_raster(data = cmip6.245.2080.pnw.df, aes(x=x, y=y, fill= wc2_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp245_2061-2080') 

#arrange plots together
temp_plots_245_2080 <- ggarrange(ssp245.2080.bio.1, ssp245.2080.bio.2, ssp245.2080.bio.3, ssp245.2080.bio.4, ssp245.2080.bio.5, ssp245.2080.bio.6, ssp245.2080.bio.7, ssp245.2080.bio.8, ssp245.2080.bio.9, ssp245.2080.bio.10, ssp245.2080.bio.11, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots_245_2080 <- ggarrange(ssp245.2080.bio.12, ssp245.2080.bio.13, ssp245.2080.bio.14, ssp245.2080.bio.15, ssp245.2080.bio.16, ssp245.2080.bio.17, ssp245.2080.bio.18, ssp245.2080.bio.19, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
#save plots
lapply(names(temp_plots_245_2080), function(x)ggsave(filename=paste(x, 'temp_plots_245_2080', '.pdf', sep =''), plot = temp_plots_245_2080[[x]], path = here('figures/world clim/ssp245_2080')))
lapply(names(perc_plots_245_2080), function(x)ggsave(filename=paste(x, 'perc_plots_245_2080', '.pdf', sep =''), plot = perc_plots_245_2080[[x]], path = here('figures/world clim/ssp245_2080')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/ssp245_2080/1temp_plots_245_2080.pdf'), here('figures/world clim/ssp245_2080/2temp_plots_245_2080.pdf'), here('figures/world clim/ssp245_2080/3temp_plots_245_2080.pdf'), 
                      here('figures/world clim/ssp245_2080/1perc_plots_245_2080.pdf'), here('figures/world clim/ssp245_2080/2perc_plots_245_2080.pdf')),
            output = here('figures/world clim/ssp245_2080/wc_245_2080.pdf'))
#####
#CMIP6 SSP585 - 2021-2040
#SSP_245_2040
#####
ssp585.2040.bio.1 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.2 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.3 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.4 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.5 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.6 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.7 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.8 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')
#can see there is a strong effect of mountain ranges - should this be dropped?

ssp585.2040.bio.9 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.10 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.11 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.12 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.13 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.14 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.15 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.16 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.17 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.18 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040')

ssp585.2040.bio.19 <- ggplot() + geom_raster(data = cmip6.585.2040.pnw.df, aes(x=x, y=y, fill= wc2_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2040') 

#arrange plots together
temp_plots_585_2040 <- ggarrange(ssp585.2040.bio.1, ssp585.2040.bio.2, ssp585.2040.bio.3, ssp585.2040.bio.4, ssp585.2040.bio.5, ssp585.2040.bio.6, ssp585.2040.bio.7, ssp585.2040.bio.8, ssp585.2040.bio.9, ssp585.2040.bio.10, ssp585.2040.bio.11, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots_585_2040 <- ggarrange(ssp585.2040.bio.12, ssp585.2040.bio.13, ssp585.2040.bio.14, ssp585.2040.bio.15, ssp585.2040.bio.16, ssp585.2040.bio.17, ssp585.2040.bio.18, ssp585.2040.bio.19, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
#save plots
lapply(names(temp_plots_585_2040), function(x)ggsave(filename=paste(x, 'temp_plots_585_2040', '.pdf', sep =''), plot = temp_plots_585_2040[[x]], path = here('figures/world clim/ssp585_2040')))
lapply(names(perc_plots_585_2040), function(x)ggsave(filename=paste(x, 'perc_plots_585_2040', '.pdf', sep =''), plot = perc_plots_585_2040[[x]], path = here('figures/world clim/ssp585_2040')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/ssp585_2040/1temp_plots_585_2040.pdf'), here('figures/world clim/ssp585_2040/2temp_plots_585_2040.pdf'), here('figures/world clim/ssp585_2040/3temp_plots_585_2040.pdf'), 
                      here('figures/world clim/ssp585_2040/1perc_plots_585_2040.pdf'), here('figures/world clim/ssp585_2040/2perc_plots_585_2040.pdf')),
            output = here('figures/world clim/ssp585_2040/wc_585_2040.pdf'))
#####
#CMIP6 SSP585 - 2041-2060
#SSP_585_2060
#####
ssp585.2060.bio.1 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.2 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.3 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.4 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.5 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.6 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.7 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.8 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')
#can see there is a strong effect of mountain ranges - should this be dropped?

ssp585.2060.bio.9 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.10 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.11 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.12 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2021-2060')

ssp585.2060.bio.13 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.14 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.15 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.16 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.17 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.18 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060')

ssp585.2060.bio.19 <- ggplot() + geom_raster(data = cmip6.585.2060.pnw.df, aes(x=x, y=y, fill= wc2_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2041-2060') 

#arrange plots together
temp_plots_585_2060 <- ggarrange(ssp585.2060.bio.1, ssp585.2060.bio.2, ssp585.2060.bio.3, ssp585.2060.bio.4, ssp585.2060.bio.5, ssp585.2060.bio.6, ssp585.2060.bio.7, ssp585.2060.bio.8, ssp585.2060.bio.9, ssp585.2060.bio.10, ssp585.2060.bio.11, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots_585_2060 <- ggarrange(ssp585.2060.bio.12, ssp585.2060.bio.13, ssp585.2060.bio.14, ssp585.2060.bio.15, ssp585.2060.bio.16, ssp585.2060.bio.17, ssp585.2060.bio.18, ssp585.2060.bio.19, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
#save plots
lapply(names(temp_plots_585_2060), function(x)ggsave(filename=paste(x, 'temp_plots_585_2060', '.pdf', sep =''), plot = temp_plots_585_2060[[x]], path = here('figures/world clim/ssp585_2060')))
lapply(names(perc_plots_585_2060), function(x)ggsave(filename=paste(x, 'perc_plots_585_2060', '.pdf', sep =''), plot = perc_plots_585_2060[[x]], path = here('figures/world clim/ssp585_2060')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/ssp585_2060/1temp_plots_585_2060.pdf'), here('figures/world clim/ssp585_2060/2temp_plots_585_2060.pdf'), here('figures/world clim/ssp585_2060/3temp_plots_585_2060.pdf'), 
                      here('figures/world clim/ssp585_2060/1perc_plots_585_2060.pdf'), here('figures/world clim/ssp585_2060/2perc_plots_585_2060.pdf')),
            output = here('figures/world clim/ssp585_2060/wc_585_2060.pdf'))
#####
#CMIP6 SSP585 - 2061-2080
#SSP_585_2080
#####
ssp585.2080.bio.1 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_1)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -120, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129, y = 59, label = 'Annual Mean Temperature', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129, y = 59, label = 'Annual Mean Temperature'), fontface = 'bold', size = 2, color = 'black') + #this does not resize when printing to pdf
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.2 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_2)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138, xmax = -124, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Mean Diurnal Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Mean Diurnal Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.3 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_3)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -135.5, xmax = -126.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -131, y = 59, label = 'Isothermality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Isothermality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.4 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_4)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C x 100') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Temperature Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Temperature Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.5 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_5)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.7, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom= 'text',x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Max Temperature of\nWarmest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.6 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_6)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -124.5, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Min Temperature of\nColdest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.7 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_7)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.6, xmax = -120.4, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -129.5, y = 59, label = 'Temperature Annual Range', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -129.5, y = 59, label = 'Temperature Annual Range'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.8 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_8)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')
#can see there is a strong effect of mountain ranges - should this be dropped?

ssp585.2080.bio.9 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_9)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.10 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_10)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.11 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_11)) +
  scale_fill_viridis(option='A') +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='°C') +
  geom_rect(aes(xmin = -138.4, xmax = -123.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Mean Temperature of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.12 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_12)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.2, xmax = -123.6, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text',x = -131, y = 59, label = 'Annual Precipitation', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -131, y = 59, label = 'Annual Precipitation'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.13 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_13)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.14 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_14)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Month', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Month'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.15 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_15)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.5, xmax = -121.5, ymin = 58.4, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -130, y = 59, label = 'Precipitation Seasonality', fontface = 'bold', size = 2, color = 'black') +
  #geom_text(aes(x = -130, y = 59, label = 'Precipitation Seasonality'), fontface = 'bold', size = 2, color = 'black') +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.16 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_16)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWettest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.17 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_17)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -127, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nDriest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.18 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_18)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nWarmest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080')

ssp585.2080.bio.19 <- ggplot() + geom_raster(data = cmip6.585.2080.pnw.df, aes(x=x, y=y, fill= wc2_19)) +
  scale_fill_viridis(option='D', direction = -1) +
  geom_sf(data=pnw.bound, fill = NA, colour = 'black', size = 1) +
  #geom_point(data=aa.occ.thin, aes(x=decimalLongitude,y=decimalLatitude), size = .5, shape=4, colour = 'grey65') + #plot thinned occurrence data
  labs(fill='mm') +
  geom_rect(aes(xmin = -138.4, xmax = -126.6, ymin = 57.1, ymax = 59.6), colour = 'black', fill = 'white', linewidth = 1, alpha = .2) +
  annotate(geom = 'text', x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter', fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  #geom_text(aes(x = -138, y = 58.5, label = 'Precipitation of\nColdest Quarter'), fontface = 'bold', size = 2, color = 'black', lineheight = 1, hjust = 0) +
  scale_y_continuous(expand = c(.005,.005)) +
  scale_x_continuous(expand = c(.005,.005)) +
  theme_maps + ggtitle('ssp585_2061-2080') 

#arrange plots together
temp_plots_585_2080 <- ggarrange(ssp585.2080.bio.1, ssp585.2080.bio.2, ssp585.2080.bio.3, ssp585.2080.bio.4, ssp585.2080.bio.5, ssp585.2080.bio.6, ssp585.2080.bio.7, ssp585.2080.bio.8, ssp585.2080.bio.9, ssp585.2080.bio.10, ssp585.2080.bio.11, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))

perc_plots_585_2080 <- ggarrange(ssp585.2080.bio.12, ssp585.2080.bio.13, ssp585.2080.bio.14, ssp585.2080.bio.15, ssp585.2080.bio.16, ssp585.2080.bio.17, ssp585.2080.bio.18, ssp585.2080.bio.19, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
#save plots
lapply(names(temp_plots_585_2080), function(x)ggsave(filename=paste(x, 'temp_plots_585_2080', '.pdf', sep =''), plot = temp_plots_585_2080[[x]], path = here('figures/world clim/ssp585_2080')))
lapply(names(perc_plots_585_2080), function(x)ggsave(filename=paste(x, 'perc_plots_585_2080', '.pdf', sep =''), plot = perc_plots_585_2080[[x]], path = here('figures/world clim/ssp585_2080')))

#combine temp and perc pdfs
pdf_combine(input = c(here('figures/world clim/ssp585_2080/1temp_plots_585_2080.pdf'), here('figures/world clim/ssp585_2080/2temp_plots_585_2080.pdf'), here('figures/world clim/ssp585_2080/3temp_plots_585_2080.pdf'), 
                      here('figures/world clim/ssp585_2080/1perc_plots_585_2080.pdf'), here('figures/world clim/ssp585_2080/2perc_plots_585_2080.pdf')),
            output = here('figures/world clim/ssp585_2080/wc_585_2080.pdf'))
#####
#group bio clim vars
blank.plot <- ggplot() + geom_blank() +theme_transparent() #make blank plot for ggarange
#BIO1 - Mean Annual Temp
bio.1.plots <- ggarrange(bio.1, blank.plot, ssp245.2040.bio.1,ssp585.2040.bio.1, ssp245.2060.bio.1, ssp585.2060.bio.1, ssp245.2080.bio.1, ssp585.2080.bio.1, 
                                 ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.1.plots), function(x)ggsave(filename=paste(x, 'bio.1.plots', '.pdf', sep =''), plot = bio.1.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.1.plots.pdf'), here('figures/world clim/wc_combined/2bio.1.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio1_combined.pdf'))
#BIO2 - Mean Diurnal Range
bio.2.plots <- ggarrange(bio.2, blank.plot, ssp245.2040.bio.2,ssp585.2040.bio.2, ssp245.2060.bio.2, ssp585.2060.bio.2, ssp245.2080.bio.2, ssp585.2080.bio.2, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.2.plots), function(x)ggsave(filename=paste(x, 'bio.2.plots', '.pdf', sep =''), plot = bio.2.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.2.plots.pdf'), here('figures/world clim/wc_combined/2bio.2.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio2_combined.pdf'))
#BIO3 - Isothermality 
bio.3.plots <- ggarrange(bio.3, blank.plot, ssp245.2040.bio.3,ssp585.2040.bio.3, ssp245.2060.bio.3, ssp585.2060.bio.3, ssp245.2080.bio.3, ssp585.2080.bio.3, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.3.plots), function(x)ggsave(filename=paste(x, 'bio.3.plots', '.pdf', sep =''), plot = bio.3.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.3.plots.pdf'), here('figures/world clim/wc_combined/2bio.3.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio3_combined.pdf'))
#BIO4 - Temperature Seasonality (standard deviation ×100)
bio.4.plots <- ggarrange(bio.4, blank.plot, ssp245.2040.bio.4,ssp585.2040.bio.4, ssp245.2060.bio.4, ssp585.2060.bio.4, ssp245.2080.bio.4, ssp585.2080.bio.4, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.4.plots), function(x)ggsave(filename=paste(x, 'bio.4.plots', '.pdf', sep =''), plot = bio.4.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.4.plots.pdf'), here('figures/world clim/wc_combined/2bio.4.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio4_combined.pdf'))
#BIO5 - Max Temperature of Warmest Month
bio.5.plots <- ggarrange(bio.5, blank.plot, ssp245.2040.bio.5,ssp585.2040.bio.5, ssp245.2060.bio.5, ssp585.2060.bio.5, ssp245.2080.bio.5, ssp585.2080.bio.5, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.5.plots), function(x)ggsave(filename=paste(x, 'bio.5.plots', '.pdf', sep =''), plot = bio.5.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.5.plots.pdf'), here('figures/world clim/wc_combined/2bio.5.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio5_combined.pdf'))
#BIO6 - Min Temperature of Coldest Month
bio.6.plots <- ggarrange(bio.6, blank.plot, ssp245.2040.bio.6,ssp585.2040.bio.6, ssp245.2060.bio.6, ssp585.2060.bio.6, ssp245.2080.bio.6, ssp585.2080.bio.6, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.6.plots), function(x)ggsave(filename=paste(x, 'bio.6.plots', '.pdf', sep =''), plot = bio.6.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.6.plots.pdf'), here('figures/world clim/wc_combined/2bio.6.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio6_combined.pdf'))
#BIO7 - Temperature Annual Range 
bio.7.plots <- ggarrange(bio.7, blank.plot, ssp245.2040.bio.7,ssp585.2040.bio.7, ssp245.2060.bio.7, ssp585.2060.bio.7, ssp245.2080.bio.7, ssp585.2080.bio.7, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.7.plots), function(x)ggsave(filename=paste(x, 'bio.7.plots', '.pdf', sep =''), plot = bio.7.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.7.plots.pdf'), here('figures/world clim/wc_combined/2bio.7.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio7_combined.pdf'))
#BIO8 - Mean Temperature of Wettest Quarter
bio.8.plots <- ggarrange(bio.8, blank.plot, ssp245.2040.bio.8,ssp585.2040.bio.8, ssp245.2060.bio.8, ssp585.2060.bio.8, ssp245.2080.bio.8, ssp585.2080.bio.8, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.8.plots), function(x)ggsave(filename=paste(x, 'bio.8.plots', '.pdf', sep =''), plot = bio.8.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.8.plots.pdf'), here('figures/world clim/wc_combined/2bio.8.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio8_combined.pdf'))
#BIO9 - Mean Temperature of Driest Quarter
bio.9.plots <- ggarrange(bio.9, blank.plot, ssp245.2040.bio.9,ssp585.2040.bio.9, ssp245.2060.bio.9, ssp585.2060.bio.9, ssp245.2080.bio.9, ssp585.2080.bio.9, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.9.plots), function(x)ggsave(filename=paste(x, 'bio.9.plots', '.pdf', sep =''), plot = bio.9.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.9.plots.pdf'), here('figures/world clim/wc_combined/2bio.9.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio9_combined.pdf'))
#BIO10 - Mean Temperature of Warmest Quarter
bio.10.plots <- ggarrange(bio.10, blank.plot, ssp245.2040.bio.10,ssp585.2040.bio.10, ssp245.2060.bio.10, ssp585.2060.bio.10, ssp245.2080.bio.10, ssp585.2080.bio.10, 
                         ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.10.plots), function(x)ggsave(filename=paste(x, 'bio.10.plots', '.pdf', sep =''), plot = bio.10.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.10.plots.pdf'), here('figures/world clim/wc_combined/2bio.10.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio10_combined.pdf'))
#BIO11 - Mean Temperature of Coldest Quarter
bio.11.plots <- ggarrange(bio.11, blank.plot, ssp245.2040.bio.11,ssp585.2040.bio.11, ssp245.2060.bio.11, ssp585.2060.bio.11, ssp245.2080.bio.11, ssp585.2080.bio.11, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.11.plots), function(x)ggsave(filename=paste(x, 'bio.11.plots', '.pdf', sep =''), plot = bio.11.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.11.plots.pdf'), here('figures/world clim/wc_combined/2bio.11.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio11_combined.pdf'))
#BIO12 - Annual Precipitation
bio.12.plots <- ggarrange(bio.12, blank.plot, ssp245.2040.bio.12,ssp585.2040.bio.12, ssp245.2060.bio.12, ssp585.2060.bio.12, ssp245.2080.bio.12, ssp585.2080.bio.12, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.12.plots), function(x)ggsave(filename=paste(x, 'bio.12.plots', '.pdf', sep =''), plot = bio.12.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.12.plots.pdf'), here('figures/world clim/wc_combined/2bio.12.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio12_combined.pdf'))
#BIO13 - Precipitation of Wettest Month
bio.13.plots <- ggarrange(bio.13, blank.plot, ssp245.2040.bio.13,ssp585.2040.bio.13, ssp245.2060.bio.13, ssp585.2060.bio.13, ssp245.2080.bio.13, ssp585.2080.bio.13, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.13.plots), function(x)ggsave(filename=paste(x, 'bio.13.plots', '.pdf', sep =''), plot = bio.13.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.13.plots.pdf'), here('figures/world clim/wc_combined/2bio.13.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio13_combined.pdf'))
#BIO14 - Precipitation of Driest Month
bio.14.plots <- ggarrange(bio.14, blank.plot, ssp245.2040.bio.14,ssp585.2040.bio.14, ssp245.2060.bio.14, ssp585.2060.bio.14, ssp245.2080.bio.14, ssp585.2080.bio.14, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.14.plots), function(x)ggsave(filename=paste(x, 'bio.14.plots', '.pdf', sep =''), plot = bio.14.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.14.plots.pdf'), here('figures/world clim/wc_combined/2bio.14.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio14_combined.pdf'))
#BIO15 - Precipitation Seasonality (Coefficient of Variation)
bio.15.plots <- ggarrange(bio.15, blank.plot, ssp245.2040.bio.15,ssp585.2040.bio.15, ssp245.2060.bio.15, ssp585.2060.bio.15, ssp245.2080.bio.15, ssp585.2080.bio.15, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.15.plots), function(x)ggsave(filename=paste(x, 'bio.15.plots', '.pdf', sep =''), plot = bio.15.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.15.plots.pdf'), here('figures/world clim/wc_combined/2bio.15.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio15_combined.pdf'))
#BIO16 - Precipitation of Wettest Quarter
bio.16.plots <- ggarrange(bio.16, blank.plot, ssp245.2040.bio.16,ssp585.2040.bio.16, ssp245.2060.bio.16, ssp585.2060.bio.16, ssp245.2080.bio.16, ssp585.2080.bio.16, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.16.plots), function(x)ggsave(filename=paste(x, 'bio.16.plots', '.pdf', sep =''), plot = bio.16.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.16.plots.pdf'), here('figures/world clim/wc_combined/2bio.16.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio16_combined.pdf'))
#BIO17 - Precipitation of Driest Quarter
bio.17.plots <- ggarrange(bio.17, blank.plot, ssp245.2040.bio.17,ssp585.2040.bio.17, ssp245.2060.bio.17, ssp585.2060.bio.17, ssp245.2080.bio.17, ssp585.2080.bio.17, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.17.plots), function(x)ggsave(filename=paste(x, 'bio.17.plots', '.pdf', sep =''), plot = bio.17.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.17.plots.pdf'), here('figures/world clim/wc_combined/2bio.17.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio17_combined.pdf'))
#BIO18 - Precipitation of Warmest Quarter
bio.18.plots <- ggarrange(bio.18, blank.plot, ssp245.2040.bio.18,ssp585.2040.bio.18, ssp245.2060.bio.18, ssp585.2060.bio.18, ssp245.2080.bio.18, ssp585.2080.bio.18, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.18.plots), function(x)ggsave(filename=paste(x, 'bio.18.plots', '.pdf', sep =''), plot = bio.18.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.18.plots.pdf'), here('figures/world clim/wc_combined/2bio.18.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio18_combined.pdf'))
#BIO19 - Precipitation of Coldest Quarter
bio.19.plots <- ggarrange(bio.19, blank.plot, ssp245.2040.bio.19,ssp585.2040.bio.19, ssp245.2060.bio.19, ssp585.2060.bio.19, ssp245.2080.bio.19, ssp585.2080.bio.19, 
                          ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1))
lapply(names(bio.19.plots), function(x)ggsave(filename=paste(x, 'bio.19.plots', '.pdf', sep =''), plot = bio.19.plots[[x]], path = here('figures/world clim/wc_combined')))
pdf_combine(input = c(here('figures/world clim/wc_combined/1bio.19.plots.pdf'), here('figures/world clim/wc_combined/2bio.19.plots.pdf')),
            output = here('figures/world clim/wc_combined/wc_bio19_combined.pdf'))
