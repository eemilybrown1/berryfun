library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(dplyr) # data management
library(lubridate) # round months to years
library(weathermetrics) # convert kelvin to celcius
library(geodata)
library(here)
library(terra)
library(tidyr)
library(tidyverse)
#basemap of pnw
bound <- c('Oregon', 'Washington', 'British Columbia')
pnw.bound <- gadm(country=c('USA', 'CAN'), level =1, path = here('data')) 
pnw.bound <- pnw.bound[pnw.bound$NAME_1 %in% bound,]
pnw.bound.exnt <- ext(pnw.bound)


#SSP data from https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0
#historical data from https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/MRI/MRI-ESM2-0/historical
#Surface temperature----
historical_nc <- nc_open("C:\\Users\\terre\\Documents\\R\\berryfun\\data\\MRI-EMS2-0\\tas_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc")


historical_lon <- ncvar_get(historical_nc, 'lon')
historical_lat <- ncvar_get(historical_nc, "lat", verbose = F)
historical_tim <- ncvar_get(historical_nc, "time")
historical_tas <- ncvar_get(historical_nc, "tas")

historical_brick <- brick(historical_tas, xmn=min(historical_lat), xmx=max(historical_lat), ymn=min(historical_lon), ymx=max(historical_lon), crs="EPSG:4326")
historical_brick <- flip(t(historical_brick), direction='y')

nc_close(historical_nc)


historical_series <- raster::as.data.frame(historical_brick, xy=T)

historical_series_pnw <- historical_series %>% dplyr::filter(between(x,221, 246)) %>% 
  dplyr::filter(between(y, 42, 60)) %>% 
  select(-c(1:2)) 


historical_series_mean <- data.frame(mean_k = apply(historical_series_pnw, 2, mean), 
                                 year_month = seq(as.Date("1850/1/1"), by = 'month', length.out = 1980)) %>%   
  mutate(mean_c = kelvin.to.celsius(mean_k)) %>% 
  group_by(year = floor_date(year_month, 'year')) %>% 
  summarise(tas_mean_c=mean(mean_c),
            tas_sd_c=sd(mean_c),
            scenario=c(paste('Historic')))


ssp245_nc <- nc_open("C:\\Users\\terre\\Documents\\R\\berryfun\\data\\MRI-EMS2-0\\tas_Amon_MRI-ESM2-0_ssp245_r1i1p1f1_gn_201501-210012.nc")


ssp245_lon <- ncvar_get(ssp245_nc, 'lon')
ssp245_lat <- ncvar_get(ssp245_nc, "lat", verbose = F)
ssp245_tim <- ncvar_get(ssp245_nc, "time")
ssp245_tas <- ncvar_get(ssp245_nc, "tas")

ssp245_brick <- brick(ssp245_tas, xmn=min(ssp245_lat), xmx=max(ssp245_lat), ymn=min(ssp245_lon), ymx=max(ssp245_lon), crs="EPSG:4326")
ssp245_brick <- flip(t(ssp245_brick), direction='y')

nc_close(ssp245_nc)


ssp245_series <- raster::as.data.frame(ssp245_brick, xy=T)

ssp245_series_pnw <- ssp245_series %>% dplyr::filter(between(x,221, 246)) %>% 
  dplyr::filter(between(y, 42, 60)) %>% 
  select(-c(1:2)) 


ssp245_series_mean <- data.frame(mean_k = apply(ssp245_series_pnw, 2, mean), 
                                 year_month = seq(as.Date("2015/1/1"), by = 'month', length.out = 1032)) %>%   
  mutate(mean_c = kelvin.to.celsius(mean_k)) %>% 
  group_by(year = floor_date(year_month, 'year')) %>% 
  summarise(tas_mean_c=mean(mean_c),
            tas_sd_c=sd(mean_c),
            scenario=c(paste('SSP245')))

ssp585_nc <- nc_open("C:\\Users\\terre\\Documents\\R\\berryfun\\data\\MRI-EMS2-0\\tas_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_201501-210012.nc")

ssp585_lon <- ncvar_get(ssp585_nc, 'lon')
ssp585_lat <- ncvar_get(ssp585_nc, "lat", verbose = F)
ssp585_tim <- ncvar_get(ssp585_nc, "time")
ssp585_tas <- ncvar_get(ssp585_nc, "tas")

ssp585_brick <- brick(ssp585_tas, xmn=min(ssp585_lat), xmx=max(ssp585_lat), ymn=min(ssp585_lon), ymx=max(ssp585_lon), crs="EPSG:4326")
ssp585_brick <- flip(t(ssp585_brick), direction='y')

nc_close(ssp585_nc)

plot(ssp585_brick$layer.1032)

ssp585_series <- raster::as.data.frame(ssp585_brick, xy=T)

ssp585_series_pnw <- ssp585_series %>% dplyr::filter(between(x,221, 246)) %>% 
                                       dplyr::filter(between(y, 42, 60)) %>% 
                                       select(-c(1:2)) 

                    
ssp585_series_mean <- data.frame(mean_k = apply(ssp585_series_pnw, 2, mean), 
                                   year_month = seq(as.Date("2015/1/1"), by = 'month', length.out = 1032)) %>%   
                        mutate(mean_c = kelvin.to.celsius(mean_k)) %>% 
                        group_by(year = floor_date(year_month, 'year')) %>% 
                        summarise(tas_mean_c=mean(mean_c),
                                  tas_sd_c=sd(mean_c),
                                  scenario=c(paste('SSP585')))
                                                      


ssp_tas_data <- rbind(historical_series_mean, ssp245_series_mean, ssp585_series_mean) 


#plot time series of mean annual surface temperature 
ggplot() +
    geom_line(data=ssp_tas_data, aes(x=year, y=tas_mean_c, colour = scenario), linewidth = 1, alpha = 0.1) +
    geom_smooth(data=ssp_tas_data, aes(x=year, y=tas_mean_c, colour = scenario), linewidth = 1, alpha = 0.3, level = 0.95, se = T, method = 'loess') +
    scale_x_date(date_breaks = '20 years', date_labels = "%Y", limits = as.Date(c("1850-01-01",'2100-01-01')), expand = c(0,0)) +
    scale_y_continuous(n.breaks = 10) +
    theme_classic() +
    theme(plot.margin = unit(c(0.2,1,0.5,0.5), 'cm')) +
    labs(x = NULL, y = 'Pacific Northwest\nMean Annual Surface Tempature (°C)', colour = 'MRI-ESM2-0 Projection') +
    scale_color_manual(values =  c('black','blue', 'red')) +
    theme(axis.text = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 14, colour = 'black'),
          legend.text = element_text(size = 12, colour = 'black'),
          legend.title = element_text(size =12, colour = 'black'),
          legend.position = c(0.175, 0.8),
          legend.background = element_rect(colour = 'black'), 
          legend.key.width = unit(1.5, 'cm')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))


#precipitation----
#https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/pr/gn
historical_pr_nc <- nc_open("C:\\Users\\terre\\Documents\\R\\berryfun\\data\\MRI-EMS2-0\\pr_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc")

historical_pr_lon <- ncvar_get(historical_pr_nc, 'lon')
historical_pr_lat <- ncvar_get(historical_pr_nc, "lat", verbose = F)
historical_pr_tim <- ncvar_get(historical_pr_nc, "time")
historical_pr_pr <- ncvar_get(historical_pr_nc, "pr")

historical_pr_brick <- brick(historical_pr_pr, xmn=min(historical_pr_lat), xmx=max(historical_pr_lat), ymn=min(historical_pr_lon), ymx=max(historical_pr_lon), crs="EPSG:4326")
historical_pr_brick <- flip(t(historical_pr_brick), direction='y')

nc_close(historical_pr_nc)

plot(historical_pr_brick$layer.1) #take a peak

historical_pr_series <- raster::as.data.frame(historical_pr_brick, xy=T) 

historical_pr_series_pnw <- historical_pr_series %>% dplyr::filter(between(x,221, 246)) %>% 
  dplyr::filter(between(y, 42, 60)) %>% 
  select(-c(1:2)) #deselect lat and lon cols


historical_pr_series_mean <- data.frame(mean_kgms = apply(historical_pr_series_pnw, 2, mean), 
                                    year_month = seq(as.Date("1850/1/1"), by = 'month', length.out = 1980)) %>%   
  mutate(mean_mm = (mean_kgms*86400*30)) %>% #covert kg*m^2*s^-1 to mm/month by multiplying by 86400 * 30 
  group_by(year = floor_date(year_month, 'year')) %>% 
  summarise(pr_mean_mm=sum(mean_mm),
            pr_sd_mm=sd(mean_mm),
            scenario=c(paste('Historic')))

#SSP245 perc
ssp245_pr_nc <- nc_open("C:\\Users\\terre\\Documents\\R\\berryfun\\data\\MRI-EMS2-0\\pr_Amon_MRI-ESM2-0_ssp245_r1i1p1f1_gn_201501-210012.nc")

ssp245_pr_lon <- ncvar_get(ssp245_pr_nc, 'lon')
ssp245_pr_lat <- ncvar_get(ssp245_pr_nc, "lat", verbose = F)
ssp245_pr_tim <- ncvar_get(ssp245_pr_nc, "time")
ssp245_pr_pr <- ncvar_get(ssp245_pr_nc, "pr")

ssp245_pr_brick <- brick(ssp245_pr_pr, xmn=min(ssp245_pr_lat), xmx=max(ssp245_pr_lat), ymn=min(ssp245_pr_lon), ymx=max(ssp245_pr_lon), crs="EPSG:4326")
ssp245_pr_brick <- flip(t(ssp245_pr_brick), direction='y')

nc_close(ssp245_pr_nc)

plot(ssp245_pr_brick$layer.1) #take a peak

ssp245_pr_series <- raster::as.data.frame(ssp245_pr_brick, xy=T) 

ssp245_pr_series_pnw <- ssp245_pr_series %>% dplyr::filter(between(x,221, 246)) %>% 
  dplyr::filter(between(y, 42, 60)) %>% 
  select(-c(1:2)) #deselect lat and lon cols


ssp245_pr_series_mean <- data.frame(mean_kgms = apply(ssp245_pr_series_pnw, 2, mean), 
                                    year_month = seq(as.Date("2015/1/1"), by = 'month', length.out = 1032)) %>%   
  mutate(mean_mm = (mean_kgms*86400*30)) %>% #covert kg*m^2*s^-1 to mm/day by multiplying by 86400 * 30 
  group_by(year = floor_date(year_month, 'year')) %>% 
  summarise(pr_mean_mm=sum(mean_mm),
            pr_sd_mm=sd(mean_mm),
            scenario=c(paste('SSP245')))

ssp585_pr_nc <- nc_open("C:\\Users\\terre\\Documents\\R\\berryfun\\data\\MRI-EMS2-0\\pr_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_201501-210012.nc")

ssp585_pr_lon <- ncvar_get(ssp585_pr_nc, 'lon')
ssp585_pr_lat <- ncvar_get(ssp585_pr_nc, "lat", verbose = F)
ssp585_pr_tim <- ncvar_get(ssp585_pr_nc, "time")
ssp585_pr_pr <- ncvar_get(ssp585_pr_nc, "pr")

ssp585_pr_brick <- brick(ssp585_pr_pr, xmn=min(ssp585_pr_lat), xmx=max(ssp585_pr_lat), ymn=min(ssp585_pr_lon), ymx=max(ssp585_pr_lon), crs="EPSG:4326")
ssp585_pr_brick <- flip(t(ssp585_pr_brick), direction='y')

nc_close(ssp585_pr_nc)

plot(ssp585_pr_brick$layer.1) #take a peak

ssp585_pr_series <- raster::as.data.frame(ssp585_pr_brick, xy=T) 

ssp585_pr_series_pnw <- ssp585_pr_series %>% dplyr::filter(between(x,221, 246)) %>% 
  dplyr::filter(between(y, 42, 60)) %>% 
  select(-c(1:2)) #deselect lat and lon cols


ssp585_pr_series_mean <- data.frame(mean_kgms = apply(ssp585_pr_series_pnw, 2, mean), 
                                 year_month = seq(as.Date("2015/1/1"), by = 'month', length.out = 1032)) %>%   
  mutate(mean_mm = (mean_kgms*86400*30)) %>% #covert kg*m^2*s^-1 to mm/day by multiplying by 86400 * 30 
  group_by(year = floor_date(year_month, 'year')) %>% 
  summarise(pr_mean_mm=sum(mean_mm),
            pr_sd_mm=sd(mean_mm),
            scenario=c(paste('SSP585')))

ssp_pr_data <- rbind(historical_pr_series_mean, ssp245_pr_series_mean, ssp585_pr_series_mean) 


#plot time series of mean annual surface temperature 
ggplot() +
  geom_line(data=ssp_pr_data, aes(x=year, y=pr_mean_mm, colour = scenario), linewidth = 1, alpha = 0.1) +
  geom_smooth(data=ssp_pr_data, aes(x=year, y=pr_mean_mm, colour = scenario), linewidth = 1, alpha = 0.3, level = 0.95, se = T, method = 'loess') +
  scale_x_date(date_breaks = '20 years', date_labels = "%Y", limits = as.Date(c("1850-01-01",'2100-01-01')), expand = c(0,0)) +
  scale_y_continuous(n.breaks = 10) +
  theme_classic() +
  theme(plot.margin = unit(c(0.2,1,0.5,0.5), 'cm')) +
  labs(x = NULL, y = 'Pacific Northwest\nMean Annual Precipitation (mm)', colour = 'MRI-ESM2-0 Projection') +
  scale_color_manual(values =  c('black','blue', 'red')) +
  theme(axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        legend.text = element_text(size = 12, colour = 'black'),
        legend.title = element_text(size =12, colour = 'black'),
        legend.position = c(0.175, 0.8),
        legend.background = element_rect(colour = 'black'), 
        legend.key.width = unit(1.5, 'cm')) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

