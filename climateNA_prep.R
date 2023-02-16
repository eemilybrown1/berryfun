########
#Terrell Roulston
#terrell.roulston@ubc.ca
#Started Feb 15, 2023
#ZOOL/RES509
#Berry Bunch Project
#Preparing for ClimateNA steps
#Pulling geolocation data for each species occurence from gbif data and adding elevation
#Export ClimateNA input ready .csv files for climatic variables
########
library(readr) #download cvs
library(dplyr) #piping and data managment
library(elevatr) #elevation data from the web
library(sp)
library(sf)
library(rgdal)

#ignore parsing error warning, it imports fine
amel.alni.df <- read_csv("data/amelanchieralnifolia.csv",  show_col_types = F)%>%  #import dataset, ask not to show coloumn types
                    filter(stateProvince == c('British Columbia', 'Washington', 'Oregon')) %>% #filter data found within PNW
                    select(gbifID, species, decimalLongitude, decimalLatitude, elevation) #note that some occurrences already contain elevation data
          
#note: gbif occurrences elevation units is in meters (m)
amel.alni.location <- amel.alni.df %>% 
                    filter(is.na(elevation)) %>% #filter by elevation values containing NA (i.e. missing elevation)
                    select(x = decimalLongitude, y = decimalLatitude) %>% 
                    as.data.frame()


#pull elevation data
wkt <- sf::st_crs(4326)[[2]] #projection for lat-long coordinates - WGS84 ellipsoid


#data source 'src' is Copernicus Digital Elevation Model (DEM) - an AMAZON Sustainability Data Initiative
#z = 12, determines ground resolution, at 45 degrees lat z = 10 is ~108 meters (https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution)
amel.alni.elev <- get_elev_point(locations = amel.alni.location, units = 'meters', prj = wkt, src = 'aws', z = 10)
save(amel.alni.elve, file='amel.alni.elev.Rda')

