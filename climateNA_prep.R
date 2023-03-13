########
#Terrell Roulston
#terrell.roulston@ubc.ca
#Started Feb 15, 2023
#ZOOL/RES509
#Berry Bunch Project
#Preparing for ClimateNA steps
#Pulling geolocation data for each species occurrence from gbif data and adding elevation
#Export ClimateNA input ready .csv files for climatic variables
########
library(readr) #download cvs
library(dplyr) #piping and data management
library(elevatr) #elevation data from the web
library(sf) #geo referencing package
library(rgdal) #map projection package

#occurrence  
#####Repeat the following process for all desired species occurrences
#ignore parsing error warning, it imports fine
occurrence.df <- read_csv("data/vacciniumparvifolium.csv",  show_col_types = F)%>%  #import occurrence data set of given species, ask not to show coloumn types
                    filter(stateProvince %in% c('British Columbia', 'Washington', 'Oregon')) %>% #filter data found within PNW
                    select(gbifID, species, decimalLongitude, decimalLatitude, elevation) #note that some occurrences already contain elevation data
          
#note: gbif occurrences elevation units is in meters (m)
occurrence.location <- occurrence.df %>% 
                    filter(is.na(elevation)) %>% #filter by elevation values containing NA (i.e. missing elevation)
                    select(x = decimalLongitude, y = decimalLatitude, gbifID = gbifID) %>% 
                    as.data.frame() #ensure it is df for next steps

#define projection 
wkt <- sf::st_crs(4326)[[2]] #projection for lat-long coordinates - WGS84 ellipsoid
#data source 'src' is Copernicus Digital Elevation Model (DEM) - an AMAZON Sustainability Data Initiative
#'z' aka zoom, determines ground resolution, at 45 degrees lat z = 9 is ~216 meters resolution (https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution)
occurrence.elevation <- get_elev_point(locations = occurrence.location, units = 'meters', prj = wkt, src = 'aws', z = 9) #~3 minute processing 

#STOP
save(occurrence.elevation, file='vacc.parv.elev.Rda') #save data frame of downloaded elevation to save time if reloading
#RENAME FILE^^^

elevation.data <- as.data.frame(occurrence.elevation@data) #elevation data with corresponding ID (gbif ID)


#format for climateNA, col1 = ID_factor, col2 = ID_factor, col3 = Latitude,  col4 = Longitude, col5 = Elevation (in m)
combined.df <- occurrence.df %>% full_join(elevation.data, by = 'gbifID') %>% #join by gbifID of occurrence
                mutate(elevation = coalesce(elevation.x, elevation.y)) %>% #coalesce values of elevation from original df and elevation df (temp_df)
                select(gbifID, species, decimalLatitude, decimalLongitude, elevation) #select columns in correct order (drop elevation units)
View(combined.df) #check 
nrow(occurrence.df) #check number of rows are the same for the occurrence df 
nrow(combined.df) #and new combined df
#STOP
write.csv(combined.df, 'C:\\Users\\terre\\Documents\\R\\berryfun\\data\\climateNA_input\\vacc_parv_input.csv', row.names = F) #specify local path and push to git. make sure to change file name as needed                                                                    
#RENAME FILE^^^

#load occurrence elevation data 
#useful for saving time if needed to restart session
load('amel.alni.elev.Rda') #Amelanchier alnifolia
load('rubu.lasi.elev.Rda') #Rubus lasiococcus
load('rubu.niva.elev.Rda') #Rubus nivalis
load('vacc.parv.elev.Rda')#Vaccinium parvifolium

####
#END
####
