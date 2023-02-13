### Data Collection for occurrence data
# jens ulrich; started oct 27, 2022

# Collect occurrence data from GBIF

library(tidyverse)
library(rgbif)

## --------------------------------------------------
# 1) GBIF DOWNLOAD
## --------------------------------------------------


## --------------------------------------------------
# Enter GBIF credentials

# Enter your GBIF user info here before proceeding

# I will delete this before pushing this file to protect my account privacy
# must be rewritten every time running this script

user="" 
pwd="" 
email=""

## --------------------------------------------------

# Enter download filters
taxonKey <- #
basisOfRecord <- c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION') # museum specimens and human obs (i.e. inaturalist)
hasCoordinate <- TRUE # excludes records without coordinates
# years <- seq(2000, 2022, 1) # limit by years?

## --------------------------------------------------
# Download data from GBIF 
# use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  #pred_in("gadm", stateProvince), 
  pred("hasCoordinate", hasCoordinate),
  #pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)

# now get the records from your GBIF downloads page
setwd("./") # set your working directory
occurrences <- occ_download_get(down_code[1], overwrite = TRUE)
occurrences_df <- occ_download_import(occurrences)


# 
# write.csv(occurrences_df, "_.csv")