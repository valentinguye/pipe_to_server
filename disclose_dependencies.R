########
# This script discloses the project required packages to renv::dependencies
########

# NOT RUN 
library(data.table)
library(dplyr)
library(plyr)
library(Hmisc)
library(sjmisc)
library(stringr)
library(readstata13)
library(foreign)
library(readxl)
library(writexl)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(gfcanalysis)
library(doParallel)
library(foreach)
library(parallel)

# RUN if a new line 'library(new.package)' has been added 
renv::snapshot()



#################################################################################################

# Specify packages to install
# packages <- c("data.table", "dplyr", "plyr", "Hmisc", "sjmisc", "stringr",
#               "readstata13", "foreign", "readxl", "writexl", 
#               "raster", "rgdal", "sp", "sf","gfcanalysis",
#               "doParallel", "foreach", "parallel")


