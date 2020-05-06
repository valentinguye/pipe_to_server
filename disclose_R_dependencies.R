########

# This script discloses the project required packages to renv::dependencies
# And describes the second-best steps to go through if renv::restore() could not install the project library 

###   THIS IS NOT WRITTEN TO BE RUN COMPLETELY  ### 

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
library(rgdal)
library(raster)
library(sp)
library(sf)
library(gfcanalysis)
library(doParallel)
library(foreach)
library(parallel)


# RUN if a line 'library(package)' has been added, or commented out:
renv::snapshot()
renv::restore()


#################################################################################################

#### SECOND-BEST SOLUTION for packages for which installation through renv::restore() fails ####
  # For instance sf https://github.com/r-spatial/sf/issues/921 
  # or magick, as a dependency of raster and rgdal. 

# /!\ Note that these packages won't necessary be installed in the same version as they were on the project initial data processing
#     and therefore reproducibility is not guaranteed in this case. 

### Follow the steps: 
# 1. In the present script, above, comment out the lines disclosing (by library(package)) the packages 
#   causing renv::restore() to fail; and save the changes made to the script. 
# 2. Run command of line 27 of the present script: renv::snapshot() to remove these packages from the renv.lock file. 
# 3. Run command of line 28 of the present script: renv::restore() to check if all remaining packages are installed properly. 
# 3. Report all the packages you had to remove from the renv.lock for renv::restore() not to fail, 
#   in this trouble_packages vector 
    trouble_packages <- c("")
#   and install them if not available
    allPackages    = c(trouble_packages %in% installed.packages()[ , "Package"]) 
    if(!all(allPackages)) {
      missingIDX = which(allPackages == FALSE)
      needed     = trouble_packages[missingIDX]
      lapply(needed, install.packages)
    }
    











# Specify packages to install
# packages <- c("data.table", "dplyr", "plyr", "Hmisc", "sjmisc", "stringr",
#               "readstata13", "foreign", "readxl", "writexl", 
#               "raster", "rgdal", "sp", "sf","gfcanalysis",
#               "doParallel", "foreach", "parallel")


