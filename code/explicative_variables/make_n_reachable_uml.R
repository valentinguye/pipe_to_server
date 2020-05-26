### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Program to compute the number of UML mills each parcel can reach within 10, 30 and 50km. 
# 
#   Inputs: UML most complete version
#           --> UML_valentin_imputed_est_year.dta
# 
#           parcel panel from previous step (i.e. making weighted averages)
#           --> pattern: wa_panel_parcels_ ; for each parcel_size and catchment_radius combination
# 
#   Outputs: parcel panel with 3 new columns: the parcel and time varying numbers of UML mills reachable within 10, 30 and 50km. 
#           --> pattern panel_parcels_reachable_uml_ ; for each parcel_size and catchment_radius combination
# 
# 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "readstata13", 
                   "rgdal", "sf")
#install.packages("sf", source = TRUE)
# library(sf)
# 
# neededPackages = c("tidyverse","data.table", "readxl","foreign", "data.table", "readstata13", "here",
#                    "rgdal", "raster", "velox","sp", "lwgeom", "rnaturalearth", 
#                    "rlist", "parallel", "foreach", "iterators", "doParallel" )
# 

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Load them 
lapply(neededPackages, library, character.only = TRUE)

# /!\/!\ IF renv::restore(neededPackages) FAILS TO INSTALL SOME PACKAGES /!\/!\ 

# For instance sf could cause trouble https://github.com/r-spatial/sf/issues/921 
# or magick, as a dependency of raster and rgdal. 

# FOLLOW THESE STEPS:
# 1. Remove these package names from neededPackages above, and rerun renv::restore(packages = neededPackages)
# 2. Write them in troublePackages below, uncomment, and run the following code chunk: 

# # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("")
# # Attempt to load packages from user's default libraries. 
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 



### NEW FOLDERS USED IN THIS SCRIPT 


# #install.packages("sf", source = TRUE)
# library(sf)
# 
# neededPackages = c("tidyverse","data.table", "readxl","foreign", "data.table", "readstata13", "here",
#                    "rgdal", "raster", "velox","sp", "lwgeom", "rnaturalearth", 
#                    "rlist", "parallel", "foreach", "iterators", "doParallel" )


### INDONESIAN CRS 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### IBS YEARS 
years <- seq(1998, 2015, 1)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



# read the most complete version of UML we have. 
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
uml <- uml[!is.na(uml$lat),]
uml <- st_as_sf(uml, coords = c("lon", "lat"), crs = 4326)
uml <- st_transform(uml, crs = indonesian_crs)

# read the sample panel of IBS geolocalized mills
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
ibs <- ibs[!is.na(ibs$lat),]
length(unique(ibs$firm_id))
class(ibs$year)
ibs_cs <- lapply(years, FUN = function(x) ibs[ibs$year == x,]) 
ibs_cs <- lapply(ibs_cs, FUN = st_as_sf, coords =  c("lon", "lat"), remove = TRUE, crs = 4326)
ibs_cs <- lapply(ibs_cs, FUN = st_transform, crs = indonesian_crs)
ibs_cs <- lapply(ibs_cs, FUN = st_geometry)

make_n_reachable_uml <- function(parcel_size, catchment_radius){
  
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  # make a spatial cross section of it (parcels' coordinates are constant over time)
  parcels_centro <- parcels[parcels$year == 1998, c("parcel_id", "lat", "lon")]
  # (lon lat are already expressed in indonesian crs)
  parcels_centro <- st_as_sf(parcels_centro, coords = c("lon", "lat"), remove = T, crs = indonesian_crs)
  
  CR <- 10000
  while(CR < 60000){
    parcels$newv_uml <- rep(0, nrow(parcels))

    for(t in 1:length(years)){
      
      # UML
      # This is not a panel, so the information on presence or not a given year is whether 
      # the establishment year is anterior. We impute NA establishment year to be older than 1998. 
      present_uml <- uml[uml$est_year_imp <= years[t] | is.na(uml$est_year_imp),]

      annual_reachable_uml <- st_is_within_distance(parcels_centro, present_uml, dist = CR)
      parcels[parcels$year == years[t], "newv_uml"] <- lengths(annual_reachable_uml)

      # IBS
      # n_reachable_ibs was already computed in wa_at_parcels.R
      # We used the year variable from the IBS panel to determine whether a firm was present in a given year. 
      # This means that when a firm has a yearly record missing, although we know it was there 
      # that year bc we observe an older and an earlier records, we do not count it as reachable 
      # by the parcel, because no IBS information would be usable that year in such a case. 
    
    }
    
    
  # IBS/UML (sample coverage)
  # any reachable ibs is also counted as a reachable uml
  nrow(parcels[parcels$newv_uml< parcels$n_reachable_ibs,])
  
  # we give to each parcel a ratio that informs on how much of the total influence it gets from all (uml) 
  # reachable mills is observed in our sample (ibs)  
  parcels$newv_ibs_uml <- rep(0, nrow(parcels))  
  parcels[parcels$newv_uml != 0, "newv_ibs_uml"] <- parcels[parcels$newv_uml != 0, "n_reachable_ibs"]/parcels[parcels$newv_uml != 0,"newv_uml"]
    
  colnames(parcels)[colnames(parcels) == "newv_uml"] <- paste0("n_reachable_uml_",CR,"km")
  colnames(parcels)[colnames(parcels) == "newv_ibs_uml"] <- paste0("sample_coverage_",CR,"km")
  
  CR <- CR + 20000
  }

  return(parcels)
}

parcel_size <- 3000
catchment_radius <- 10000
while(catchment_radius < 60000){
  
  make_n_reachable_uml(parcel_size, catchment_radius) %>% 
    saveRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_reachable_uml_",
                             parcel_size/1000,"km_",
                             catchment_radius/1000,"CR.rds")))
  
  catchment_radius <- catchment_radius + 20000
}





# parcels_centro_tmp <- readRDS(paste0("./outcome_variables/dataframes/panel_Indonesia_",parcel_size/1000,"km_",catchment_radius/1000,"CR.rds"))
# # keep only one cross-section
# parcels_centro_tmp <- filter(parcels_centro_tmp, year == 2001)
# parcels_centro_tmp <- dplyr::select(parcels_centro_tmp, parcel_id, lat, lon)
# parcels_centro_tmp <- dplyr::arrange(parcels_centro_tmp, parcel_id)
# 
# 
# make_n_reachable_uml <- function(parcel_size, catchment_radius){
#   
#   # read the parcel panel
#   parcels <- readRDS(here(paste0("/build/output/wa_panel_parcels_",
#                                  parcel_size/1000,"km_",
#                                  catchment_radius/1000,"CR.rds")))
#   
#   parcels_centro <- left_join(x = parcels, y = parcels_centro_tmp,  by = "parcel_id")
# 









