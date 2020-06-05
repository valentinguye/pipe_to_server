### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# This script adds variables to the parcel panel data frame. 
# In particular, it compute the number of UML mills each parcel can reach within 10, 30 and 50km. 
# And it adds geographic variables (island and districts). 
#
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
  
  CR <- 10 # in km 
  while(CR < 60){
    parcels$newv_uml <- rep(0, nrow(parcels))

    for(t in 1:length(years)){
      
      # UML
      # This is not a panel, so the information on presence or not a given year is whether 
      # the establishment year is anterior. We impute NA establishment year to be older than 1998. 
      present_uml <- uml[uml$est_year_imp <= years[t] | is.na(uml$est_year_imp),]

      annual_reachable_uml <- st_is_within_distance(parcels_centro, present_uml, dist = CR*1000)
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
  parcels[parcels$newv_uml != 0, "newv_ibs_uml"] <- 100*(parcels[parcels$newv_uml != 0, "n_reachable_ibs"]/parcels[parcels$newv_uml != 0,"newv_uml"])
    
  colnames(parcels)[colnames(parcels) == "newv_uml"] <- paste0("n_reachable_uml_",CR,"km")
  colnames(parcels)[colnames(parcels) == "newv_ibs_uml"] <- paste0("sample_coverage_",CR,"km")
  
  CR <- CR + 20
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




##### ADD GEGRAPHIC VARIABLES ##### 

parcel_size <- 3000
catchment_radiuseS <- c(1e4, 3e4, 5e4)
for(catchment_radius in catchment_radiuseS){
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_reachable_uml_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  parcels <- st_as_sf(parcels, coords = c("lon", "lat"), crs = indonesian_crs)
  
  # ISLAND variable
  island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
  names(island_sf)[names(island_sf)=="island"] <- "shape_des"
  
  island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)
  
  parcels$island <- rep("", nrow(parcels))
  
  # make the operation faster by using island bbox (other wise the island polygons make 
  # the computation very long)
  # (and this also includes parcel centroids in the sea)
  island_sf_prj_bbox <- sapply(island_sf_prj$geometry, function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
  
  sgbp <- st_within(parcels$geometry, island_sf_prj_bbox)
  # the bboxes of Sumatra and Kalimantan intersect a bit, so we check that no parcel falls 
  # in the intersection, this is the case for catchment_radius = 50km, for 2538 parcels, 
  # these parcel centroids belong to Kalimantan (after visual check). 
  # intersect <- st_intersection(island_sf_prj_bbox[1], island_sf_prj_bbox[3])
  # plot(island_sf_prj_bbox[[1]])
  # plot(island_sf_prj, add = TRUE)
  # plot(parcels$geometry[parcels$island==4], col = "red", add = TRUE)
  
  sgbp[lengths(sgbp)==2] <- 3
  
  # island_sf_prj features are in this order : 1 Sumatra; 2 Papua; 3 Kalimantan
  unique(unlist(sgbp))
  parcels$island <- unlist(sgbp)
  
  parcels$island <- replace(parcels$island, parcels$island == 1, "Sumatra")
  unique(parcels$island)
  parcels$island <- replace(parcels$island, parcels$island == 2, "Papua")
  parcels$island <- replace(parcels$island, parcels$island == 3, "Kalimantan")
  
  
  
  # DISTRICT variable
  district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
  
  district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)
  
  dlm <- st_contains(district_sf_prj, parcels$geometry, sparse = FALSE)
  
  nrow(dlm)
  row.names(dlm) <- district_sf_prj$d__2000
  
  districts <- lapply(c(1:ncol(dlm)), FUN = function(col){row.names(dlm)[dlm[,col]==T]})
  
  districts[lengths(districts)==0] %>% length()
  
  # 2952, 73044 and (for 10, 30 and 50km CR resp.) grid cells have their centroids in the sea. 
  # We did not want to discard them because they may have interesting patterns in the cell's 
  # part that is on ground. 
  districts[lengths(districts)==0] <- "sea"
  
  parcels$district <- unlist(districts)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_more_variables_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"CR.rds")))
}












