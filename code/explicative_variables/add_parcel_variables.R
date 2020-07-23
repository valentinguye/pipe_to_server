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
#           island polygons
#           --> temp_data/processed_indonesia_spatial/island_sf
#   
#           Province polygons
#           --> input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
#
#           district polygons and names (prepare in prepare_crosswalks.do)
#           --> input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp
#           --> temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
# 
#           baseline forest extents variables in cross-section (prepare in prepare_2000_forest_extents.R)
#           --> pattern: temp_data/processed_parcels/baseline_fc_cs_", for each parcel_size and IBS catchment_radius combination 
#
#
#   Outputs: parcel panel with new columns: the parcel and time varying numbers of UML mills reachable within 10, 30 and 50km. 
#                                           the island, the district, and the pixelcounts and areas (ha) of 2000 forest extents
#                                           (for 30% tree canopy density outside indsutrial plantations and total primary forest. 
#           --> pattern temp_data/processed_parcels/parcels_panel_reachable_uml_
#                       temp_data/processed_parcels/parcels_panel_geovars_
#                       temp_data/processed_parcels/parcels_panel_final_ 
#                       ; for each parcel_size and catchment_radius combination
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
neededPackages = c("dplyr", "readstata13", "foreign",
                   "rgdal", "sf", 
                   "DataCombine")
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

# 3. If the troubling packages could not be loaded ("there is no package called ...) 
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

### GRID CELL SIZE
parcel_size <- 3000

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# This takes ~1h 
#### ADD N REACHABLE UML AND SAMPLE COVERAGE ####

# The sample coverage has to be computed as the number of reachable IBS-UML matched sample
# related to the number of reachable UML mills, with the all the former being included in the latter. 
# We cannot compute a ratio of all mills in the analysis sample - i.e. also including IBS not matched with 
# UML but with a desa centroid - because we would not know whether a mill from the latter group 
# is an additional mill that is not in UML or if it is in UML but it was not matched. 

# read the sample panel of IBS geolocalized mills
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# keep only those that are matched with UML. 
ibsuml <- ibs[ibs$uml_matched_sample==1,]
# make it a cross section
ibsuml <- ibsuml[!duplicated(ibsuml$firm_id),]
ibsuml <- ibsuml[!is.na(ibsuml$lat),]
ibsuml <- st_as_sf(ibsuml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
ibsuml <- st_transform(ibsuml, crs = indonesian_crs)

# read the most complete version of UML we have. 
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
uml <- uml[!is.na(uml$lat),]
uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
uml <- st_transform(uml, crs = indonesian_crs)


# so we do not do that: 
# # all mills are IBS-UML matched + IBS not matched with UML but with desa centroid)
# all_mills <- bind_rows(uml[,c("trase_code", "year","lon","lat")], 
#                        ibs[ibs$analysis_sample==1 & ibs$uml_matched_sample == 0,c("firm_id", "year", "lon", "lat")])
# 
# ibs[ibs$analysis_sample==1 & ibs$uml_matched_sample == 0,"firm_id"] %>% unique() %>% length()
# length(unique(uml$trase_code))
# nrow(unique(dplyr::select(all_mills, -year)))
# break it down to cross sections just like for IBS above. 
# class(all_mills$year)
# all_mills_cs <- lapply(years, FUN = function(x) all_mills[all_mills$year == x,]) 
# all_mills_cs <- lapply(all_mills_cs, FUN = st_as_sf, coords =  )
# all_mills_cs <- lapply(all_mills_cs, FUN = st_transform, crs = indonesian_crs)
# all_mills_cs <- lapply(all_mills_cs, FUN = st_geometry)
### ### ###

make_n_reachable_uml <- function(parcel_size, catchment_radius){
  
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  # make a spatial cross section of it (parcels' coordinates are constant over time)
  parcels_centro <- parcels[parcels$year == 1998, c("parcel_id", "lat", "lon")]
  # (lon lat are already expressed in indonesian crs)
  parcels_centro <- st_as_sf(parcels_centro, coords = c("lon", "lat"), remove = T, crs = indonesian_crs)
  
  parcels$newv_uml <- rep(0, nrow(parcels))
  parcels$newv_ibsuml <- rep(0, nrow(parcels))
    
  for(t in 1:length(years)){
      
      # UML
      # This is not a panel, so the information on presence or not a given year is whether 
      # the establishment year is anterior. We impute NA establishment year to be older than 1998. 
      present_uml <- uml[uml$est_year_imp <= years[t] | is.na(uml$est_year_imp),]

      annual_reachable_uml <- st_is_within_distance(parcels_centro, present_uml, dist = catchment_radius)
      parcels[parcels$year == years[t], "newv_uml"] <- lengths(annual_reachable_uml)

      # IBS-UML
      present_ibsuml <- ibsuml[ibsuml$est_year_imp <= years[t] | is.na(ibsuml$est_year_imp),]
      
      annual_reachable_ibsuml <- st_is_within_distance(parcels_centro, present_ibsuml, dist = catchment_radius)
      parcels[parcels$year == years[t], "newv_ibsuml"] <- lengths(annual_reachable_ibsuml)
      
      # Note that n_reachable_ibs was already computed in wa_at_parcels.R, 
      # but this includes ibs that are not matched with uml, which we do not want to count here. 
      # We used the year variable from the IBS panel to determine whether a firm was present in a given year. 
      # This means that when a firm has a yearly record missing, although we know it was there 
      # that year bc we observe an older and an earlier records, we do not count it as reachable 
      # by the parcel, because no IBS information would be usable that year in such a case.
      # But when the mill has a record line in IBS this year, but some or all information is missing we still 
      # count the mill as one more being reachable, altough we use no info from it.
    
    }
    
    
  # IBS/all mills (IBS-UML matched + IBS not matched with UML but with desa centroid) -> sample coverage
  # any reachable ibs is also counted as a reachable uml : YES if ibs is defined as ibs[ibs$uml_matched_sample==1,]
  nrow(parcels[parcels$newv_uml< parcels$newv_ibsuml,])
  
  # we give to each parcel a ratio that informs on the share of the total influence (from all possible mills known)
  # that is catched by our sample of analysis, i.e. geo-localized palm oil mills. 
  parcels$ratio <- rep(0, nrow(parcels))  
  parcels[parcels$newv_uml != 0, "ratio"] <- 100*(parcels[parcels$newv_uml != 0, "newv_ibsuml"]/parcels[parcels$newv_uml != 0,"newv_uml"])
    
  colnames(parcels)[colnames(parcels) == "newv_ibsuml"] <- paste0("n_reachable_ibsuml")
  colnames(parcels)[colnames(parcels) == "newv_uml"] <- paste0("n_reachable_uml")
  colnames(parcels)[colnames(parcels) == "ratio"] <- paste0("sample_coverage")
  
  return(parcels)
}

catchment_radius <- 10000
while(catchment_radius < 60000){
  
  make_n_reachable_uml(parcel_size, catchment_radius) %>% 
    saveRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_reachable_uml_",
                             parcel_size/1000,"km_",
                             catchment_radius/1000,"CR.rds")))
  
  catchment_radius <- catchment_radius + 20000
}





#### ADD GEGRAPHIC VARIABLES AND BASELINE FOREST EXTENT VARIABLES ####
catchment_radiuseS <- c(1e4, 3e4, 5e4)#
for(catchment_radius in catchment_radiuseS){
  # this is prepared in this script's previous part (add n reachable uml... )
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_reachable_uml_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  parcels <- st_as_sf(parcels, coords = c("lon", "lat"), crs = indonesian_crs, remove = FALSE)
  
  
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
  
  
  # PROVINCE variable
  
  # Work with a cross section for province and district attribution
  parcels_cs <- parcels[!duplicated(parcels$parcel_id),]

  
  province_sf <- st_read(file.path("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp"))
  province_sf <- dplyr::select(province_sf, NAME_1)
  province_sf_prj <- st_transform(province_sf, crs = indonesian_crs)
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_prov_idx <- st_nearest_feature(parcels_cs, province_sf_prj)
  
  parcels_cs$province <- province_sf_prj$NAME_1[nearest_prov_idx]
  
  # DISTRICT variable
  district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
  district_names <- read.dta13(file.path("temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta"))
  district_names <- district_names[!duplicated(district_names$bps_),]
  district_sf$d__2000 <- district_sf$d__2000 %>% as.character()
  district_names$bps_ <- district_names$bps_ %>% as.character()
  district_sf <- left_join(x = district_sf, y = district_names[,c("name_", "bps_")], 
                       by = c("d__2000" = "bps_"),
                       all = FALSE, all.x = FALSE, all.y = FALSE)
  
  
  district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_dstr_idx <- st_nearest_feature(parcels_cs, district_sf_prj)
  
  # 4 parcels are closest to district with no name (NA) 
  parcels_cs$district <- district_sf_prj$name_[nearest_dstr_idx]
  
  parcels <- merge(st_drop_geometry(parcels),
                   st_drop_geometry(parcels_cs[,c("parcel_id", "province", "district")]),
                   by = "parcel_id")
 
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"CR.rds")))
}

  # parcels_list[[match(catchment_radius, catchment_radiuseS)]] <- parcels




#### BASELINE FOREST EXTENT VARIABLES AND TIME DYNAMICS VARIABLES ####
catchment_radiuseS <- c(1e4, 3e4, 5e4)#
for(catchment_radius in catchment_radiuseS){ 
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                                  parcel_size/1000,"km_",
                                                  catchment_radius/1000,"CR.rds")))
  
  # this is a cross section, computed in prepare_2000_forest_extents.R
  bfe <- readRDS(file.path(paste0("temp_data/processed_parcels/baseline_fc_cs_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"km_IBS_CR.rds")))
  bfe <- dplyr::select(bfe, parcel_id, lat, lon, everything())
  
  # merge it with our parcel panel
  parcels <- base::merge(parcels, 
                         dplyr::select(bfe, -lat, -lon), 
                         by = "parcel_id")
  
  # test that the parcel_id have been attributed to parcels equally in the two processes 
  # (prepare_lucpfip.R and prepare_2000_forest_extents.R)
  # 
  # parcels2 <- base::merge(parcels, bfe, by = c("parcel_id", "lat","lon"))
  # setorder(parcels1, parcel_id, year)
  # setorder(parcels2, parcel_id, year)
  # row.names(parcels1) <- NULL
  # row.names(parcels2) <- NULL
  # all.equal(st_drop_geometry(parcels1[,c("parcel_id", "lon","lat")]), 
  #           st_drop_geometry(parcels2[,c("parcel_id", "lon","lat")]))
  # 
  # all(names(parcels1)==names(parcels2))
  # all.equal(st_drop_geometry(parcels1), st_drop_geometry(parcels2))
  # rm(parcels2)
  
  
  
#### TIME DYNAMICS VARIABLES #### 
  #parcels <- st_drop_geometry(parcels)
    
  ### Simple lags and leads on a large set of variables
  variables <- c("wa_ffb_price_imp1", "wa_ffb_price_imp2", 
                "wa_cpo_price_imp1", "wa_cpo_price_imp2", "wa_prex_cpo_imp1","wa_prex_cpo_imp2",       
                "wa_pko_price_imp1",       "wa_pko_price_imp2",       "wa_prex_pko_imp1",        "wa_prex_pko_imp2",       
                "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp",  "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",     
                #"wa_concentration_10",     "wa_concentration_30", "wa_concentration_50",     
                "n_reachable_ibs", "n_reachable_uml", "n_reachable_ibsuml",   "sample_coverage")
  
  for(voi in variables){
    ## lags
    for(lag in c(1:5)){
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = voi, 
                                    TimeVar = "year",
                                    GroupVar = "parcel_id",
                                    NewVar = paste0(voi,"_lag",lag),
                                    slideBy = -lag, 
                                    keepInvalid = TRUE)
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      
    }
      
    # ## leads                               
    # for(lag in c(1:5)){
    #   parcels <- dplyr::arrange(parcels, parcel_id, year)
    #   parcels <- DataCombine::slide(parcels,
    #                                 Var = voi, 
    #                                 TimeVar = "year",
    #                                 GroupVar = "parcel_id",
    #                                 NewVar = paste0(voi,"_lead",lag),
    #                                 slideBy = lag, 
    #                                 keepInvalid = TRUE) 
    #   parcels <- dplyr::arrange(parcels, parcel_id, year)
    # } 
  }
  
  parcels1 <- parcels

  ### Operations relating contemporaneous to past information - on prices only
  variables <- c("wa_ffb_price_imp1", "wa_ffb_price_imp2", 
                 "wa_cpo_price_imp1", "wa_cpo_price_imp2",        
                 "wa_pko_price_imp1", "wa_pko_price_imp2")
  
  for(voi in variables){
    
    ## Past-year averages (2, 3 and 4 years) - LONG RUN MEASURE - 
    
    for(py in c(2,3,4)){
      parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_lag",c(1:py))], na.rm = TRUE)
      parcels[is.nan(parcels$newv),"newv"] <- NA
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_",py,"pya")
      
      ## and absolute deviation - SHORT RUN MEASURE -
      parcels <- mutate(parcels,
                        !!as.symbol(paste0(voi,"_dev_",py,"pya")) := !!as.symbol(paste0(voi)) - 
                                                                     !!as.symbol(paste0(voi,"_",py,"pya")))
      # # and relative deviation
      # parcels <- mutate(parcels,
      #                   !!as.symbol(paste0(voi,"_rdev_",py,"pya")) := (!!as.symbol(paste0(voi)) - 
      #                                                               !!as.symbol(paste0(voi,"_",py,"pya"))) /
      #  
      
      # Lag these deviations by one year
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                  Var = paste0(voi,"_dev_",py,"pya"), 
                                  TimeVar = "year",
                                  GroupVar = "parcel_id",
                                  NewVar = paste0(voi,"_dev_",py,"pya_lag1"),
                                  slideBy = -1, 
                                  keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      
      # lag also just the past year average (useful if we use those per se. as measures of LR price signal)
      # note that 3pya_lag1 is different from 4pya. 
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_",py,"pya"), 
                                    TimeVar = "year",
                                    GroupVar = "parcel_id",
                                    NewVar = paste0(voi,"_",py,"pya_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, parcel_id, year)
    }
    
    
    
    
    ### Contemporaneous and past-year averaged year-on-year growth 
    
    ## contemporaneous yoyg - SHORT RUN MEASURE - (invalid for at least the first record of each parcel_id)
    parcels <- mutate(parcels,
                      !!as.symbol(paste0(voi,"_yoyg")) := 100*(!!as.symbol(paste0(voi)) - 
                                                               !!as.symbol(paste0(voi,"_lag1"))) /
                                                               !!as.symbol(paste0(voi,"_lag1")))
                        
    ## Lagged yoyg (this is only a step)
    # (the first lag is invalid for at least two first records of each parcel_id;    
    # the fourth lag is invalid for at least 5 first records)
    for(lag in c(1:4)){
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_yoyg"), 
                                    TimeVar = "year",
                                    GroupVar = "parcel_id",
                                    NewVar = paste0(voi,"_yoyg_lag",lag),
                                    slideBy = -lag, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, parcel_id, year)
    }
    
    
    ## past-years averaged yoyg - LONG RUN MEASURE - 
    for(py in c(2,3,4)){
      parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_yoyg_lag",c(1:py))], na.rm = TRUE)
      # treat NaNs that arise from means over only NAs when na.rm = T 
      parcels[is.nan(parcels$newv),"newv"] <- NA
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_yoyg_",py,"pya")
    
    
    # Lag by one year
    parcels <- dplyr::arrange(parcels, parcel_id, year)
    parcels <- DataCombine::slide(parcels,
                                  Var = paste0(voi,"_yoyg_",py,"pya"), 
                                  TimeVar = "year",
                                  GroupVar = "parcel_id",
                                  NewVar = paste0(voi,"_yoyg_",py,"pya_lag1"),
                                  slideBy = -1, 
                                  keepInvalid = TRUE)  
    parcels <- dplyr::arrange(parcels, parcel_id, year)
    }
  }  

  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))  
}

# 
# voi <- "wa_ffb_price_imp1"
# View(parcels[parcels$parcel_id==16500,c("parcel_id", "year",
#                                         voi,#
#                                         paste0(voi,"_lag",c(1:3)),#
#                                         paste0(voi,"_3pya"),#
#                                         paste0(voi,"_3pya_lag1"),#
#                                         paste0(voi,"_dev_3pya"),#
#                                         paste0(voi,"_dev_3pya_lag1"),#
#                                         paste0(voi,"_yoyg"),#
#                                         paste0(voi,"_yoyg_lag", c(1:3)),#
#                                         paste0(voi,"_yoyg_3pya"),#
#                                         paste0(voi,"_yoyg_3pya_lag1"))])#


  





















#   dlm <- st_contains(district_sf_prj, parcels$geometry, sparse = FALSE)
#   
#   nrow(dlm)
#   class(dlm[,1])
#   row.names(dlm) <- district_sf_prj$name_
#   
#   # to each column of dlm, i.e. to each grid cell, attribute the name of the district 
#   # it is contained by. 
#   districts <- lapply(c(1:ncol(dlm)), FUN = function(col){row.names(dlm)[dlm[,col]==T]})
#   
#   districts[lengths(districts)==0] %>% length()
#   
#   # 6264, 101250 and (for 10 and 30 CR resp.) grid cells have their centroids in the sea. 
#   # We did not want to discard them because they may have informative patterns in the cell's 
#   # part that is on ground. 
#   
#   
#   
#   districts[lengths(districts)==0] <- "sea"
#   
#   parcels$district <- unlist(districts)
# 
#   nearest_dstr_idx <- st_nearest_feature(parcels[parcels$district=="sea", "geometry"], district_sf_prj)
#   # with this operation, some parcels get NA for an index of nearest district (from district_sf_prj)
#   # I don't know where it comes from. It is 72 records (4 parcels). We give them no name (NA). 
#   parcels1[parcels1$district=="sea" & is.na(nearest_dstr_idx1), "district"] <- ""
# 
#   parcels1$district[parcels1$district=="sea"][!is.na(nearest_dstr_idx1)]  <- district_sf_prj$name_[nearest_dstr_idx1][!is.na(nearest_dstr_idx1)]  
# all.equal(district_sf_prj$name_[nearest_dstr_idx1][!is.na(nearest_dstr_idx1)] , district_sf_prj$name_[!is.na(nearest_dstr_idx1)] )
#   
#   nearest_dstr_idx1[is.na(district_sf_prj$name_[nearest_dstr_idx1])]
#   
#   parcels1[parcels1$seadis == "Kab. Aceh Barat", "geometry"] %>% plot()
#   parcels1[parcels1$seadis == "sea" & parcels1$district == "Kab. Aceh Barat", "geometry"] %>% plot(col = "red", add = T)
#   district_sf_prj[district_sf_prj$name_=="Kab. Aceh Barat",]$geometry %>% plot(add =T)