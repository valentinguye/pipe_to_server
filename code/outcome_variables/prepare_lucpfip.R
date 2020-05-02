### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# 
#       BUILDING PANEL DATAFRAMES OF LAND USE CHANGE FROM PRIMARY FOREST TO INDUSTRIAL PLANTATION (LUCPFIP)
# 
#   Inputs:   * Island polygons 
#             ---> temp_data/processed_indonesia_spatial/island_sf 
#
#             * Global Forest Change (Hansen et al. 2013) 
#             ---> tiles downloaded from internet from this script ///!!!\\\ THIS SCRIPT REQUIRES INTERNET CONNEXION ///!!!\\\
#         
#             * 2015 oil palm plantations (Austin et al. 2017) for Sumatra, Kalimantan and Papua,
#             ---> oilpalm_2015_WGS1984.tif 
# 
#             * Primary forest extent in 2000 (Margono et al. 2014) 
#             ---> timeseq_change00_12.tif
#      
#             * Georeferenced mills (from georeferencing works)                                           
#             ---> IBS_UML_cs.dta
#
#   Outputs:  panel dataframes of LUCPFIP pixel-event count in parcels of a given size, from 2001 to 2018,  
#             for the whole Indonesia (Sumatra, Kalimantan, Papua "row-binded"), 
#             for 3 forest definitions (30% canopy closure in intact, degraded, and intact or degraded (total) primary forest).
#             
#             There is one such dataframe for each combination of parcel size (only 3x3km for now) and catchment radius (10, 30, 50km)
#             ---> lucpfip_panel_3km_10CR.rds 
#             ---> lucpfip_panel_3km_30CR.rds 
#             ---> lucpfip_panel_3km_50CR.rds
# 
#   
#   Actions:  This script consists of mainly three functions.  
#             0. load needed packages; set working directory; set raster options; define the crs used throughout the script. 
#                 /// !!! \\\ the chunksize and maxmemory raster options should be set accordingly with the machine used,  
#                             considering that this script executes parallel functions using parallel::detectCores() - 1
#                             For instance, here we set chunksize to 1Go so that our 3 working cores processed 3Go together. 
#
#             1. prepare_pixel_lucpfip(island)
# 
#             2. aggregate_lucpfip(island, parcel_size)
# 
#             3. to_panel_within_CR(island, parcel_size, catchment_radius)
#
#             Finally, functions are run and there outputs are merged across islands and primary forest types.   
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###

# sf need to be installed from source for lwgeom te be installed from source. 
if (!require(sf)) install.packages("sf", source = TRUE)
#"plyr", 
neededPackages = c("dplyr", "raster", "foreign", "sp", "lwgeom", "rnaturalearth", "data.table",
                   "rgdal", "readstata13",
                   "rlist", "velox", "parallel", "foreach", "iterators", "doParallel", "readxl", "here")
allPackages    = c(neededPackages %in% installed.packages()[ , "Package"]) 

# Install packages (if not already installed) 
if(!all(allPackages)) {
  missingIDX = which(allPackages == FALSE)
  needed     = neededPackages[missingIDX]
  lapply(needed, install.packages)
}

# Load all defined packages
lapply(neededPackages, library, character.only = TRUE)
library(sf)

# install other packages not from source.
if (!require(devtools)) install.packages("devtools")
library(devtools)

# package tictoc
install_github("jabiru/tictoc")
library(tictoc)

#install.packages("sf", source = TRUE)
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("r-spatial/lwgeom")
#library(lwgeom)

#INSTALL GFC ANALYSIS
# Install the snow package used to speed spatial analyses
if (!require(snow)) install.packages("snow")
library(snow)

# Install Alex's gfcanalysis package
if (!require(gfcanalysis)) install.packages('gfcanalysis')
library(gfcanalysis)


### NEW FOLDERS USED IN THIS SCRIPT 

# this is where GFC tiles are stored 
dir.create("temp_data/GFC_tiles")

# this is where will be stored temporary raster function outputs. 
dir.create("temp_data/processed_lu")

# annual layers will be stored in a subfolder just for the sake of tidyness in processed_lu
dir.create("temp_data/processed_lu/annual_maps")

# dataframes of parcels are stored here 
dir.create("temp_data/processed_parcels")


### RASTER OPTIONS ### 
# Do not change chunksize, as it is not safe in combinatin with clusterR and machine specific. 
rasterOptions(timer = TRUE, 
              tmpdir = "temp_data/raster_tmp")

### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Applied to GFC data, this makes a resolution of 27.6 ; 27.8 meters.
# Applied to GFC data, the crs from Margono et al. 2014 
# (+proj=sinu +lon_0=140 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs )
# makes a resolution of 27.7 ; 27.8 meters. 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# preliminary: read the island shapefile
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "island_name"

##### 1. PREPARE 30m PIXEL-LEVEL MAPS OF LUCFIP ##### 

prepare_pixel_lucpfip <- function(island){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Define area of interest (AOI) ####
  
    aoi <- island_sf[island_sf$island_name == island,]
    aoi <- st_as_sfc(st_bbox(aoi))
    aoi_sp <- as(aoi, "Spatial")
    rm(aoi)
    
    # it is not projected # 
    
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Prepare loss layer from GFC data  ####
  
    # Download, extract, threshold, project forest loss layer. 
    
    ### Download GFC data
    
    #define where all tiles are going to be stored
    data_folder <- file.path("temp_data/GFC_tiles")

    # Calculate tiles needed to cover the AOI
    tiles <- calc_gfc_tiles(aoi_sp)
    # length(tiles) 

    # version of GFC used here.
    gfc_version <- "GFC-2018-v1.6"
    # //!\\ this script is not written flexibly to adjust for other versions of GFC. One should check every "18" entries in this script for instance. 

    #download tiles - with all layers otherwise later extract_gfc does not work
    download_tiles(tiles, 
                   output_folder = data_folder, 
                   images = c("treecover2000", "lossyear", "gain", "datamask"), 
                   dataset = gfc_version)

    rm(tiles)
    
    
    ### extract GFC data 
    
    # (can only extract all layers with default stack=change)
    # to better understand extract_gfc see https://rdrr.io/cran/gfcanalysis/src/R/extract_gfc.R
    extract_gfc(aoi_sp, data_folder,
                stack = "change",
                to_UTM = FALSE,
                dataset = gfc_version,
                filename = file.path(paste0("temp_data/processed_lu/gfc_data_",island,".tif")),
                overwrite = TRUE )


    
    ### Threshold GFC data and select the loss layer
    
    # The forest loss layer is based on a forest definition. This requires to specify a 
    # pixel-level canopy cover percentage as the threshold between forest and non-forest state in 2000.
    # Here, this is 30%. 
    # function threshold_gfc is already defined in package gfc_analysis
    
    gfc_data <- brick(file.path(paste0("temp_data/processed_lu/gfc_data_",island,".tif")))
    
    threshold_gfc(gfc_data,
                  forest_threshold=30,
                  filename= file.path(paste0("temp_data/processed_lu/gfc_data_",island,"_30th.tif")),
                  overwrite = TRUE)
    
    rm(gfc_data)
    
    ## select the loss layer (15 ad 40 are arbitrary, the max value of loss layer is an integer corresponding 
    # to the latest year after 2000 when loss is observed. We need a GFC version with this year being at least 2015. 
    # and we do not want to select a layer with percentage and values up to 100. 
    
    thed_gfc_data <- brick(file.path(paste0("temp_data/processed_lu/gfc_data_",island,"_30th.tif"))) 

    loss <- thed_gfc_data[[which(thed_gfc_data@data@max > 15 & thed_gfc_data@data@max < 40)]]
   
    rm(thed_gfc_data)
    
    ### Project loss layer
    
    # This is necessary because we will need to make computations on this map within mills' catchment *areas*.
    # If one does not project maps, then catchment areas all have different areas while being defined with a common buffer.
    
    beginCluster() # this uses by default detectCores() - 1 
    
    projectRaster(loss,
                  method = "ngb",
                  crs = indonesian_crs, 
                  filename = file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")), 
                  dataType = "INT1U",
                  overwrite = TRUE)
    
    endCluster()
  
    rm(loss) 
    removeTmpFiles(h=0)
    
    print(paste0("complete ", "temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif"))
      
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Prepare plantation maps ####

  ### Austin et al. (2017) 2015 industrial oil palm plantation map - ioppm2015 
  ioppm2015 <- raster(file.path("input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_2015_WGS1984.tif")) 
  
  # As such, these raw data are lat-lon (WGS1984)
  # extent is Indonesia  
  # resolution is 0.002277, 0.002277
  # while gfc_data resolution is 0.00030, 0.00025
  # Cells have value 1 for oil palm plantation, NA else. 
    
  # first adjust ioppm2015 maps to the same extent as gfc_data, i.e. the island bounding box. 
  # intermediary step to reclassify NA to 0 (not time consuming)
  # then projectRaster to match crs res - disaggregate before is not necessary (yields the same result)
    

  ### Adjust to roughly the same extent as GFC loss, i.e. island bbox.  
  # First, crop from Indonesia wide to island extent (in the longitude mainly)
  # then, extend it (in the latitude mainly) because raw data do not cover northern part of Sumatra. 
  
  cropped_ioppm2015 <- raster::crop(ioppm2015, y = aoi_sp)
  
  raster::extend(cropped_ioppm2015, y = aoi_sp, 
                 value = NA, # the default
                 filename = file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,".tif")),
                 datatype = "INT1U",
                 overwrite = TRUE)
  
  rm(cropped_ioppm2015, ioppm2015)
  
  
  ### Reclassify NA into 0 
  ioppm2015 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,".tif")))
  
  raster::reclassify(ioppm2015, 
                     rcl = cbind(NA,0), 
                     filename = file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,"_reclassified.tif")), 
                     overwrite=TRUE, 
                     datatype = "INT1U")

   
  ### Align to GFC loss crs, resolution and exact extent.  

  # read the ioppm2015 map (with roughly the island extent and reclassified) 
  ioppm2015 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,"_reclassified.tif")))
  
  # Read the target GFC loss layer
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))

  beginCluster() # this uses by default detectCores() - 1
  
  projectRaster(from = ioppm2015, to = loss,
                method = "ngb",
                filename = file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,"_aligned.tif")),
                datatype = "INT1U",
                overwrite = TRUE )
  endCluster()
  # ~4200s. 

  rm(ioppm2015, loss)
  removeTmpFiles(h=0)  
  
  print(paste0("complete ", "temp_data/processed_lu/austin_ioppm_2015_",island,"_aligned.tif"))
  
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Prepare primary forest map ####    
  
  ### Primary forest extent map is downloaded from https://glad.umd.edu/dataset/primary-forest-cover-loss-indonesia-2000-2012
  pf <- raster(file.path("input_data/margono_primary_forest/timeseq_change00_12.tif"))
  
  # Note: there is something strange with this raster, its min value is 1, while it should be 0, according to the readme attached to the file.
  # However, 
  # spf <- sampleRegular(pf, size = 1e6)
  # spf %>% as.vector() %>%  unique() 
  # shows that 0 is a value from the raster. 
  
  ### Crop to gfc-island extent first
  aoi_sp_prj <- spTransform(aoi_sp, crs(pf))

  crop(pf, y = aoi_sp_prj,
       filename = file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,".tif")),
       datatype = "INT1U",
       overwrite = TRUE)
  
  # no need to expend here, as for ioppm, because pf already completely covers aoi_sp_prj *
  rm(aoi_sp_prj)
  
  ### Reclassify primary forest map
  
  # read the island-croped pf map 
  pf <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,".tif")))
  
  # The values in Margono et al. data have the following interpretation:
  # 0   - Out of area study
  # 1   - No change of primary degraded forest from 2000-2012 
  # 2   - No change of primary intact forest from 2000-2012
  # 3   - No change of non-primary from 2000-2012
  # 4   - Primary intact, cleared 2005
  # 5   - Primary intact, cleared 2010
  # 6   - Primary intact, cleared 2012 
  # 7   - Primary intact, degraded 2005
  # 8   - Primary intact, degraded 2010
  # 9   - Primary intact, degraded 2012
  # 10 - Primary degraded, cleared 2005
  # 11 - Primary degraded, cleared 2010 
  # 12 - Primary degraded, cleared 2012 
  # 13 - Primary intact degraded 2005, cleared 2010
  # 14 - Primary intact degraded 2005, cleared 2012
  # 15 - Primary intact degraded 2010, cleared 2012
  # 
  # Here, we do not care when each type of forest was degraded. 
  # The only classes that are not primary forest in 2000 are
  # 0 - Out of area study 
  # 3   - No change of non-primary from 2000-2012
  
  # The intact classes are 2, 4-9, 13-15
  # The degraded classes are 10-12
  
  # Therefore the reclassification matrix is 
  m <- c(0,0,0, 
         0,1,2,
         1,2,1,
         2,3,0, 
         3,9,1,
         9,12,2,
         12,15,1)
  rclmat <- matrix(m, ncol = 3, byrow = TRUE)
  
  # Reclassify primary forest classes in a cluster framework
  # The two first columns of rclmat are intervals, values within which should be converted in the third column's value.
  # The right = TRUE (the default) means that intervals are open on the left and closed on the right. 
  # include.lowest means that the lowest interval is closed on the left. 
  
  beginCluster() # this uses by default detectCores() - 1
  
  clusterR(pf,
           fun = reclassify,
           args = list(rcl = rclmat, right = T, include.lowest = TRUE), 
           filename = file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_reclassified.tif")), 
           datatype = "INT1U",
           overwrite = TRUE)
  
  endCluster()

  
  ### Align to GFC loss crs, resolution and exact extent.  
  
  # Read the pf map (with roughly the island extent and reclassified) 
  pf <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_reclassified.tif")))
  
  # Read the target GFC loss layer
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
 
  beginCluster() # this uses by default detectCores() - 1
  
  projectRaster(from = pf, to = loss,
                method = "ngb",
                filename = file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif")),
                datatype = "INT1U",
                overwrite = TRUE)
  endCluster()
  
  rm(pf, loss, m, rclmat)
  removeTmpFiles(h=0)
  
  
  print(paste0("complete", "temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif"))
  


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Overlay forest loss and oil palm plantation maps ####
  
  # We want to keep forest loss pixels only within 2015 plantations in order to induce forest conversion to plantation,
  # BUT within primary forest only. In particular, in order not to count plantation renewals as forest conversion to plantation.
  # po maps are binary with 1 meaning plantation in 2015
  
  # We repeat this for three primary forest types: intact only, degraded only, and intact or degraded. 
  # Thus, we define three overlay functions.
  overlay_intact   <- function(rs){rs[[1]]*rs[[2]]*(rs[[3]] == 1)}
  overlay_degraded <- function(rs){rs[[1]]*rs[[2]]*(rs[[3]] == 2)}
  overlay_total    <- function(rs){rs[[1]]*rs[[2]]*(rs[[3]] != 0)}
  # multiplies a cell of forest loss (rs[[1]]) by 0 if it is not a plantation in 2015 (rs[[2]])
  # or if it something else than intact primary forest in 2000 (overlay_intact)
  # or if it something else than degraded primary forest in 2000 (overlay_degraded)
  # or if it something else than either intact or degraded primary forest in 2000 (overlay_total)
  
  ## Read necessary layers and stack them 
  
  # GFC loss layer (rs[[1]])
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
  
  # plantations (rs[[2]])
  ioppm2015 <- raster(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,"_aligned.tif"))
  
  # primary forest (rs[[3]])
  pf <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif")))
  
  # stack is necessary for clusterR
  rs <- stack(loss, ioppm2015, pf)
  
  # note that using calc below is equivalent to using overlay but more appropriate to the input being a stack, 
  # which is necessary to pass several raster layers to the first argument of clusterR
  
  # run the computation in parallel with clusterR, as cells are processed one by one independently.
  beginCluster() # uses by default detectedCores() - 1
  
  # For intact primary forest
  clusterR(rs,
           fun = calc, #
           args = list(overlay_intact),
           filename = file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_intact.tif")),
           datatype = "INT1U",
           overwrite = TRUE )
  
  # For degraded primary forest
  clusterR(rs,
           fun = calc, 
           args = list(overlay_degraded),
           filename = file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_degraded.tif")),
           datatype = "INT1U",
           overwrite = TRUE )
  
  # For total primary forest
  clusterR(rs,
           fun = calc, 
           args = list(overlay_total),
           filename = file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_total.tif")),
           datatype = "INT1U",
           overwrite = TRUE )
  
  endCluster()
  
  rm(loss, ioppm2015, pf, overlay_intact, overlay_degraded, overlay_total)
  removeTmpFiles(h=0)
    

  print(paste0("complete", "temp_data/processed_lu/lucpfip_",island,"_total.tif"))
  
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Split LUCFP to annual layers ####
  
  ### Function description
  # parallel_split has for input a single lucpfip layer where each pixel has a value corresponding to the year when a lucpfip event occured;
  # it outputs annual layers in each of which pixels are either 1 if a lucpfip event occured that year, and 0 else.
  # the tasks are year specific and independent across years, therefore they are executed parallely over years.
  parallel_split <- function(pf_type, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the forest type level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the forest type level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_split <- function(time){
      # define process (island, pf_type) we are in 
      process <- file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_",pf_type,".tif"))
      
      # #set temp directory
      dir.create(paste0(process,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(process,"_Tmp")))
      
      # read in the input
      lucpfip_prj <- raster(process)
      
      # define output file name 
      output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",island,"_",pf_type,"_", years[time],".tif"))
      
      # split it into annual binary layers
      calc(lucpfip_prj,
           fun = function(x){if_else(x == time, true = 1, false = 0)},
           filename = output_filename,
           datatype = "INT1U",
           overwrite = TRUE )
      # remove process temporary files
      unlink(file.path(paste0(process,"_Tmp")), recursive = TRUE)
    }
    
    ## register cluster
    registerDoParallel(cores = ncores)
    
    ## define foreach object.
    foreach(t = 1:length(years),
            # .combine combine the outputs as a mere character list (by default)
            .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
            .multicombine = TRUE,
            .export = c("island"),
            .packages = c("dplyr", "raster", "rgdal")
    ) %dopar%  annual_split(time = t)
  }
  
  ### Execute it for each primary forest type

  pf_typeS <- c("intact", "degraded", "total")
  for(pf_type in pf_typeS){
    
    parallel_split(pf_type = pf_type, detectCores() - 1) # ~500 seconds / annual layer
    
    removeTmpFiles(h=0)
  }  
    
  rasterOptions(tmpdir = "temp_data/raster_tmp")  
  
  return(print("complete prepare_lucpfip ", island))
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 2. AGGREGATE THE PIXELS TO A GIVEN PARCEL SIZE. #####

aggregate_lucpfip <- function(island, parcel_size){
  
  ### Function description
  # The function has for inputs annual layers of lucfp events at the pixel level.
  # It aggregates these pixels to a parcel size defined by parcel_size (in meters).
  # The aggregation operation is the sum of the pixel lucfp events.
  # Each annual aggregation is tasked in parallel.
  parallel_aggregate <- function(pf_type, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the pf_type level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the pf_type level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_aggregate <- function(time){
      # Define which process (island, pf_type, and year) we are in:
      processname <- file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",island,"_",pf_type,"_", years[time],".tif"))
      
      #set temp directory
      dir.create(paste0(processname,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(processname,"_Tmp")))
      
      # read in the input.
      lucpfip_annual <- raster(processname)
      
      # define output file name
      output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,"_",years[time],".tif"))
      
      # aggregate it from the ~30m cells to parcel_size cells with mean function.
      raster::aggregate(lucpfip_annual, fact = c(parcel_size/res(lucpfip_annual)[1], parcel_size/res(lucpfip_annual)[2]),
                        expand = FALSE,
                        fun = sum,
                        na.rm = FALSE, # NA cells are in margins, see the NOTES part. If FALSE, aggregations at margins that use NA 
                        # are discarded because the sum would be spurious as it would count all NA as 0s while it is not necessary the case.
                        filename = output_filename,
                        datatype = "INT4U", # because the sum may go up to ~ 10 000 with parcel_size = 3000,
                        # but to more than 65k with parcel_size = 10000 so INT4U will be necessary;
                        overwrite = TRUE)
      #removes entire temp directory without affecting other running processes (but there should be no temp file now)
      unlink(file.path(paste0(processname,"_Tmp")), recursive = TRUE)
      #unlink(file.path(tmpDir()), recursive = TRUE)
      # return the path to this parcels file
      #return(output_filename)
    }
    
    ## register cluster
    registerDoParallel(cores = ncores)
    
    ##  define foreach object.
    foreach(t = 1:length(years),
            # .combine combine the outputs as a mere character list (by default)
            .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
            .multicombine = TRUE,
            .export = c("island", "parcel_size"),
            .packages = c("raster", "rgdal")
    ) %dopar% annual_aggregate(time = t)
  }

    
  ### Execute the function to compute the RasterBrick object of 18 annual layers for each primary forest type
  
  pf_typeS <- c("intact", "degraded", "total")
  for(pf_type in pf_typeS){
    # run the computation, that writes the layers 
    parallel_aggregate(pf_type = pf_type, ncores = detectCores() - 1)
    
    # brick the layers together and write the brick
    rasterlist <- list.files(path = "temp_data/processed_lu/annual_maps", 
                             pattern = paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,"_"), 
                             full.names = TRUE) %>% as.list()
    
    parcels_brick <- brick(rasterlist)
    
    writeRaster(parcels_brick,
                filename = file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,".tif")),
                datatype = "INT4U",
                overwrite = TRUE)
    
    rm(rasterlist, parcels_brick)
    removeTmpFiles(h=0)
  }
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")    
  
  print(paste0("complete aggregate_lucpfip ",island," ",parcel_size/1000, "km"))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 3. CONVERT TO A DATAFRAME THE PARCELS WITHIN A GIVEN CATCHMENT RADIUS  #####

to_panel_within_CR <- function(island, parcel_size, catchment_radius){
  
  ### Function description
  # raster_to_df converts the raster bricks of annual layers of parcels to a panel dataframe.
  # This is executed for each pf_type, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment radius.
  # 3. reshaping the values in these parcels to a long format panel dataframe
  raster_to_df <- function(pf_type){
    
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## 1. Masking.
    # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
    
    # Make the mask
    mills <- read.dta13(file.path("temp_data/processed_mill_geolocalization/IBS_UML_cs.dta"))
    mills <- mills[mills$island_name == island,]
    #turn into an sf object.
    mills <- st_as_sf(mills,	coords	=	c("lon",	"lat"), crs=4326)
    # keep only the geometry, we do not need mills attributes here.
    mills <- st_geometry(mills)
    # set CRS and project
    mills_prj <- st_transform(mills, crs = indonesian_crs)
    #define big catchment areas to have a large AOI.
    mills_ca <- st_buffer(mills_prj, dist = 60000)
    # work with squares rather than with circles
    for(i in 1:length(mills_ca)){
      mills_ca[i] <- st_as_sfc(st_bbox(mills_ca[i]))
    }
    total_ca <- st_union(st_geometry(mills_ca))
    # coerce to a SpatialPolygon
    total_ca_sp <- as(total_ca, "Spatial")
    # keep mills_prj we need it below
    rm(total_ca, mills_ca, mills)
    
    # Mask
    parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type)
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
    
    mask(x = parcels_brick, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(parcels_brick, total_ca_sp)
    
    
    ## 2. Selecting parcels within a given distance to a mill at least one year
    # (i.e. the parcel is present in the dataframe in all years even if it is within say 50km of a mill only since 2014)
    
    # Turn the masked raster to a sf dataframe
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_masked.tif")))
    
    m.df_wide <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)

    m.df_wide <- m.df_wide %>% dplyr::rename(lon = x, lat = y)
    m.df_wide <- st_as_sf(m.df_wide, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
    
    # Remove here parcels that are not within the catchment area of a given size (defined by catchment radius)
    # coordinates of all mills (crs is indonesian crs, unit is meter)
    within <- st_is_within_distance(m.df_wide, mills_prj, dist = catchment_radius)
    m.df_wide <- m.df_wide %>% dplyr::filter(lengths(within) >0)
    m.df_wide <- m.df_wide %>% st_drop_geometry()
    
    rm(within, parcels_brick)
    
    
    ## 3. Reshaping to long format
    # make parcel id
    island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
    m.df_wide$parcel_id <- paste0(island_id, c(1:nrow(m.df_wide))) %>% as.numeric()
    
    # vector of the names in the wide format of our time varying variables
    # the column names are the layer names in the parcels_brick + the layer index, separated by "." see examples in raster::as.data.frame
    varying_vars <- paste0(parcels_brick_name, "_masked.", seq(from = 1, to = 18))
    
    # reshape to long
    m.df <- stats::reshape(m.df_wide,
                           varying = varying_vars,
                           v.names = paste0("lucpfip_pixelcount_",pf_type),
                           sep = ".",
                           timevar = "year",
                           idvar = c("parcel_id"), # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                           ids = "parcel_id",
                           direction = "long",
                           new.row.names = seq(from = 1, to = nrow(m.df_wide)*length(years), by = 1))
    
    rm(varying_vars, m.df_wide)
    # replace the indices from the raster::as.data.frame with actual years.
    m.df <- mutate(m.df, year = years[year])
   
    m.df <- setorder(m.df, parcel_id, year)
    saveRDS(m.df,
            file = file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR_",
                                    pf_type,".rds")))
  }

  ### Execute it
  pf_typeS <- c("intact", "degraded", "total")
  for(pf_type in pf_typeS){
    
    raster_to_df(pf_type = pf_type)
    
    removeTmpFiles(h=0)
  }  

  
  print(paste0("complete to_panel_within_CR ",island," ",parcel_size/1000,"km ",catchment_radius/1000,"CR"))

}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### EXECUTE FUNCTIONS AND MERGE THE OUTPUTS #####

#### Execute the functions ####
# Only if their outputs have not been already computed

### Prepare a 30m pixel map of lucfp for each Island
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",Island,"_total_2018.tif")))){
    
    prepare_pixel_lucpfip(Island)
  }
}

### Aggregate this Island map to a chosen parcel size (3km, 6km and 9km for instance)
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",Island,"_",PS/1000,"km_total.tif")))){
    
    aggregate_lucpfip(island = Island,
                      parcel_size = PS)
  }
}

### For that Island and for each aggregation factor, extract panels of parcels within different catchment area sizes 
# (radius of 10km, 30km and 50km)
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
    if(!file.exists(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"CR_total.rds")))){
      
      to_panel_within_CR(island = Island,
                         parcel_size = PS,
                         catchment_radius = CR)
      
    } # only the function execution is conditioned to the file existance, not the loop incrementation
    
    CR <- CR + 20000
  }
}

#### Gather the lucfp variables for each parcel_size and catchment_radius combinations. ####
PS <- 3000  

CR <- 10000 # i.e. 10km radius
while(CR < 60000){
  
  # For each Island, join columns of lucfp variable for different forest thresholds. 
  df_list <- list()
  IslandS <- c("Sumatra", "Kalimantan", "Papua")
  for(Island in IslandS){

    df_intact   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"CR_intact.rds")))
    df_degraded <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"CR_degraded.rds")))
    df_total    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"CR_total.rds")))
    
    df_degraded <- dplyr::select(df_degraded, -lon, -lat)
    df <- inner_join(df_intact, df_degraded, by = c("parcel_id", "year"))
    
    df_total <- dplyr::select(df_total, -lon, -lat)
    df_list[[match(Island, IslandS)]] <- inner_join(df, df_total, by = c("parcel_id", "year"))
  }
  
  # stack the three Islands together
  indo_df <- bind_rows(df_list)
  
  saveRDS(indo_df, file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",PS/1000,"km_",CR/1000,"CR.rds")))
  
  rm(indo_df, df_list)
  CR <- CR + 20000
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### NOTES #####
# Not necessarily relevant 

# note on rm in function: removes object in the function "frame" i.e. environment without
# having to specify something (I tested it)


### ON EXTRACT_GFC
# extract télécharge les tiles qui couvrent notre AOI
# et ensuite pour ces tiles là il n'y a aucun NA, dans le résultat du extract, même
# pour les pixels en dehors de l'AOI. En revanche il y a des NAs
# dans la partie du raster couverte pas un tile qui ne couvre pas l'AOI.
# donc si on prend un bbox comme aoi on n'a pas de NA.

# to extract and project in the same time (AOI_sp should be projected) but did not work, I don't know why.
# extract_gfc(AOI_sp, data_folder, stack = "change", to_UTM = TRUE,
# dataset = "GFC-2017-v1.5", filename = "gfc_data_prj.tif", overwrite = TRUE)
# defo_e <- raster("gfc_data_prj.tif")
# crs(defo_e)

#   We don't use gfc_stats because we want to keep information at the pixel level and not at the aoi's in order 
#   to overlay it with plantations.




### GFC mask 
# we do not mask with water mask from gfc because there are many other reasons why 
# deforestation cannot occur (cities, mountains) and these are taken into account in FE 
# and the distributional effect of this on our outcome variabel will be captured 
# by the poisson model. 

# Note also that For Sumatra at least, the northern part of the most north CA is not covered 
# by Austin data. Those NAs are simply not converted to the dataframe. This is not an 
# issue because the unit of observation is not the CA (one of which would not be the same size)
# but the parcel. In other words we don't bother that some CA are not fully observed in our maps. 




### NAs
# The forest loss and aligned plantations maps don't have NAs
# But once overlaid, some NAs are produced at the margins (because they don't perfectly align)
# These NAs are not removed (trimmed) or reclassified. 
# They disappear when 
#lucfp30 <- raster("lucfp_30th.tif")
# clusterR(lucfp30, 
#          fun = reclassify, 
#          args = list(rcl = cbind(NA,0)),
#          filename = "lucfp_30th.tif",
#          datatype = "INT1U",
#          overwrite = TRUE )




#  brick after a normal aggregation (no foreach)
# # read in the 18 layers and brick them
# layers_paths <- list.files(path = "./annual_parcels", pattern = paste0("parcels_",PS/1000,"km_",threshold, "th_"), full.names = TRUE) %>% as.list()
#   #layers_paths <- list.files(path = "./annual_maps", pattern = paste0("defo_",threshold, "th_"), full.names = TRUE) %>% as.list()
#   parcels_brick <- layers_paths %>% brick()
#   
#   #write the brick
#   writeRaster(parcels_brick, 
#               filename = paste0("./bricked_parcels/parcels_",PS/1000,"km_",threshold,"th.tif"), 
#               overwrite = TRUE)
#   rm(parcels_brick)



### rationale if aggregation is with mean: 
# Since in each annual map, original pixels are either 1 or 0 valued (conversion from forest to op plantation remotely sensed 
# for that pixel that year or not) and each pixel is the same area, the mean value of a group of pixels is 
# the ratio of the area deforested over the parcel area. This is refered to as the percentage of deforestation. 



### MASK ? is not useful because reclassifying NAs to 0s does not make the file lighter (NA is stored in a large value 
# see NAvalue() and ?datatype)

### RECLASSIFY PRIMARY FOREST  
# sample testing
# pf <- raster(file.path("input_data/margono_primary_forest/timeseq_change00_12.tif"))
# 
# provinces <- st_read("C:/Users/GUYE/Desktop/opalval/analysis/input/IDN_adm/IDN_adm1.shp")
# provinces <- dplyr::select(provinces, NAME_1)
# # papua sample 
# papua_bbox <- st_bbox(provinces[provinces$NAME_1 == "Papua" |
#                                    provinces$NAME_1 == "Irian Jaya Barat", ]) %>% st_as_sfc()
# papua_bbox <- st_transform(papua_bbox, crs = crs(pf))
# papua_bbox <- as(papua_bbox, "Spatial")
# 
# # Riau sample 
# riau_bbox <- st_bbox(provinces[provinces$NAME_1 == "Riau",]) %>% st_as_sfc()
# riau_bbox <- st_transform(riau_bbox, crs = crs(pf))
# riau_bbox <- as(riau_bbox, "Spatial")
# 
# dir.create("temp_data/test")
# 
# 
# crop(pf, papua_bbox, filename = "temp_data/test/papua_margono_primary_forest.tif",
#      overwrite = TRUE,
#      progress = "text")
# 
# crop(pf, riau_bbox, filename = "temp_data/test/riau_margono_primary_forest.tif",
#      overwrite = TRUE,
#      progress = "text")
# 
# papuapf.tif <- raster("temp_data/test/papua_margono_primary_forest.tif")
# papuapf.gri <- raster("temp_data/test/papua_margono_primary_forest.gri")
# 
# riaupf.tif <- raster("temp_data/test/riau_margono_primary_forest.tif")
# riaupf.gri <- raster("temp_data/test/riau_margono_primary_forest.gri")
# 
# riaupf.tif %>% getValues %>% unique()
# riaupf.gri %>% getValues %>% unique()
# 
# plot(papuapf)
# rasterOptions(maxmemory = 5e9)
# riaupf <- raster("temp_data/test/riau_margono_primary_forest")
# 
# riaupf %>% getValues() %>% unique()
# 
# spf <- sampleRegular(pf, size = 1e6)
# 
# spf %>% as.vector() %>%  unique() %>% length()
# sriaupf <- sampleRegular(riaupf, size = 1e6)
# 
# spf %>% class()
# length(spf)
# 
# 
# # test reclassify
# r <- raster(ncols = 5, nrow = 3)
# values(r) <- 0:14
# m <- c(0,0,0,
#        0,1,2,
#        1,2,1,
#        2,3,0,
#        3,9,1,
#        9,12,2,
#        12,14,1)
# rclmat <- matrix(m, ncol = 3, byrow = TRUE)
# rr <- reclassify(r, rclmat, right = TRUE, include.lowest = T)
# 
# # cluster reclassify on Papua
# beginCluster()
# clusterR(papuapf,
#          fun = reclassify,
#          args = list(rcl = rclmat, right = T, include.lowest = TRUE),
#          filename = "temp_data/test/papua_reclassified_pf.tif",
#          datatype = "INT1U",
#          overwrite = TRUE,
#          progress = "window")
# 
# endCluster()
# 
# riaupf_rec <- raster("temp_data/test/papua_reclassified_pf.tif")
# 
# 

### Frmaework foreach dopar for projectRaster lucpfip

# ### Function description
# # Project to Indonesian crs, in parallel, for each threshold definition. 
# parallel_project_lucfp <- function(ncores){
#   
#   ## sequence over which to execute the task
#   lucpfip_files <-list.files(path = file.path("temp_data/processed_lu"), 
#                              pattern = paste0("lucpfip_",island,"_"))
#   
#   ## read the input to the task
#   # it is task specific
#   
#   ## define the task
#   project_lucpfip <- function(lucpfip_file_index){
#     lucpfip <- raster(lucpfip_files[[lucpfip_file_index]])
#     
#     projectRaster(from = lucpfip,
#                   crs = indonesian_crs,
#                   method = "ngb",
#                   filename = file.path(paste0("temp_data/processed_lu/",names(lucpfip),"_prj.tif")),
#                   datatype = "INT1U",
#                   overwrite = TRUE )
#   }
#   
#   ## register cluster
#   registerDoParallel(cores = ncores)
#   
#   ## define foreach object
#   foreach(i = 1:length(lucpfip_files),
#           # .combine combine the outputs as a mere character list (by default)
#           .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
#           .multicombine = TRUE,
#           .export = c("island", "indonesian_crs"),
#           .packages = c("raster", "rgdal")
#   ) %dopar% project_lucfp(lucpfip_file_index = i)
# }
# 
# ### Execute it
# parallel_project_lucfp(detectCores() - 1)



### Project lucpfip

# 
# # This is necessary because we will need to make computations on this map within mills' catchment *areas*.
# # If one does not project this map, then catchment areas all have different areas while being defined with a common buffer.
# 
# # Read in rasters to reproject
# lucpfip_intact <- raster(file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_intact.tif")))
# lucpfip_degraded <- raster(file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_degraded.tif")))
# lucpfip_total <- raster(file.path(paste0("temp_data/processed_lu/lucpfip_",island,"_total.tif")))
# 
# beginCluster()
# 
# projectRaster(from = lucpfip_intact,
#               crs = indonesian_crs,
#               method = "ngb",
#               filename = file.path(paste0("temp_data/processed_lu/",names(lucpfip_intact),"_prj.tif")),
#               datatype = "INT1U",
#               overwrite = TRUE )
# 
# projectRaster(from = lucpfip_degraded,
#               crs = indonesian_crs,
#               method = "ngb",
#               filename = file.path(paste0("temp_data/processed_lu/",names(lucpfip_degraded),"_prj.tif")),
#               datatype = "INT1U",
#               overwrite = TRUE )
# 
# projectRaster(from = lucpfip_total,
#               crs = indonesian_crs,
#               method = "ngb",
#               filename = file.path(paste0("temp_data/processed_lu/",names(lucpfip_total),"_prj.tif")),
#               datatype = "INT1U",
#               overwrite = TRUE )
# 
# endCluster()
# 
# rm(lucpfip_intact, lucpfip_degraded, lucpfip_total)
# removeTmpFiles(h=0)