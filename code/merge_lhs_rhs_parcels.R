### Merge panel data frames of parcels with LHS and RHS variables 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "foreign")

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
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# There are several merges to do because several data frames
# Each data frame is a different set of parcels, depending on parcel size and catchment radius. 

# For each set of parcels, there are also different sets of outcome variables 
# (lucfip, lucpfip, emissions, lucpfsmp...)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


merge_lhs_rhs <- function(outcome_variables, parcel_size, catchment_radius){

#  outcome variables 
LHS <- readRDS(file.path(paste0("temp_data/processed_parcels/",
                                outcome_variables,"_panel_",
                                parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))

# keep only year before 2015 (after they mean nothing since we plantation data are from 2015)
LHS <- LHS[LHS$year<=2015,] # now runs from 2001-1998
# remove coordinates, they are already in RHS
LHS <- dplyr::select(LHS, -lat, -lon)

# explicative variables (runs from 1998-2015)
RHS <-  readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_reachable_uml_",
                                  parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))

# MERGE
# years 1998 - 2000 will thus not match, but we want to keep them, hence all = TRUE
final <- base::merge(LHS, RHS, by = c("parcel_id", "year"), all = TRUE)  

# some arrangements
final <- dplyr::arrange(final, parcel_id, year)
row.names(final) <- seq(1,nrow(final))

saveRDS(final, file.path(paste0("temp_data/panel_parcels_final_",
                                outcome_variables,"_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))  

write.dta(final, file.path(paste0("temp_data/panel_parcels_final_",
                                outcome_variables,"_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.dta")))  

# this does not write the largest data frames (catchment_radius of 50km)
# (Erreur : Error in libxlsxwriter: 'Worksheet row or column index out of range.')
# run it to get the two other ones in xlsx format still. 
write_xlsx(final, file.path(paste0("temp_data/panel_parcels_final_",
                                  outcome_variables,"_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"CR.xlsx")))  

}

OV <- "lucpfip"
PS <- 3000 
catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
for(CR in catchment_radiuseS){
  merge_lhs_rhs(outcome_variables = OV, 
                parcel_size = PS, 
                catchment_radius = CR)
}


#### TESTING ZONE #### 
# names(LHS)  
# names(RHS)  
# 
# length(unique(LHS$parcel_id))
# length(unique(RHS_2001$parcel_id))
# 
# # LHS was not ordered as expected, hence 
# RHS_2001 <- RHS[RHS$year >= 2001,]
# all.equal(RHS_2001$parcel_id, LHS$parcel_id) # returns FALSE
# 
# LHS_ordered <- dplyr::arrange(LHS, parcel_id, year)
# all.equal(LHS_ordered, LHS)
# # RHS was ordered as expected
# RHS_ordered <- dplyr::arrange(RHS, parcel_id, year)
# all.equal(RHS_ordered, RHS)
# 
# # once reordered, we indeed have the same set of parcels in both data frames. 
# all.equal(RHS_2001$parcel_id, LHS_ordered$parcel_id)
# all.equal(RHS_2001[,c("lat", "lon")], LHS_ordered[,c("lat", "lon")]) 
# coordinates match too (only the row names are different). 



