
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("plyr", "dplyr", "data.table", 
                   "foreign", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "fixest")
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
# troublePackages <- c() 
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 

### NEW FOLDERS USED IN THIS SCRIPT 



### Set dictionary for names of variables to display in regression tables 
setFixest_dict()



catchment_radius <- 3e4
island <- c("Sumatra",",", "Kalimantan")
outcome_variable <- "lucfip_pixelcount_30th"

compare_estimators <- function(catchment_radius, island, outcome_variable){
  
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))
  names(d)
  
  length(unique(d$parcel_id)) # length cross sections: 47648
  unique(d$year) # years from 2001 to 2015
  
  ## Island subsampling
  d <- d[d$island %in% island,]
  
  ## Base LU subsampling
  # df30 <- d[d$any_fc2000_30th,]
  # length(unique(df30$parcel_id)) # 41844
  # dpf <- d[d$any_pfc2000_total,]
  # length(unique(dpf$parcel_id)) # 24959 
  if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
    d <- d[d$any_pfc2000_total,]
  }
  if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
     outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
     outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
    d <- d[d$any_fc2000_30th,]
  }

!!as.symbol(outcome_variable)

poiglm_est <- fixest::feglm(lucfip_pixelcount_30th~wa_ffb_price_imp1_yoyg_3pya | parcel_id, data = d, 
                        family = "poisson")

qpoiglm_est <- fixest::feglm(lucfip_pixelcount_30th~wa_ffb_price_imp1_yoyg_3pya | parcel_id, data = d, 
                       family = "quasipoisson")

poimlm_est <- fixest::femlm(lucfip_pixelcount_30th~wa_ffb_price_imp1_yoyg_3pya | parcel_id, data = d, 
                       family = "poisson")

nb_est <- fixest::fenegbin(lucfip_pixelcount_30th~wa_ffb_price_imp1_yoyg_3pya | parcel_id, data = d)

gauglm_est <- fixest::feglm(lucfip_pixelcount_30th~wa_ffb_price_imp1_yoyg_3pya | parcel_id, data = d, 
                       family = "gaussian")


return(etable(poiglm_est, qpoiglm_est, poimlm_est, nb_est, gauglm_est, 
       se = "standard", 
       tex = TRUE,
       file = file.path(paste0("outputs/regressions/tables/est_comparisons_",
                               island,"_",
                               parcel_size/1000,"km_",
                               catchment_radius/1000,"km_",
                               outcome_variable)), 
       title = paste0("Estimator comparisons, ",island, 
                      " catchment radius of ", catchment_radius/1000,"km, "),
       subtitles = c("Poisson GLM", "Quasi-Poisson GLM", "Poisson MLM", "Negative Binomial", "Gaussian GLM"),
       coefstat = "confint",
       sdBelow = TRUE,
       dict = TRUE)  
)

  
  
  
  
  
  d$wa_ffb_price_imp1_yoyg_3pya
  
  
  
  
  
}



# # give it a panel class for lag analyses
# pd <- panel(d, panel.id = c("parcel_id", "year"))

# on the sample of parcels covered with a positive extent of forest of 30% tree canopy density in 2000.
# using the area in hectare of LUCFIP on forest of 30% tree canopy density (outside industrial plantations in 2000)
est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1), data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+district, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+district+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island^year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+district^year, data = df30_01)
etable(est_f30)



### CPO
# on the sample of parcels covered with a positive extent of forest of 30% tree canopy density in 2000.
# using the area in hectare of LUCFIP on forest of 30% tree canopy density (outside industrial plantations in 2000)
est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1), data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+district, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+district+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island^year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+district^year, data = df30_01)
etable(est_f30)













