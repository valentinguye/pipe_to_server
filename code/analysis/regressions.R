
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table", 
                   "foreign", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "fixest", 
                   "modelsummary", 
                   "ggplot2")
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

### PARCEL SIZE
parcel_size <- 3000


### Set dictionary for names of variables to display in regression tables 
setFixest_dict(c(parcel_id = "grid cell",
                 lucpfip_ha_total = "LUCPFIP (ha)", 
                 lucpfip_pixelcount_total = "LUCPFIP (pixels)", 
                 lucfip_ha_30th = "LUCFIP (30 pct. canopy density, ha)",
                 lucfip_ha_60th = "LUCFIP (60 pct. canopy density, ha)",
                 lucfip_ha_90th = "LUCFIP (90 pct. canopy density, ha)",
                 lucfip_pixelcount_30th = "LUCFIP (30 pct. canopy density, pixels)",
                 lucfip_pixelcount_60th = "LUCFIP (60 pct. canopy density, pixels)",
                 lucfip_pixelcount_90th = "LUCFIP (90 pct. canopy density, pixels)",
                 # No time dynamics FFB variables 
                 wa_ffb_price_imp1_3ya = "FFB price signal, 3 year average",
                 wa_ffb_price_imp1_3ya_lag1 = "FFB price signal, 3 year average (lagged)",
                 wa_ffb_price_imp1_4ya = "FFB price signal, 4 year average",
                 wa_ffb_price_imp1_4ya_lag1 = "FFB price signal, 4 year average (lagged)",
                 wa_ffb_price_imp1_5ya = "FFB price signal, 5 year average",
                 wa_ffb_price_imp1_5ya_lag1 = "FFB price signal, 5 year average (lagged)",
                 wa_ffb_price_imp1_yoyg_3ya = "FFB price signal y-o-y growth rate, 3 year average",
                 wa_ffb_price_imp1_yoyg_3ya_lag1 = "FFB price signal y-o-y growth rate, 3 year average (lagged)",
                 wa_ffb_price_imp1_yoyg_4ya = "FFB price signal y-o-y growth rate, 4 year average",
                 wa_ffb_price_imp1_yoyg_4ya_lag1 = "FFB price signal y-o-y growth rate, 4 year average (lagged)",
                 wa_ffb_price_imp1_yoyg_5ya = "FFB price signal y-o-y growth rate, 5 year average",
                 wa_ffb_price_imp1_yoyg_5ya_lag1 = "FFB price signal y-o-y growth rate, 5 year average (lagged)",
                 # No time dynamics CPO variables 
                 wa_cpo_price_imp1_3ya = "CPO price signal, 3 year average",
                 wa_cpo_price_imp1_3ya_lag1 = "CPO price signal, 3 year average (lagged)",
                 wa_cpo_price_imp1_4ya = "CPO price signal, 4 year average",
                 wa_cpo_price_imp1_4ya_lag1 = "CPO price signal, 4 year average (lagged)",
                 wa_cpo_price_imp1_5ya = "CPO price signal, 5 year average",
                 wa_cpo_price_imp1_5ya_lag1 = "CPO price signal, 5 year average (lagged)",
                 wa_cpo_price_imp1_yoyg_3ya = "CPO price signal y-o-y growth rate, 3 year average",
                 wa_cpo_price_imp1_yoyg_3ya_lag1 = "CPO price signal y-o-y growth rate, 3 year average (lagged)",
                 wa_cpo_price_imp1_yoyg_4ya = "CPO price signal y-o-y growth rate, 4 year average",
                 wa_cpo_price_imp1_yoyg_4ya_lag1 = "CPO price signal y-o-y growth rate, 4 year average (lagged)",
                 wa_cpo_price_imp1_yoyg_5ya = "CPO price signal y-o-y growth rate, 5 year average",
                 wa_cpo_price_imp1_yoyg_5ya_lag1 = "CPO price signal y-o-y growth rate, 5 year average (lagged)",
                 # SR FFB variables
                 wa_ffb_price_imp1 = "FFB price signal",
                 wa_ffb_price_imp1_lag1 = "FFB price signal (lagged)",
                 wa_ffb_price_imp1_yoyg = "FFB price signal y-o-y growth rate",
                 wa_ffb_price_imp1_yoyg_lag1 = "FFB price signal y-o-y growth rate (lagged)",
                 wa_ffb_price_imp1_dev_2pya = "FFB price signal deviation from 2 past year average",
                 wa_ffb_price_imp1_dev_2pya_lag1 = "FFB price signal deviation from 2 past year average (lagged)",
                 wa_ffb_price_imp1_dev_3pya = "FFB price signal deviation from 3 past year average",
                 wa_ffb_price_imp1_dev_3pya_lag1 = "FFB price signal deviation from 3 past year average (lagged)",
                 wa_ffb_price_imp1_dev_4pya = "FFB price signal deviation from 4 past year average",
                 wa_ffb_price_imp1_dev_4pya_lag1 = "FFB price signal deviation from 4 past year average (lagged)",
                 # LR FFB variables
                 wa_ffb_price_imp1_2pya = "FFB price signal, 2 past year average",
                 wa_ffb_price_imp1_2pya_lag1 = "FFB price signal, 2 past year average (lagged)",
                 wa_ffb_price_imp1_yoyg_2pya = "FFB price signal y-o-y growth rate, 2 past year average",
                 wa_ffb_price_imp1_yoyg_2pya_lag1 = "FFB price signal y-o-y growth rate, 2 past year average (lagged)",
                 wa_ffb_price_imp1_3pya = "FFB price signal, 3 past year average",
                 wa_ffb_price_imp1_3pya_lag1 = "FFB price signal, 3 past year average (lagged)",
                 wa_ffb_price_imp1_yoyg_3pya = "FFB price signal y-o-y growth rate, 3 past year average",
                 wa_ffb_price_imp1_yoyg_3pya_lag1 = "FFB price signal y-o-y growth rate, 3 past year average (lagged)",
                 wa_ffb_price_imp1_4pya = "FFB price signal, 4 past year average",
                 wa_ffb_price_imp1_4pya_lag1 = "FFB price signal, 4 past year average (lagged)",
                 wa_ffb_price_imp1_yoyg_4pya = "FFB price signal y-o-y growth rate, 4 past year average",
                 wa_ffb_price_imp1_yoyg_4pya_lag1 = "FFB price signal y-o-y growth rate, 4 past year average (lagged)",
                 # SR CPO variables
                 wa_cpo_price_imp1 = "CPO price signal",
                 wa_cpo_price_imp1_lag1 = "CPO price signal (lagged)",
                 wa_cpo_price_imp1_yoyg = "CPO price signal y-o-y growth rate",
                 wa_cpo_price_imp1_yoyg_lag1 = "CPO price signal y-o-y growth rate (lagged)",
                 wa_cpo_price_imp1_dev_2pya = "CPO price signal deviation from 2 past year average",
                 wa_cpo_price_imp1_dev_2pya_lag1 = "CPO price signal deviation from 2 past year average (lagged)", 
                 wa_cpo_price_imp1_dev_3pya = "CPO price signal deviation from 3 past year average",
                 wa_cpo_price_imp1_dev_3pya_lag1 = "CPO price signal deviation from 3 past year average (lagged)", 
                 wa_cpo_price_imp1_dev_4pya = "CPO price signal deviation from 4 past year average",
                 wa_cpo_price_imp1_dev_4pya_lag1 = "CPO price signal deviation from 4 past year average (lagged)", 
                 # LR CPO variables
                 wa_cpo_price_imp1_2pya = "CPO price signal, 2 past year average",
                 wa_cpo_price_imp1_2pya_lag1 = "CPO price signal, 2 past year average (lagged)",
                 wa_cpo_price_imp1_yoyg_2pya = "CPO price signal y-o-y growth rate, 2 past year average",
                 wa_cpo_price_imp1_yoyg_2pya_lag1 = "CPO price signal y-o-y growth rate, 2 past year average (lagged)",
                 wa_cpo_price_imp1_3pya = "CPO price signal, 3 past year average",
                 wa_cpo_price_imp1_3pya_lag1 = "CPO price signal, 3 past year average (lagged)",
                 wa_cpo_price_imp1_yoyg_3pya = "CPO price signal y-o-y growth rate, 3 past year average",
                 wa_cpo_price_imp1_yoyg_3pya_lag1 = "CPO price signal y-o-y growth rate, 3 past year average (lagged)",
                 wa_cpo_price_imp1_4pya = "CPO price signal, 4 past year average",
                 wa_cpo_price_imp1_4pya_lag1 = "CPO price signal, 4 past year average (lagged)",
                 wa_cpo_price_imp1_yoyg_4pya = "CPO price signal y-o-y growth rate, 4 past year average",
                 wa_cpo_price_imp1_yoyg_4pya_lag1 = "CPO price signal y-o-y growth rate, 4 past year average (lagged)",
                 ## controls
                 lucpfip_pixelcount_total_lag1 = "LUCPFIP (pixels, lagged)",
                 n_reachable_uml = "# reachable UML mills",
                 n_reachable_uml_lag1 = "# reachable UML mills (lagged)",
                 wa_pct_own_cent_gov_imp = "Local government mill ownership (pct.)",
                 wa_pct_own_cent_gov_imp_lag1 = "Local government mill ownership (pct., lagged)",
                 wa_pct_own_loc_gov_imp = "Local government mill ownership (pct.)",
                 wa_pct_own_loc_gov_imp_lag1 = "Local government mill ownership (pct., lagged)",
                 wa_pct_own_nat_priv_imp = "Domestic private mill ownership (pct.)",
                 wa_pct_own_nat_priv_imp_lag1 = "Domestic private mill ownership (pct., lagged)",
                 wa_pct_own_for_imp = "Foreign mill ownership (pct.)",
                 wa_pct_own_for_imp_lag1 = "Foreign mill ownership (pct., lagged)"
                 ))


catchment_radius <- 5e4
island <- c("Sumatra", "Kalimantan", "Papua")
outcome_variable <- "lucpfip_pixelcount_total"
commo <- c("ffb","cpo")
dynamics <- "both"

compare_estimators <- function(catchment_radius, island, outcome_variable, commo, dynamics){
  
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))
  # names(d)
  # 
  # length(unique(d$parcel_id)) # length cross sections: 47648
  # unique(d$year) # years from 2001 to 2015
  
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
  

  ## Specifications
  if(outcome_variable == "lucpfip_pixelcount_total"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucpfip_pixelcount_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucpfip_pixelcount_total ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucpfip_pixelcount_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    if(commo == "both" & dynamics == "SR"){
      specification <- lucpfip_pixelcount_total ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
        
    if(commo == "both" & dynamics == "both"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  if(outcome_variable == "lucpfip_ha_total"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucpfip_ha_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucpfip_ha_total ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucpfip_ha_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "SR"){
      specification <- lucpfip_ha_total ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "both"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  if(outcome_variable == "lucfip_pixelcount_30th"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucfip_pixelcount_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucfip_pixelcount_30th ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucfip_pixelcount_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "SR"){
      specification <- lucfip_pixelcount_30th ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "both"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  if(outcome_variable == "lucfip_ha_30th"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucfip_ha_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucfip_ha_30th ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucfip_ha_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "SR"){
      specification <- lucfip_ha_30th ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "both"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  
  ## Regressions
  poiglm_est <- fixest::feglm(specification, data = d, family = "poisson", notes = FALSE)
  
  qpoiglm_est <- fixest::feglm(specification, data = d, family = "quasipoisson", notes = FALSE)
  
  poimlm_est <- fixest::femlm(specification, data = d, family = "poisson", notes = FALSE)
  
  nb_est <- fixest::fenegbin(specification, data = d, notes = FALSE)
  
  gauglm_est <- fixest::feglm(specification, data = d, family = "gaussian", notes = FALSE)
  
  # title
  if(length(island) == 1){
    table_title <- paste0("Estimator comparisons, ",island, 
                        " catchment radius of ", catchment_radius/1000,"km, ") 
  }else{
    table_title <- paste0("Estimator comparisons, all islands, catchment radius of ", 
                          catchment_radius/1000,"km, ") 
  }
  
  # # file 
  # if(length(island) == 1){
  #   table_file <- file.path(paste0("outputs/regressions/tables/est_comparisons_",
  #                                  island,"_",
  #                                  parcel_size/1000,"km_",
  #                                  catchment_radius/1000,"km_",
  #                                  outcome_variable, 
  #                                  ".tex")) 
  # }else{
  #   table_file <- file.path(paste0("outputs/regressions/tables/est_comparisons_all_",
  #                                  parcel_size/1000,"km_",
  #                                  catchment_radius/1000,"km_",
  #                                  outcome_variable, 
  #                                  ".tex")) 
  # }


  return(etable(poiglm_est, qpoiglm_est, poimlm_est, nb_est, gauglm_est, 
         cluster = ~ parcel_id, 
         tex = TRUE,
         # file = table_file, 
         # replace = TRUE,
         title = table_title,
         subtitles = c("Poisson GLM", "Quasi-Poisson GLM", "Poisson MLM", "Negative Binomial", "Gaussian GLM"),
         family = FALSE,
         coefstat = "confint",
         sdBelow = TRUE,
         dict = TRUE))  
  
 rm(d)
}


OV <- "lucpfip_pixelcount_total"
COMMO <- "ffb"
DYN <- "both"
CR <- 1e4

# "All" islands
ISL <- c("Sumatra", "Kalimantan", "Papua")
for(CR in c(1e4, 3e4, 5e4)){
    for(COMMO in c("ffb", "cpo", "both")){
      for(DYN in c("SR", "LR", "both")){
        for(OV in c("lucpfip_pixelcount_total", "lucpfip_ha_total", 
                    "lucfip_pixelcount_30th", "lucfip_ha_30th")){
        compare_estimators(catchment_radius = CR, 
                     island = ISL, 
                     outcome_variable = OV, 
                     commo = COMMO, 
                     dynamics = DYN)
      }
    }
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### COMPARE FIXED-EFFECTS WITH QUASIPOISSON #### 
catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
outcome_variable <- "lucpfip_pixelcount_total"
dynamics <- FALSE
commo <- c("ffb", "cpo")
yoyg <- FALSE
short_run <- "unt level" # from c("unt_level", "dev", "yoyg")
imp <- 1
x_pya <- 2 # 2, 3 or 4
lag_or_not <- "_lag1" # from c("_lag1", "")
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
oneway_cluster <- ~parcel_id


compare_fe_all_islands <- function(catchment_radius, # c(1e4, 3e4, 5e4)
                                   island, # c("Sumatra", "Kalimantan", "Papua", c("Sumatra", "Kalimantan", "Papua"))
                                   outcome_variable, # c("lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucfip_pixelcount_intact", "lucfip_pixelcount_degraded", "lucfip_pixelcount_total",)
                                   dynamics, # TRUE or FALSE
                                   commo, # c("ffb", "cpo", c("ffb", "cpo"))
                                   yoyg, # TRUE or FALSE 
                                   short_run, # sub or full vector from c("unt level", "yoyg", "dev") 
                                   imp, # c(1,2)
                                   x_pya, # c(2, 3, 4)
                                   lag_or_not, # c("_lag1", "")
                                   controls, # character vectors of base names of controls (don't specify in their names)
                                   weights,
                                   oneway_cluster # formula if clusters are to be specified (see argument cluster in ?fixest::etable)
                                   ){
  
  
  # DATA
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))
  
  # subsampling
  d <- d[d$island %in% island,]
  
  if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
    d <- d[d$any_pfc2000_total,]
  }
  if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
     outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
     outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
    d <- d[d$any_fc2000_30th,]
  }
  
  ### Specifications
  if(dynamics == FALSE){ 
    # Only the overall effect is investigated then, no short vs. long run
    # In this case, the variables have name element _Xya_ with X the number of years over which the mean has been 
    # computed, always including the contemporaneous record. See add_parcel_variables.R
    #short_run <- ""
    
    # if we omit one commodity 
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not)
      }else{
        regressors <- paste0("wa_",commo,"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)
      }
    }
    # if we don't omit a commodity. 
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not))

      }else{
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)) 
      }
    }
  }
  
  # if, on the other hand, we want to disentangle the short run effects from the long run's. 
  if(dynamics == TRUE){ 
    
    if(length(commo) == 1){
      if(yoyg == TRUE){
          regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                          paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
      }
      if(yoyg == FALSE){
        if(short_run == "unt level"){
          regressors <- c(paste0("wa_",commo,"_price_imp",imp,lag_or_not), # SR measure
                          paste0("wa_",commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
        }
        if(short_run == "dev"){
          regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not), # SR measure
                          paste0("wa_",commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
        }  
      }
    }
  
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0("wa_",commo[1],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not))    
      }
      if(yoyg == FALSE){
        if(short_run == "unt level"){
          regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,lag_or_not),# FFB SR measure
                          paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not), # FFB LR measure 
                          paste0("wa_",commo[2],"_price_imp",imp,lag_or_not),# CPO SR measure
                          paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # CPO LR measure 
        }
        if(short_run == "dev"){
          regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                          paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not),
                          paste0("wa_",commo[2],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                          paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not))
        }
      }
    }
  }
  
  # list to be filled with formulae of different fixed-effect models
  fe_model_list <- list() 
  
  # character vector of fixed effect specifications to be compared in a single regression table. 
  # Should be given in a fixest format.
  if(length(island) == 1){
    fixed_effects <- c("parcel_id", 
                       "parcel_id + year", 
                       "parcel_id + province^year", 
                       "parcel_id + district^year")
  }else{
    fixed_effects <- c("parcel_id", 
                       "parcel_id + year", 
                       "parcel_id + island^year", 
                       "parcel_id + province^year", 
                       "parcel_id + district^year") 
  } 
  
  # list model formulae with different fixed effects
  for(fe in fixed_effects){
    fe_model_list[[match(fe, fixed_effects)]] <- as.formula(paste0(outcome_variable,
                                                                " ~ ",
                                                                paste0(regressors, collapse = "+"),
                                                                " + ",
                                                                paste0(paste0(controls,lag_or_not), collapse = "+"),
                                                                " | ",
                                                                fe))
  }
  
# run quasipoisson regression for each fixed-effect model 
  if(weights == TRUE){
    var_weights <- d$sample_coverage_lag1/100 
    fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                          data = d, 
                          family = "quasipoisson", 
                          notes = FALSE, 
                          weights = var_weights)
  }else{
    fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                          data = d, 
                          family = "quasipoisson", 
                          notes = FALSE)
  }

  
  
  ### Return tables
  # title
  if(length(island) == 1){
    table_title <- paste0("Fixed-effect comparisons, ",island, 
                          " catchment radius of ", catchment_radius/1000,"km ") 
  }else{
    table_title <- paste0("Fixed-effect comparisons, all islands, catchment radius of ", 
                          catchment_radius/1000,"km ") 
  }
  
  # esttable(fe_reg_list,            
  #          se = "cluster", 
  #          drop = c("own", "reachable"))
  
  ## ONE WAY CLUSTERING
  if(length(island)>1){
    etable(fe_reg_list, 
           #cluster = oneway_cluster,
           se = "cluster",
           tex = TRUE,
           # file = table_file, 
           # replace = TRUE,
           title = table_title,
           # subtitles = c("FE: grid cell", "FE: year", "FE: grid cell + year", "FE: grid cell + island*year", "FE: grid cell + province*year", "FE: grid cell + district*year"),
           family = TRUE,
           drop = c("own", "reachable"),
           coefstat = "confint",
           sdBelow = FALSE,
           yesNoFixef = "X",
           fitstat = c("sq.cor"),
           dict = TRUE)
  }else{
    etable(fe_reg_list,
           #cluster = oneway_cluster,
           se = "cluster",
           tex = TRUE,
           # file = table_file, 
           # replace = TRUE,
           title = table_title,
           # subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
           # subtitles = c("FE: grid cell", "FE: year", "FE: grid cell + year", "FE: grid cell + island*year", "FE: grid cell + province*year", "FE: grid cell + district*year"),
           family = TRUE,
           drop = c("own", "reachable"),
           coefstat = "confint",
           sdBelow = FALSE,
           yesNoFixef = "X",
           fitstat = c("sq.cor"),
           dict = TRUE)
  }
  
  
  # ## TWO WAY CLUSTERING
  # if(length(island)>1){
  #   etable(qpoiglm_est_lucfip_ufe, 
  #          qpoiglm_est_lucfip_tfe, 
  #          qpoiglm_est_lucfip_twfe, 
  #          qpoiglmspec_lucfip_iyfe, 
  #          qpoiglmspec_lucfip_dyfe, 
  #          cluster = twoway_cluster,
  #          se = "twoway",
  #          tex = TRUE,
  #          # file = table_file, 
  #          # replace = TRUE,
  #          title = table_title,
  #          subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
  #          family = FALSE,
  #          drop = c("own", "reachable"),
  #          coefstat = "confint",
  #          sdBelow = TRUE,
  #          dict = TRUE)
  # }else{
  #   etable(qpoiglm_est_lucfip_ufe, 
  #          qpoiglm_est_lucfip_tfe, 
  #          qpoiglm_est_lucfip_twfe, 
  #          qpoiglmspec_lucfip_dyfe, 
  #          cluster = twoway_cluster,
  #          se = "twoway",
  #          tex = TRUE,
  #          # file = table_file, 
  #          # replace = TRUE,
  #          title = table_title,
  #          subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
  #          family = FALSE,
  #          drop = c("own", "reachable"),
  #          coefstat = "confint",
  #          sdBelow = TRUE,
  #          dict = TRUE)
  # }
  rm(d)
}



# Run on all islands
# "All" islands
ISL <- c("Sumatra", "Kalimantan", "Papua")
OV <- "lucpfip_pixelcount_total"
CR <- 3e4
YOYG <- FALSE

ISL <- "Kalimantan"

DYN <- TRUE
XPYA <- 2

#for(YOYG in c(0, 1)){ # put this first because these are not comaparable measures and hence coeff. 
#for(CR in c(3e4, 5e4)){
  #for(SR in c("unt level", "dev")){
for(ISL in c("Sumatra", "Kalimantan")){
for(DYN in c(0,1)){
  for(XPYA in c(2, 3, 4)){
    compare_fe_all_islands(catchment_radius = CR, 
                           island = ISL,
                           outcome_variable = OV,
                           dynamics = DYN,
                           commo = c("ffb", "cpo"), 
                           yoyg = YOYG,
                           short_run = "unt level", # does not matter if dynamics == FALSE
                           imp = 1,
                           x_pya = XPYA,
                           lag_or_not = "_lag1", # bien vérifier ça !  
                           controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml"),
                           weights = FALSE,
                           oneway_cluster = ~parcel_id)
  }
}
} 
#}
#}





##### DEMAND FOR LUCFP #####
catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
island <- "Sumatra"
outcome_variable <- "lucfip_pixelcount_30th"
outcome_variable <- "lucpfip_pixelcount_total"
x_pya <- 3
lag_or_not <- "_lag1"
fixed_effects <- "parcel_id + district^year"
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
weights <- TRUE

### Prepare the data ### 

  # DATA
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))
  
  # subsampling
  d <- d[d$island %in% island,]
  
  if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
    d <- d[d$any_pfc2000_total,]
  }
  if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
     outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
     outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
    d <- d[d$any_fc2000_30th,]
  }

  ### compute the aggregation factor
  length(unique(d$parcel_id))*length(unique(d$year)) == nrow(d)
  
  years <- unique(d$year)
  n_parcels_within_cr <- c() 
  for(t in years){
  n_parcels_within_cr[[match(t, years)]] <- d[d$year == t & d$n_reachable_uml > 0,] %>% nrow()
  }
  aggr_factor <- sum(n_parcels_within_cr)
  aggr_factor < nrow(d)
  
  ### Run the regression ### 
  
  # there is no dynamics (distinction btw SR and LR effects) here as we aim at 
  # estimating the effect of a tax that remains in time. 
  regressors <- c(paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not),
                  paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)) 
  
  
  fml <- as.formula(paste0(outcome_variable,
                           " ~ ",
                           paste0(regressors, collapse = "+"),
                           " + ",
                           paste0(paste0(controls,lag_or_not), collapse = "+"),
                           " | ",
                           fixed_effects))
  
  if(weights == TRUE){
    var_weights <- d$sample_coverage_lag1/100 
    results <- fixest::feglm(fml,
                             data = d, 
                             family = "quasipoisson", 
                             notes = FALSE, 
                             weights = var_weights)
  }else{
    results <- fixest::feglm(fml,
                             data = d, 
                             family = "quasipoisson", 
                             notes = TRUE)
  }
  
  ## save coefficients 
  coeff_ffb <- results$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
  coeff_cpo <- results$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
  # names(mult_effect_ffb) <- NULL
  # names(mult_effect_cpo) <- NULL
  
  # this is the total number of cells within CR and island, not the the number of cells used for estimation. 
  # this is matter of external validity.
  # therefore we read in the dataframes of grid cells within particular catchment radii, 
  # we adjust for the same restrictions (island and positive baseline forest extent) 
  # and count the number of parcels
  
  # n_cells <- length(unique(d$parcel_id))*length(unique(d$year))
  #names(n_cells) <- "n_cells"
  # then, multiply by 15, the number of years, if we want to compare with observed accumulated amouts
  # or if we want to say about the whole period. 
  # But keep the cross-sectional length (n_cells) if we want to say something like "each year", or "annually". 
  
  
  ## FITTED VALUE AT THE AVERAGE 
  covariate_means <- c(sapply(regressors, FUN = function(var){mean(d[,var], na.rm = TRUE)}),
                       sapply(paste0(controls, lag_or_not), FUN = function(var){mean(d[,var], na.rm = TRUE)}))
  
  linear_predictors_atavg <- c()
  for(i in 1:length(covariate_means)){
    linear_predictors_atavg[i] <- covariate_means[i]*results$coefficients[i]
  }
  # the issue with this is that covariate means are for the whole sample and not just the obs. 
  # used in the regression. (but is it very different, since those not used in the regression
  # are not precisely because they are missing and hence don't count in the mean)
  fitted_value_atavg <- exp(sum(linear_predictors_atavg) + mean(results$sumFE))
  fitted_value_atavg <- fitted_value_atavg*(27.8*27.6)/(1e7)

  ## FITTED VALUE AT THE MEDIAN 
  covariate_medians <- c(sapply(regressors, FUN = function(var){median(d[,var], na.rm = TRUE)}),
                       sapply(paste0(controls, lag_or_not), FUN = function(var){median(d[,var], na.rm = TRUE)}))
  
  linear_predictors_atmed <- c()
  for(i in 1:length(covariate_medians)){
    linear_predictors_atmed[i] <- covariate_medians[i]*results$coefficients[i]
  }

  fitted_value_atmed <- exp(sum(linear_predictors_atmed) + median(results$sumFE))
  fitted_value_atmed <- fitted_value_atmed*(27.8*27.6)/(1e7)
  
  ## AVERAGE FITTED VALUE 
  avg_fitted_value <- mean(results$fitted.values) # this is equal to mean(exp(results$linear.predictors))
  # it is to say that linear.predictors already encompass the fixed effects.
  # It's in pixelcount, convert it to thousand hectares
  avg_fitted_value <- avg_fitted_value*(27.8*27.6)/(1e7)
  
  ## AVERAGE ACTUAL OUTCOME
  ## if we use actual lucpfip and not fitted values
  avg_lucpfip <- mean(d[,outcome_variable])
  avg_lucpfip <- avg_lucpfip*(27.8*27.6)/(1e7)
  
  ## MARGINAL EFFECTS 
  # for a premium / tax of USD10/ton FFB (~1/3 of a std.dev) (and rescale to hectares by *1e3)
  ffb_marginal_effect_avg <- coeff_ffb*avg_fitted_value*1e3*10 # average partial effect
  ffb_marginal_effect_atavg <- coeff_ffb*fitted_value_atavg*1e3*10 # partial effect at average 
  ffb_marginal_effect_atmed <- coeff_ffb*fitted_value_atmed*1e3*10 # partial effect at average 
  ffb_marginal_effect_actual <- coeff_ffb*avg_lucpfip*1e3*10
  cpo_marginal_effect_avg <- coeff_cpo*avg_fitted_value*1e3*10 # average partial effect
  cpo_marginal_effect_atavg <- coeff_cpo*fitted_value_atavg*1e3*10 # partial effect at average 
  cpo_marginal_effect_atmed <- coeff_cpo*fitted_value_atmed*1e3*10 # partial effect at average 
  cpo_marginal_effect_actual <- coeff_cpo*avg_lucpfip*1e3*10
  
  ffb_marginal_effect_avg 
  ffb_marginal_effect_atavg 
  ffb_marginal_effect_atmed 
  ffb_marginal_effect_actual 
  cpo_marginal_effect_avg 
  cpo_marginal_effect_atavg 
  cpo_marginal_effect_atmed 
  cpo_marginal_effect_actual 
  

  ######### test manual computation of partial effects at average with mfx package ###################
  # library(mfx)
  # fml2 <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_4ya_lag1 + wa_cpo_price_imp1_4ya_lag1 + 
  #   wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + 
  #   wa_pct_own_for_imp_lag1 + n_reachable_uml_lag1
  # 
  # mfx <- poissonmfx(formula = fml2, 
  #            data = d,
  #            atmean = TRUE)
  # 
  # results <- fixest::feglm(fml2,
  #                          data = d, 
  #                          family = "poisson", 
  #                          notes = TRUE)
  # 
  # # now reproduce exactly the same, manually: 
  # 
  # 
  # coeff_ffb_mfx <- mfx$fit$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
  # coeff_cpo_mfx <- mfx$fit$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
  # coeff_ffb_fixest <- results$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
  # coeff_cpo_fixest <- results$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
  # 
  # # they are not exactly equal 
  # coeff_ffb_mfx == coeff_ffb_fixest
  # round(coeff_ffb_mfx, 6) == round(coeff_ffb_fixest,6)
  # 
  # ## at average, manually
  # covariate_means <- c(sapply(regressors, FUN = function(var){mean(d[,var], na.rm = TRUE)}),
  #                      sapply(paste0(controls, lag_or_not), FUN = function(var){mean(d[,var], na.rm = TRUE)}))
  # 
  # # based on mfx (glm)
  # linear_predictors_mfx <- c()
  # for(i in 1:length(covariate_means)){
  #   linear_predictors_mfx[i] <- covariate_means[i]*mfx$fit$coefficients[i+1] # +1 bc of the intercept
  # }
  # fitted_value_atavg_mfx <- exp(sum(linear_predictors_mfx) + mfx$fit$coefficients[1]) # linear predictors don't integrate the intercept.
  # ffb_marginal_effect_atavg_mfx <- coeff_ffb_mfx*fitted_value_atavg_mfx
  # cpo_marginal_effect_atavg_mfx <- coeff_cpo_mfx*fitted_value_atavg_mfx
  # 
  # # based on fixest coeffs
  # linear_predictors_fixest <- c()
  # for(i in 1:length(covariate_means)){
  #   linear_predictors_fixest[i] <- covariate_means[i]*results$coefficients[i+1]
  # }
  # 
  # fitted_value_atavg_fixest <- exp(sum(linear_predictors_fixest) + results$coefficients[1])
  # ffb_marginal_effect_atavg_fixest <- coeff_ffb_fixest*fitted_value_atavg_fixest
  # cpo_marginal_effect_atavg_fixest <- coeff_cpo_fixest*fitted_value_atavg_fixest
  # 
  # # because coeffs are not exactly equal, linear predictors are not either, nor fitted values and marginal effects
  # linear_predictors_mfx == linear_predictors_fixest
  # round(linear_predictors_mfx, 6) == round(linear_predictors_fixest, 6)
  # fitted_value_atavg_mfx == fitted_value_atavg_fixest
  # ffb_marginal_effect_atavg_mfx == ffb_marginal_effect_atavg_fixest
  # round(ffb_marginal_effect_atavg_mfx, 6) == round(ffb_marginal_effect_atavg_fixest, 6)
  # 
  # # but this does not explain the gap with the marginal effect directly computed by mfx
  # mfx$mfxest[1,1] == ffb_marginal_effect_atavg_mfx
  # round(mfx$mfxest[1,1], 6) == round(ffb_marginal_effect_atavg_mfx,6) # is FALSE
  # 
  # mfx <- poissonmfx(formula = fml2, 
  #            data = d,
  #            atmean = FALSE)
  # 
  # # recover the fitted value at the average from mfx marginal effect
  # # bc we know the marginal effect at average is the product of the coeff and the fitted value at the average
  # mfx$mfxest[1,1]/coeff_ffb_mfx == fitted_value_atavg_mfx
  # # is false, therefore it's really the fitted value at the average that is not computed exactly the same in 
  # # mfx and manually 
  # 
  # ## AVERAGE FITTED VALUE 
  # avg_fitted_value <- mean(results$fitted.values)
  # ffb_marginal_effect_avg <- coeff_ffb*avg_fitted_value # average partial effect
  # cpo_marginal_effect_avg <- coeff_cpo*avg_fitted_value # average partial effect
  # 
  # poissonmfx(formula = fml2, 
  #                   data = d,
  #                   atmean = FALSE)
  # 
  # # so for average partial effect, results are equal. 
  # # NOT for partial effect at the average... NOW WITH THE INTERCEPT ADDED IT'S VERY SIMILAR, BUT STILL NOT EQUAL !!!
  # # --> try to compute manually on a balanced panel (remove records that have at least one rhs missing)
  # # or build a basic balanced df. 
  # # it's not coeff that differ since for average partial effect it works. 
  # # voir ce qu'est glm$effects (mfx$fit$effects)
  # 
  # 
  # ### and now with {margins}
  # library(margins)
  #### end of test ####
  
  
  ## FFB DEMAND FUNCTIONS
  ffb_demand_avg <- function(tax){
    avg_fitted_value*aggr_factor/exp(coeff_ffb*tax)
  }
  ffb_demand_atavg <- function(tax){
    fitted_value_atavg*aggr_factor/exp(coeff_ffb*tax)
  }
  ffb_demand_atmed <- function(tax){
    fitted_value_atmed*aggr_factor/exp(coeff_ffb*tax)
  }
  ffb_demand_actual <- function(tax){
    avg_lucpfip*aggr_factor/exp(coeff_ffb*tax)
  }
  
  ffb_demand_avg(0)
  ffb_demand_atavg(0)
  ffb_demand_atmed(0) 
  ffb_demand_actual(0) 
  
  avg_lucpfip*nrow(d) # we get exactly the same amount as the actual accumulated LUCPFIP in table 2. 
    
 
  
  # for plotting purpose: 
  ffb_inv_demand_avg <- function(CF){
    (1/coeff_ffb)*(log(avg_fitted_value*aggr_factor) - log(CF))
  }
  
  ffb_inv_demand_atavg <- function(CF){
    (1/coeff_ffb)*(log(fitted_value_atavg*aggr_factor) - log(CF))
  }
  
  ffb_inv_demand_atmed <- function(CF){
    (1/coeff_ffb)*(log(fitted_value_atmed*aggr_factor) - log(CF))
  }
  
  ffb_inv_demand_actual <- function(CF){
    (1/coeff_ffb)*(log(avg_lucpfip*aggr_factor) - log(CF))
  }
  
  ggplot(data = data.frame(CF=c(0, ffb_demand_avg(0))), 
          aes(x=CF)) + 
          stat_function(fun=ffb_inv_demand_avg)

  ggplot(data = data.frame(CF=c(0, ffb_demand_atavg(0))), 
         aes(x=CF)) + 
    stat_function(fun=ffb_inv_demand_atavg)
  
  ggplot(data = data.frame(CF=c(0, ffb_demand_atmed(0))), 
         aes(x=CF)) + 
    stat_function(fun=ffb_inv_demand_atmed)
  
  ggplot(data = data.frame(CF=c(0, ffb_demand_actual(0))), 
         aes(x=CF)) + 
    stat_function(fun=ffb_inv_demand_actual) 
  
## CPO DEMAND FUNCTIONS
  cpo_demand_avg <- function(tax){
    avg_fitted_value*aggr_factor/exp(coeff_cpo*tax)
  }
  cpo_demand_atavg <- function(tax){
    fitted_value_atavg*aggr_factor/exp(coeff_cpo*tax)
  }
  cpo_demand_atmed <- function(tax){
    fitted_value_atmed*aggr_factor/exp(coeff_cpo*tax)
  }
  cpo_demand_actual <- function(tax){
    avg_lucpfip*aggr_factor/exp(coeff_cpo*tax)
  }
  
  cpo_demand_avg(0)
  cpo_demand_atavg(0)
  cpo_demand_atmed(0)
  cpo_demand_actual(0) # we get exactly the same amount as the actual accumulated LUCPFIP in table 2. 
  
  
  
  # for plotting purpose: 
  cpo_inv_demand_avg <- function(CF){
    (1/coeff_cpo)*(log(avg_fitted_value*aggr_factor) - log(CF))
  }
  
  cpo_inv_demand_atavg <- function(CF){
    (1/coeff_cpo)*(log(fitted_value_atavg*aggr_factor) - log(CF))
  }
  
  cpo_inv_demand_atmed <- function(CF){
    (1/coeff_cpo)*(log(fitted_value_atmed*aggr_factor) - log(CF))
  }
  
  cpo_inv_demand_actual <- function(CF){
    (1/coeff_cpo)*(log(avg_lucpfip*aggr_factor) - log(CF))
  }
  
  ggplot(data = data.frame(CF=c(0, cpo_demand_avg(0))), 
         aes(x=CF)) + 
    stat_function(fun=cpo_inv_demand_avg) 
    
    ggplot(data = data.frame(CF=c(0, cpo_demand_atavg(0))), 
           aes(x=CF)) + 
    stat_function(fun=cpo_inv_demand_atavg)
  
  ggplot(data = data.frame(CF=c(0, cpo_demand_atmed(0))), 
         aes(x=CF)) + 
    stat_function(fun=cpo_inv_demand_atmed)
  
  ggplot(data = data.frame(CF=c(0, cpo_demand_actual(0))), 
         aes(x=CF)) + 
    stat_function(fun=cpo_inv_demand_actual) 
  

# effect of a premium of 25$/ton FFB
ffb_demand_atavg(0) - ffb_demand_atavg(25)
# effect of a premium of 100$/ton CPO
cpo_demand_atavg(0) - cpo_demand_atavg(100)



 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## test whether the FE are in linear.predictors 
# linear.predictors, sumFE and fitted.values all have the same length: the number of obs. used in the reg
length(results$linear.predictors)==length(results$sumFE)
length(results$sumFE)==length(results$fitted.values)

# the fitted.values are equal to the exponential of the linear.predictors (here)
results$linear.predictors[1] 
results$sumFE[1]
results$fitted.values[4] == exp(results$linear.predictors[4])
results$linear.predictors[1] + results$sumFE[1]

fml_nofe <- as.formula(paste0(outcome_variable,
                              " ~ ",
                              paste0(regressors, collapse = "+"),
                              " + ",
                              paste0(paste0(controls,lag_or_not), collapse = "+")))

results_nofe <- fixest::feglm(fml_nofe,
                              data = d, 
                              family = "quasipoisson", 
                              notes = TRUE)
length(results_nofe$linear.predictors) - length(results$linear.predictors)

any(results$linear.predictors %in% results_nofe$linear.predictors)

any(results_nofe$linear.predictors %in% results$linear.predictors)

results_nofe$fitted.values[4] == exp(results_nofe$linear.predictors[4])

# so linear predictors are different when there are FE and where there is no. 
# which means that linear predictors INCLUDE FE. 

if(commo == "ffb" & dynamics == "both-yoyg" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- as.symbol(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "ffb" & dynamics == "both-3pya" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-yoyg" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1  +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-3pya" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year 
}
if(commo == "both" & dynamics == "both-yoyg" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
      wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "both" & dynamics == "both-3pya" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
      wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "ffb" & dynamics == "both-yoyg" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "ffb" & dynamics == "both-3pya" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-yoyg" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya  +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-3pya" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year 
}
if(commo == "both" & dynamics == "both-yoyg" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
      wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "both" & dynamics == "both-3pya" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
      wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}


### Regressions
qpoiglm_est_lucfip_ufe <- fixest::feglm(spec_lucfip_ufe, data = d, family = "quasipoisson", notes = FALSE)
qpoiglm_est_lucfip_tfe <- fixest::feglm(spec_lucfip_tfe, data = d, family = "quasipoisson", notes = FALSE)
qpoiglm_est_lucfip_twfe <- fixest::feglm(spec_lucfip_twfe, data = d, family = "quasipoisson", notes = FALSE)
if(length(island) > 1){
  qpoiglmspec_lucfip_iyfe <- fixest::feglm(spec_lucfip_iyfe, data = d, family = "quasipoisson", notes = FALSE)
}
qpoiglmspec_lucfip_dyfe <- fixest::feglm(spec_lucfip_dyfe, data = d, family = "quasipoisson", notes = FALSE)
qpoiglmspec_lucfip_pyfe <- fixest::feglm(spec_lucfip_pyfe, data = d, family = "quasipoisson", notes = FALSE)

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







# PREMIER SET 126 -31
# 2EME SET 41.85 - 12
# 3EME SET 219.15 - 66
36+123.2+39.6+15.4+4.95




