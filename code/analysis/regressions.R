
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
                   "modelsummary")
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
outcome_variable <- "lucfip_pixelcount_30th"
commo <- "both"
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
# commo <- c("ffb", "cpo")
# lag_or_not <- "_lag1" # from c("_lag1", "")
# dynamics <- FALSE
# short_run <- "unt level" # from c("unt_level", "dev", "yoyg")
# island <- c("Sumatra", "Kalimanta", "Papua")
# catchment_radius <- 3e4
# outcome_variable <- "lucpfip_pixelcount_total"
# fixed_effects <- c("parcel_id", "year", "parcel_id + year", "parcel_id + island^year", "parcel_id + province^year", "parcel_id + district^year")
# x_pya <- 4 # 2, 3 or 4
# imp <- 1
# controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
# oneway_cluster <- ~parcel_id


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
                                   fixed_effects, # character vector of fixed effect specifications to be compared in a single regression table. Should be given in a fixest format. e.g.:  c("parcel_id", "year", "parcel_id + year", "parcel_id + island^year", "parcel_id + province^year", "parcel_id + district^year")
                                   oneway_cluster # formula if clusters are to be specified (see argument cluster in ?fixest::etable)
                                   ){
  
  
  ### 30% canopy density forest 
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
    short_run <- ""
    
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
  fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                        data = d, 
                        family = "quasipoisson", 
                        notes= FALSE)
  
  
  ### Return tables
  # title
  if(length(island) == 1){
    table_title <- paste0("Fixed-effect comparisons, ",island, 
                          " catchment radius of ", catchment_radius/1000,"km, ") 
  }else{
    table_title <- paste0("Fixed-effect comparisons, all islands, catchment radius of ", 
                          catchment_radius/1000,"km, ") 
  }
  
  ## ONE WAY CLUSTERING
  if(length(island)>1){
    etable(fe_reg_list, 
           #cluster = oneway_cluster,
           se = "cluster",
           tex = TRUE,
           # file = table_file, 
           # replace = TRUE,
           title = table_title,
           # subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
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
ISL <- c("Sumatra", "Kalimanta", "Papua")
OV <- "lucpfip_pixelcount_total"
CR <- 5e4
for(YOYG in c(0, 1)){ # put this first because these are not comaparable measures and hence coeff. 
  for(CR in c(3e4, 5e4)){
   #for(SR in c("unt level", "dev")){
     for(XPYA in c(2, 3, 4)){
      compare_fe_all_islands(catchment_radius = CR, 
                         island = ISL,
                         outcome_variable = OV,
                         dynamics = FALSE,
                         commo = c("ffb"), 
                         yoyg = YOYG,
                         short_run = SR, # does not matter if dynamics == FALSE
                         imp = 1,
                         x_pya = XPYA,
                         lag_or_not = "", # bien vérifier ça !  
                         controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml"),
                         fixed_effects = c("parcel_id", "year", "parcel_id + year", "parcel_id + island^year", "parcel_id + province^year", "parcel_id + district^year"),
                         oneway_cluster = ~parcel_id)
     }
   }
 #}
}

# Run per island
for(ISL in c("Sumatra", "Kalimanta", "Papua")){
  for(CR in c(1e4, 3e4, 5e4)){
    compare_fe_all_islands(catchment_radius = CR, 
                           island = ISL)
  }
}












### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


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













