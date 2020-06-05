### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
#             EXPLICATIVE VARIABLE DESCRIPTIVE STATISTICS 
# 
#   input   - panel data frame of parcels with the information on number of reachable uml mills
#           as outputted from make_n_reachable_uml.R
#           ---> temp_data/processed_parcels/wa_panel_parcels_reachable_uml_PS_CR.rds
# 
#           - panel data frame of IBS 
#           ---> temp_data/IBS_UML_panel_final.dta
# 
# 
# 
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
# PB: on n'arrive pas à installed leaflet, kableExtra (qui n'est peut être pas nécessaire) et velox (qui n'est pas nécessaire)
neededPackages = c("dplyr", "data.table",  "stringr", "sjmisc", 
                   "foreign", "readxl", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "parallel", "foreach","doParallel")
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
troublePackages <- c("plyr", "kableExtra", "ExPanDaR", 
                     "tidyr", "dplyr", "ggplot2", "corrplot", "lfe", "multiwayvcov", "lmtest", "stargazer", "scales", "shiny", "DT", "openssl", "tictoc", "shinycssloaders", "kableExtra", "rio", "zip", "rlang")
# Attempt to load packages from user's default libraries.
lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 

### NEW FOLDERS USED IN THIS SCRIPT 

### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"


##### IBS MILL PANEL STAT DES ##### 
# all IBS oil palm related establishments. 
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))

# make indicator of our sample of interest (geo-referenced mills) 
ibs$geo_sample <- !is.na(ibs$lat)

# make indicator of being a mill 
is_mill <- ddply(ibs, "firm_id", summarise, 
                 is_mill = sum(out_ton_cpo, na.rm = TRUE) > 0 | 
                           sum(out_val_cpo, na.rm = TRUE) > 0 | 
                           sum(out_ton_pko, na.rm = TRUE) > 0 | 
                           sum(out_val_pko, na.rm = TRUE) > 0 | 
                           sum(in_ton_ffb, na.rm = TRUE) > 0 | 
                           sum(in_val_ffb, na.rm = TRUE) > 0 )
ibs <- merge(ibs, is_mill, by = "firm_id")
rm(is_mill)
length(unique(ibs$firm_id))
length(unique(ibs$firm_id[ibs$is_mill])) 
length(unique(ibs$firm_id[ibs$geo_sample & ibs$is_mill]))
# ibs[ibs$geo_sample & !ibs$is_mill,1:40]
# out of 1473 establishments initially extracted from IBS,
# 1004 are involved at least one year in some FFB processing or CPO or PKO producing. 
# out of which 468 have been geolocalized.
# 2 additional IBS firms have been matched with UML and geolocalized but have no sign of FFB, PKO or CPO activity. 
# (firm_id 2036 and 55630)


# Now we also remove those that have some FFB-CPO-PKO activity but have been identified as refineries. 
# But why would we remove refineries? They may be different but if they input FFB, they have an impact on proximate LUC. 
# Because we excluded refineries from geo_sample, and here the purpose is to compare this sample to the population of mills. 
# let's first see the comparative stat des without removing refineries. 


#### DESCRIPTIVE TABLES #### 

# We want to produce, for a set of ibs variables, the mean, median, std.dev. min and max statistics, 
# for the geo_sample that we are going to use, and the is_mill sample. 

variables <- c("min_year", 
               "ffb_price_imp1", "in_ton_ffb_imp1",
               "cpo_price_imp1", "out_ton_cpo_imp1",
               "pko_price_imp1", "out_ton_pko_imp1",
               "prex_cpo_imp1", 
               "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp")

statistics <- c("mean", "median", "min", "max", "std.dev")

make_des_table_island <- function(ibs_isl){
  
  ## Matrix for geo_sample
  rhs_des_sample <- matrix(NA, nrow = length(variables), ncol = length(statistics))
  row.names(rhs_des_sample) <- variables
  colnames(rhs_des_sample) <- c(statistics)
  
  for(var in variables){
    rhs_des_sample[var,statistics] <- summarise(dplyr::filter(ibs_isl, geo_sample == TRUE),
                                        mean = mean(get(var), na.rm=T),
                                        median = median(get(var), na.rm= TRUE), 
                                        min = min(get(var), na.rm= TRUE),
                                        max = max(get(var), na.rm= TRUE),
                                        std.dev. = sd(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      formatC(digits =4, format = "fg")
  }

  # number of establishments in sub-group
  # N_sample_row <- c(length(unique(ibs_isl[ibs_isl$geo_sample==TRUE, "firm_id"])), 
  #                           rep("", length(statistics))) 
  # rhs_des_sample <- rbind(N_sample_row, rhs_des_sample)
  
  ## Matrix for all ibs mills
  rhs_des_pop <- matrix(NA, nrow = length(variables), ncol = length(statistics))
  row.names(rhs_des_pop) <- variables
  colnames(rhs_des_pop) <- c(statistics)
  
  for(var in variables){
    rhs_des_pop[var,statistics] <- summarise(dplyr::filter(ibs_isl, is_mill == TRUE),
                                         mean = mean(get(var), na.rm=T),
                                         median = median(get(var), na.rm= TRUE), 
                                         min = min(get(var), na.rm= TRUE),
                                         max = max(get(var), na.rm= TRUE),
                                         std.dev. = sd(get(var), na.rm= TRUE)) %>% 
      as.matrix()%>% 
      formatC(digits =4, format = "fg")
  }
  
  # # number of establishments in sub-group
  # N_pop_row <- c(length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"])), 
  #                   rep("", length(statistics))) 
  # rhs_des_pop <- rbind(N_pop_row, rhs_des_pop)
         
  # bind two groups together
  rhs_des <- cbind(rhs_des_sample, rhs_des_pop)    
  
  # add the t-test column
  t_test <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(t_test) <- variables
  colnames(t_test) <- "t test"
  
  for(var in variables){
   test <- t.test(x = ibs_isl[ibs_isl$geo_sample==TRUE,var],
         y = ibs_isl[ibs_isl$is_mill==TRUE,var], 
         alternative = "two.sided", 
         mu = 0, 
         paired = F, 
         var.equal = FALSE)
   
   t_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
# interpretation: if the 95% CI includes 0, then the difference in means between two samples 
# is not statistically different from 0. Hence, the two samples are "similar" in means 
# with respect to the variable tested. 
# (In other words, we cannot reject the null hypothesis that the difference is null
# -i.e. the two samples are alike - with 95% confidence) 
  
  # add the Smirnov test
  ks_test <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(ks_test) <- variables
  colnames(ks_test) <- "KS test"
  
  for(var in variables){
    test <- ks.test(x = ibs_isl[ibs_isl$geo_sample==TRUE,var],
                   y = ibs_isl[ibs_isl$is_mill==TRUE,var], 
                   alternative = "two.sided", 
                   exact = FALSE)
    
    ks_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  
  # intepretation: if p-value < 0.05 we cannot reject with 95% confidence that 
  # two distributions are equal, implying that they are different. 
  
  rhs_des <- cbind(rhs_des, t_test, ks_test)
  
  
  
  # row names
  row.names(rhs_des) <- c("First year in IBS", "FFB muv (USD/ton)", "FFB input (ton)", 
                          "CPO muv (USD/ton", "CPO output (ton)", 
                          "PKO muv (USD/ton)", "PKO output (ton)", 
                          "CPO export share (%)", 
                          "Central government ownership (%)", 
                          "Local government ownership (%)", 
                          "National private ownership (%)", 
                          "Foreign ownership (%)")
  
  return(rhs_des)
}


#### Print the LateX table code SUMATRA #### 

ibs_isl <- ibs[ibs$island_name=="Sumatra",]
des_table <- make_des_table_island(ibs_isl)

# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$geo_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, Sumatra") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 367 mills" = 5,
                     "All IBS palm oil mills \n n = 677 mills" = 5, 
                     "t-test" = 1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2:11),
              width = "4em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

#### Print the LateX table code KALIMANTAN #### 

ibs_isl <- ibs[ibs$island_name=="Kalimantan",]
des_table <- make_des_table_island(ibs_isl)
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$geo_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, Kalimantan") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "p-value"=1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 83 mills" = 5,
                     "All IBS palm oil mills \n n = 200 mills" = 5, 
                     "t-test"=1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2:11),
              width = "4em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 


#### Print the LateX table code OTHERS #### 
ibs_isl <- ibs[!(ibs$island_name %in% c("Sumatra", "Kalimantan")),]
des_table <- make_des_table_island(ibs_isl)
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$geo_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, other islands") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 51 mills" = 5,
                     "All IBS palm oil mills \n n = 169 mills" = 5, 
                     "t-test" = 1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2:11),
              width = "4em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

#### Print the LateX table code ALL #### 
ibs_isl <- ibs
des_table <- make_des_table_island(ibs_isl)
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$geo_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, all islands") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "mean" = 1, "median" = 1, 
                     "min" = 1, "max" = 1, "std.dev." = 1, 
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 470 mills" = 5,
                     "All IBS palm oil mills \n n = 1004 mills" = 5, 
                     "t-test" = 1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2:11),
              width = "4em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 




##### PARCEL STAT DES ##### 

parcels30 <- readRDS(file.path("temp_data/panel_parcels_ip_final_3km_30CR.rds"))


make_des_table_parcels_island <- function(parcels_isl){
  
  
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### FIGURES ### ### ### ### ### 

any(duplicated(ibs[,c("firm_id", "year")]))

# pre-select variables before exploring the dataset with ExPanD 
# (this is the same set as the one selected in wa_at_parcels.R + geo_sample and is_mill + geographic variables)
ibs <- ibs[, c("firm_id", "year", "is_mill", "geo_sample", 
               "trase_code", "uml_id", "mill_name", "parent_co", "lat", "lon",
               "island_name", "district_name", "kec_name", "village_name", 
               "min_year","est_year", "est_year_imp", "max_year", 
               "ffb_price_imp1", "ffb_price_imp2", "in_ton_ffb_imp1", "in_ton_ffb_imp2", "in_val_ffb_imp1", "in_val_ffb_imp2",
               "cpo_price_imp1","cpo_price_imp2", "out_ton_cpo_imp1", "out_ton_cpo_imp2", "out_val_cpo_imp1", "out_val_cpo_imp2",
               "prex_cpo_imp1", "prex_cpo_imp2",
               "pko_price_imp1","pko_price_imp2", "out_ton_pko_imp1", "out_ton_pko_imp2", "out_val_pko_imp1", "out_val_pko_imp2",
               "prex_pko_imp1", "prex_pko_imp2",
               "export_pct_imp", "revenue_total", "workers_total_imp3",
               "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp", 
               "iv2_imp1", "iv2_imp2", "iv3_imp1", "iv3_imp2", "iv4_imp1", "iv4_imp2", 
               "concentration_10", "concentration_30", "concentration_50")]
ExPanD(df = ibs, cs_id = "firm_id", ts_id = "year")



# PARCEL PANEL - LUCPFIP
parcels10 <- readRDS(file.path("temp_data/panel_parcels_final_lucpfip_3km_10CR.rds"))
parcels30 <- readRDS(file.path("temp_data/panel_parcels_final_lucpfip_3km_30CR.rds"))
parcels50 <- readRDS(file.path("temp_data/panel_parcels_final_lucpfip_3km_50CR.rds"))

parcels <- list(parcels10, parcels30, parcels50)
# convert pixel counts to hectares, rearrange variables and remove pixelcount.
pixel_area <- (27.8*27.6)/(1e4)
# intact
parcels <- lapply(parcels, function(df){mutate(df, lucpfip_ha_intact = lucpfip_pixelcount_intact*pixel_area)}) 
parcels <- lapply(parcels, function(df){dplyr::select(df,parcel_id, year, lucpfip_ha_intact, everything(), -lucpfip_pixelcount_intact)})
# degraded
parcels <- lapply(parcels, function(df){mutate(df, lucpfip_ha_degraded = lucpfip_pixelcount_degraded*pixel_area)}) 
parcels <- lapply(parcels, function(df){dplyr::select(df,parcel_id, year, lucpfip_ha_degraded, everything(), -lucpfip_pixelcount_degraded)})
# total
parcels <- lapply(parcels, function(df){mutate(df, lucpfip_ha_total = lucpfip_pixelcount_total*pixel_area)}) 
parcels <- lapply(parcels, function(df){dplyr::select(df,parcel_id, year, lucpfip_ha_total, everything(), -lucpfip_pixelcount_total)})

summary(parcels10$lucpfip_pixelcount_total)
summary(parcels[[1]]$lucpfip_ha_total)

# there are NAs in outcome variables because they are not observed in 1998-2000. 

# category of parcels that experienced at least one lucpfip pixel event 
# (and hence is relevant to our analysis, it is not a city or water)
for(i in c(1:3)){

  parcels[[i]] <- merge(parcels[[i]],
                        ddply(parcels[[i]], "parcel_id", summarize, 
                             relevant = sum(lucpfip_ha_total, na.rm = TRUE) >0 ), 
                        by = "parcel_id")
}

ExPanD(df = parcels[[2]][parcels[[2]]$relevant,], cs_id = "parcel_id", ts_id = "year")





