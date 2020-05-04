IF USER HAS Stata/MP: 
The complete data work can be reproduced running master.do with Stata/MP > 14.0
Notes on outside-Stata programs: 
	- GEE: Google Earth Engine scripts are not called from master.do ; they can be provided upon request (they extract descriptive statistics). 
	- R: master.do calls R scripts; it transfers them its current working directory, and each R script invokes paths relative to this. 


IF USER DOES NOT HAVE Stata/MP and wants to run only R scripts. 
The R_project_for_individual_runs.Rproj file should be opened first 
(in particular, for each R script to run in the correct working directory, and with the correct packages and packages' versions). 
Then R scripts can be sourced within this project environment;
There is no master R script because it would make no sense to run all R scripts without running inserted Stata scripts.
How Stata and R scripts are interrelated is described below, as a paste of master.do. 


DATA: 
The original input data are stored in the read only directory ~/LUCFP/data_processing/input_data
This is a several Gb folder, that is not available from Github, but is made available upon request. 
Some data are also downloaded within scripts (e.g. GFC tiles, from prepare_lucpfip.R) and internet connexion is therefore required for those scripts and for the complete replication run. 

The temporary data prepared within this project are stored in the write and read directory  ~/LUCFP/data_processing/temp_data
This folder is also large and hence is not available from Github either. 


RUN TIME 
Several scripts can take a long time to be executed, this is the case in particular for: 
prepare_lucpfip.R
wa_at_parcels.R
Indicative run times should be given in the master.do where relevant. 


PACKAGES
Stata: 


R: 
This project's data processing working directory is home for a renv.lock file and a renv folder. 
The first thing any R script from this project does is to install the renv package in user's default library if it cannot be loaded, 
and install and load the packages used in this project, listed in the renv.lock file and hosted in ~/LUCFP/data_processing/renv/library. 
