*****************************************************************************************************************************************
* Master do file for the LUCFP project. R scripts are called from here, but this requires user to specify the local path to R executable
*****************************************************************************************************************************************

/////!!!!!!\\\\\/////!!!!!!\\\\\/////!!!!!!\\\\\
*** SPECIFY YOUR LOCAL PATH TO YOUR R PROGRAM 
global Rterm_path "C:\Program Files\R\R-3.6.1\bin\R.exe"
/////!!!!!!\\\\\/////!!!!!!\\\\\/////!!!!!!\\\\\

******* NOTHING IS USER SPECIFIC BELOW THIS LINE ********

*** THE PROJECT WORKING DIRECTORY DOES NOT NEED TO BE MANUALLY SPECIFIED, AS IT IS HOME OF THE PRESENT master.do, i.e. ~/LUCFP/data_processing/    
/* 
However, note that if you run the present .do file from a code editor like Sublime Text 3, 
or if you set a working directory prior to running it (in your profile.do for instance), you need to set the current working directory to 
your local absolute path to the present master.do 
*/
cd "$LUCFP"

* This is where all the processed data are stored through the following script execution. 
* Subsequent directories are generated within scripts.  
capture mkdir "temp_data"

*** PACKAGES 
ssc install rsource, replace
* set permanent R options to be applied at each call of an R script from Stata. 
* vanilla is necessary for the R script call to work, but it prevents R from reading its .Rprofile file at start up 
* see https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html 
* and https://stackoverflow.com/questions/12540138/how-can-i-make-my-r-session-vanilla
global Rterm_options "--vanilla"

* Note on rsource: the working directory returned by getwd() at the start of R scripts called from rsource 
* is the current working directory in Stata. 


/* 
*INRA 
global base_path_wd "C:/Users/GUYE/Desktop/opalval"
*MCC
global base_path_wd "C:/Users/guyv/ownCloud/opalval" */


**** explicative variables 

	*** extract palm oil mills from IBS main input/output datasets (IBS_IO) - and clean measurement unit and reshape. 

	do code/explicative_variables/IBS_inputs_preparation.do
	* input: input_data/IBS_IO/IBS_inputs    downloaded here from C:/Users/guyv/ownCloud/opal (2)/download/output/IBS_IO
	* output: temp_data/processed_IBS/prepared_IO/IBS_inputs_prep.dta

	do code/explicative_variables/IBS_outputs_preparation.do
	* input: input_data/IBS_IO/IBS_outputs   downloaded here from C:/Users/guyv/ownCloud/opal (2)/download/output/IBS_IO
	* output: temp_data/processed_IBS/prepared_IO/IBS_outputs_prep.dta

	** merge them together and with main IBS. 
	do code/explicative_variables/merge_IBSmain_IBSIO.do
	* input: temp_data/processed_IBS/prepared_IO/IBS_outputs_prep.dta
	*		 temp_data/processed_IBS/prepared_IO/IBS_inputs_prep.dta
	*		 input_data/si_panel.dta  		  (Rothenberg data) uploaded here from C:/Users/guyv/ownCloud/opal (2)/download/input/IBS_full/si_panel.dta
	* 		 input_data/IBS_final_panel.dta  (Sebastian cleaned data) uploaded here from C:/Users/guyv/ownCloud/opal (2)/build/output/IBS_final_panel.dta

	* output: temp_data/processed_IBS/IBS_PO_98_15.dta

	*** clean IBS mill dataset. 
	** some preparatory work with district and village crosswalks
	do code/explicative_variables/prepare_crosswalks.do
	* input 	input_data/indonesia_spatial/District-Proliferation-Crosswalk_complete.csv
	*           input_data/indonesia_spatial/desa_crosswalk_1998_2014.csv 

	* output:   temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
	*		    temp_data/processed_indonesia_spatial/desa_code_names_98_2014.dta
	
	** clean IBS_PO_98_15 (including geographic variables)
	do code/explicative_variables/cleaning_IBS.do
	* input: temp_data/processed_IBS/IBS_PO_98_15.dta
	*		 input_data/IBS_kraus/IBS_final_panel.dta
	*        input_data/IBS_kraus/IBS_base_prelu.dta
	*        temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
	*        temp_data/processed_indonesia_spatial/desa_code_names_98_2014.dta

	* output: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta

**** mill geolocalization 

	*** automatic matching with Heilmayr mill list. 

		** keep mill obs. of most recent year with a valid desa_id
		do code/mill_geolocalization/keep_valid_recent_desa.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		* output: temp_data/processed_mill_geolocalization/IBSmills_valid_desa.dta

		** give IBS mills a desa geometry
		rsource using "code/mill_geolocalization/make_IBSmills_desageom.R"
		* input: temp_data/processed_mill_geolocalization/IBSmills_valid_desa.dta
		* 		 input_data/indonesia_spatial/village_shapefiles/desa_1998_crosswalked.shp 
		*		 ... 
		*        input_data/indonesia_spatial/village_shapefiles/desa_2009_crosswalked.shp
		* 		 input_data/indonesia_spatial/village_shapefiles/desa_map_2010/indo_by_desa_2010.shp
		
		* output: temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata

		** perform spatial matching and output distinct subsets that will need different processing
		rsource using "code/mill_geolocalization/IBS_UML_matching.R"
		* input: temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata
		*  		 input_data/uml/traseMills_capEstyear.xlsx
		
		* output: temp_data/processed_mill_geolocalization/pre2011_bad_desa_id.dta
		* 		  temp_data/processed_mill_geolocalization/unreferenced_mill_desa.shp (driver geojson)
		*		  temp_data/processed_mill_geolocalization/ibs_unref.dta
		* 		  temp_data/processed_mill_geolocalization/oto.dta
		* 		  temp_data/processed_mill_geolocalization/noto.dta


	*** manual matching / geo-localization

		** Take sub-dataset of mills with observations only since 2011 included and thus surely no desa_id.
		do code/mill_geolocalization/georeferencing_IBSmills_post2010.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta	 
		
		* output: temp_data/processed_mill_geolocalization/mills_to_georeference_post2010.xls

		** Take sub-datasets of mills with at least one observation before 2011 but that cannot match automatically (no valid desa_id, or conflicting matches - "noto"). 
		do code/mill_geolocalization/georeferencing_IBSmills_pre2011.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		* 		 temp_data/processed_mill_geolocalization/pre2011_bad_desa_id.dta
		* 		 temp_data/processed_mill_geolocalization/noto.dta

		* output: temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011.xls
		* 		  temp_data/processed_mill_geolocalization/noto.xls  

		** MANUAL WORK; **NO CODE** **NOTE THAT OUTPUT IS IN READ-ONLY INPUT DATA**
		* input: temp_data/processed_mill_geolocalization/mills_to_georeference.xls
		* 		 temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011.xls
		* 		 temp_data/processed_mill_geolocalization/noto.xls 

		* output:  input_data/manually_matched_ibs_uml/mills_to_georeference_post2010_done.xls
		* 		   input_data/manually_matched_ibs_uml/mills_to_georeference_pre2011_done.xlsx
		* 		   input_data/manually_matched_ibs_uml/noto_done.xls

	*** Automatic matching of remaining unreferenced mills. 
		rsource using "code/mill_geolocalization/match_unref.R"
		* input:	input_data/manufacturing_directories/direktori_industri_merged_cleaned.xlsx
		* 			temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		*			temp_data/processed_mill_geolocalization/ibs_unref.Rdata
		*			input_data/manually_matched_ibs_uml/matching_unref/wtn_cfl_done.xlsx
		*			input_data/manually_matched_ibs_uml/matching_unref/btw_cfl_done.xlsx

		* output:	temp_data/processed_mill_geolocalization/wtn_cfl.xlsx
		*			temp_data/processed_mill_geolocalization/btw_cfl.xlsx
		*			temp_data/processed_mill_geolocalization/named_unref.xlsx		

		* This named_unref.xlsx was sent to Jason Benedict who matched these MD names to UML names and outputed md_millMatching.xlsx
		* Some modifications were made manually to md_millMatching.xlsx, resulting in md_millMatching_modified.xlsx


	*** merge automatic and manual works - and resolve remaining conflicts.  
		do code/mill_geolocalization/merging_geolocalization_works.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		*		 input_data/manually_matched_ibs_uml/mills_to_georeference_post2010_done.xls 	 sheet("mills to georef")
		*  	     input_data/manually_matched_ibs_uml/mills_to_georeference_pre2011_done.xlsx	 sheet("Mills to georef")
		*  	     input_data/manually_matched_ibs_uml/noto_done.xls 	sheet("Sheet1")
		*  		 input_data/uml/traseMills_capEstyear.xlsx 	sheet("traseMills_capEstyear")
		*		 input_data/uml/mills_20200129.xlsx
		*		 input_data/manually_matched_ibs_uml/matching_unref/md_millMatching_modified.xlsx
		* 		 input_data/manually_matched_ibs_uml/overall_btw_conflicts_done.xlsx


		* output: temp_data/processed_mill_geolocalization/mills_to_georeference_post2010_done.dta 	 
		*  	      temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011_done.dta	 
		*  	      temp_data/processed_mill_geolocalization/noto_done.dta 	
		* 		  temp_data/processed_UML/traseMills_capEstyear_modified.dta
		*		  temp_data/processed_UML/mills_20200129_modified.dta
		* 		  temp_data/processed_mill_geolocalization/matching_unref/md_millMatching_modified.dta
		* 		  temp_data/processed_mill_geolocalization/merge_geoloc_works_temp.dta

		*  	      temp_data/processed_mill_geolocalization/IBS_UML_panel.dta
		*  	      temp_data/processed_mill_geolocalization/IBS_UML_panel.xlsx
		*  		  temp_data/processed_mill_geolocalization/IBS_UML_cs.dta


**** Make additional modifications to IBS, for which IBS-UML matching was necessary. 

	** impute establishment year from est_year variable in UML and min_year variable in IBS (requires both have been matched already)
	do code/explicative_variables/impute_UML_est_year.do
	* input		input_data/uml/mills_estyear_clean.xlsx
	*			temp_data/processed_UML/traseMills_capEstyear_modified.dta
	*			temp_data/processed_UML/mills_20200129_modified.dta
	*			temp_data/processed_mill_geolocalization/IBS_UML_cs.dta
	
	* output 	temp_data/processed_UML/mills_estyear_clean_modified.dta
	*			temp_data/processed_UML/UML_valentin_imputed_est_year.dta

	** compute concentration variable 
	rsource using "code/explicative_variables/mill_concentration.R"
	* input 	temp_data/processed_UML/UML_valentin_imputed_est_year.dta

	* output 	temp_data/processed_UML/UML_panel_valentin.dta

	** add international and domestic prices; compute log variables; add mill concentration and imputed establishment year variables 
	do code/explicative_variables/add_variables.do
	* input 	input_data/macro/prices_exp.xlsx		 
	*  	      	temp_data/processed_mill_geolocalization/IBS_UML_panel.dta
	*			temp_data/processed_UML/UML_panel_valentin.dta

	* output 	temp_data/processed_macro/prices_exp.dta
	*			temp_data/IBS_UML_panel_final.dta




**** build outcome variables 

	*** prepare annual maps of deforestation according to different definitions
		rsource using "code/outcome_variables/prepare_lucpfip.R" 
		* input:	temp_data/processed_indonesia_spatial/island_sf
		*			GFC tiles downloaded on the internet from the script
		* 			input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_2015_WGS1984.tif
		*			input_data/margono_primary_forest/timeseq_change00_12.tif
		*			temp_data/processed_mill_geolocalization/IBS_UML_cs.dta


		* output: 
		*		  TEMPORARY OUTPUTS in temp_data/processed_lu

		*		  /build/input/outcome_variables/annual_maps          Here are 54 .tif files: 2001-18 years for each three forest definition
		*         /build/input/outcome_variables/annual_parcels 	  For each PS, there are 54 .tif files: 2001-18 years for each three forest definition
		* 		  /build/input/outcome_variables/bricked_parcels   	  Years are bricked together so there is one .tif file for each forest definition. 

		*		MAIN OUTPUTS - each one has columns for 3 primary forest types (intact, degraded, total)
		* 		temp_data/processed_parcels/lucpfip_panel_3km_10CR.rds 
		* 		temp_data/processed_parcels/lucpfip_panel_3km_30CR.rds 
		* 		temp_data/processed_parcels/lucpfip_panel_3km_50CR.rds



**** build explicative variables at the parcel level

		/* Distribute geolocalized IBS mill variables to square parcels with distance-weighted averages.  
		 Parcels are grouped in dataframes according to their size (3x3km only currently - i.e. PS = 3000 (meters)) and 
		 how many they are (what is the catchment radius - CR = 10000, 30000, 50000 (meters)).  */
	
	**	
	wa_at_parcels.R
		* input:  /build/input/outcome_variables	take a panel_parcels file for arbitrary forest definition (say 30%) for a given PS and CR
		* 		  /build/output/IBS_UML_panel_final.dta

		* output:  /build/input/explanatory_variables/temp_cs_wa_explanatory/cs_wa_explanatory_   for each combination of year, PS and CR 	.rds
		* output:  /build/input/explanatory_variables/wa_panel_parcels_   	 for each combination of PS and CR 	.rds

	** add variables of the number of UML mills that are reachable from each parcel in each year within 10, 30 and 50km. 
	make_n_reachable_uml.R
		* input:	/build/output/UML_valentin_imputed_est_year.dta
		*			/build/input/explanatory_variables/wa_panel_parcels_   	 for each combination of PS and CR 	.rds

		*output: 	/build/input/explanatory_variables/panel_parcels_reachable_uml_ for each combination of PS and CR 	.rds

		* MERGE WITH OUTCOME VARIABLES !!!!!!!! 
		* in a dedicated script. 
		

		
* il reste à merger les RHS and LHS

///// analysis
