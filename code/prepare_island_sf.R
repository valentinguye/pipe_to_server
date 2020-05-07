##### prepare shapefile of Indonesian Islands. 


### PACKAGES ###

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "sf")

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Load them
lapply(neededPackages, library, character.only = TRUE)

### /!\ IF renv::restore() FAILS TO INSTALL SOME PACKAGES FROM neededPackages /!\ ### 

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


### Aggregate Indonesian province shapes in islands. 
provinces <- st_read("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp")
provinces <- dplyr::select(provinces, NAME_1)

provinces$island <- NA

provinces$island[provinces$NAME_1 == "Aceh" |
                        provinces$NAME_1 == "Bangka-Belitung" |
                        provinces$NAME_1 == "Bengkulu" |
                        provinces$NAME_1 == "Jambi" |
                        provinces$NAME_1 == "Kepulauan Riau" |
                        provinces$NAME_1 == "Lampung" |
                        provinces$NAME_1 == "Riau" |
                        provinces$NAME_1 == "Sumatera Barat" |
                        provinces$NAME_1 == "Sumatera Selatan" |
                        provinces$NAME_1 == "Sumatera Utara" ] <- "Sumatra"

provinces$island[provinces$NAME_1 == "Kalimantan Barat" |
                        provinces$NAME_1 == "Kalimantan Selatan" |
                        provinces$NAME_1 == "Kalimantan Tengah" |
                        provinces$NAME_1 == "Kalimantan Timur" |
                        provinces$NAME_1 == "Kalimantan Utara" ] <- "Kalimantan"

provinces$island[provinces$NAME_1 == "Papua" |
                        provinces$NAME_1 == "Irian Jaya Barat" ] <- "Papua"

island_sf <- provinces[!is.na(provinces$island),c("island", "geometry")]

IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  island_sf$geometry[island_sf$island == Island] <- st_union(island_sf$geometry[island_sf$island == Island])
}

island_sf <- island_sf[!duplicated(island_sf$island),]

## export 
dir.create("temp_data/processed_indonesia_spatial/island_sf")
st_write(island_sf, "temp_data/processed_indonesia_spatial/island_sf", driver = "ESRI Shapefile", delete_dsn = TRUE)

