# @knitr regions
#### This is the final script run after all app data and metadata files are compiled.
#### This is for final formatting related to regions and cities.

topDir <- "/workspace/Shared/Users/mfleonawicz/AR4_AR5_Comparisons"

region.gcm.stats.path <- list.files(file.path(topDir, "region_files_GCM/stats"), full=T)
region.gcm.stats.files <- lapply(region.gcm.stats.path, list.files, full=T)
region.cru.stats.path <- list.files(file.path(topDir, "region_files_CRU/stats"), full=T)
region.cru.stats.files <- lapply(region.cru.stats.path, list.files, full=T)

region.gcm.samples.path <- list.files(file.path(topDir, "region_files_GCM/samples"), full=T)
region.gcm.samples.files <- lapply(region.gcm.samples.path, list.files, full=T)
region.cru.samples.path <- list.files(file.path(topDir, "region_files_CRU/samples"), full=T)
region.cru.samples.files <- lapply(region.cru.samples.path, list.files, full=T)
names(region.cru.stats.files) <- names(region.gcm.stats.files) <- names(region.cru.samples.files) <- names(region.gcm.samples.files) <- basename(region.gcm.stats.path)

# @knitr cities
city.gcm.files.path <- file.path(topDir, "city_files_GCM")
city.gcm.files <- list.files(city.gcm.files.path, full=T)
city.cru.files.path <- file.path(topDir, "city_files_CRU")
city.cru.files <- list.files(city.cru.files.path, full=T)

# city names to appear in app menu
city.names <- gsub("FSLASH", "/", gsub("PER", "\\.", gsub("APOS", "\\'", gsub("--", ", ", sapply(strsplit(basename(city.gcm.files), "__"), "[[", 1)))))

# subset cities metadata data frame to match final files
cities.meta.file <- "/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/cities_meta_akcan2km.RData"
load(cities.meta.file)
cities.meta.akcan2km <- subset(cities.meta.akcan2km, Location %in% city.names)
save(cities.meta.akcan2km, file=cities.meta.file) # okay to save over original file
rm(cities.meta.file, topDir)

# @knitr save_metadata
load("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/meta.RData")
save.image(file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/meta.RData"))
