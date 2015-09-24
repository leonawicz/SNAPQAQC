# @knitr regions
#### This is the final script run after all app data and metadata files are compiled.
#### This is for final formatting related to regions and cities.

load("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final/meta.RData")

topDir <- "/workspace/Shared/Users/mfleonawicz/qaqcAppData"

region.gcm.stats.path <- list.files(file.path(topDir, "region_files_GCM/stats"), full=T)
region.gcm.stats.files <- lapply(region.gcm.stats.path, list.files, full=T)
region.cru.stats.path <- list.files(file.path(topDir, "region_files_CRU32/stats"), full=T)
region.cru.stats.files <- lapply(region.cru.stats.path, list.files, full=T)

region.gcm.samples.path <- list.files(file.path(topDir, "region_files_GCM/samples"), full=T)
region.gcm.samples.files <- lapply(region.gcm.samples.path, list.files, full=T)
region.cru.samples.path <- list.files(file.path(topDir, "region_files_CRU32/samples"), full=T)
region.cru.samples.files <- lapply(region.cru.samples.path, list.files, full=T)
names(region.cru.stats.files) <- names(region.gcm.stats.files) <- names(region.cru.samples.files) <- names(region.gcm.samples.files) <- basename(region.gcm.stats.path)

# @knitr cities
city.gcm.files.path <- file.path(topDir, "city_files_GCM")
city.gcm.files.2km <- list.files(file.path(city.gcm.files.path, "akcan2km"), full=T)
city.gcm.files.10min <- list.files(file.path(city.gcm.files.path, "world10min"), full=T)
city.cru.files.path <- file.path(topDir, "city_files_CRU32")
city.cru.files.2km <- list.files(file.path(city.cru.files.path, "akcan2km"), full=T)
city.cru.files.10min <- list.files(file.path(city.cru.files.path, "world10min"), full=T)

# city names to appear in app menu
# use one group of cities from above, all four sets are identical
city.names <- gsub("FSLASH", "/", gsub("PER", "\\.", gsub("APOS", "\\'", gsub("--", ", ", sapply(strsplit(basename(city.gcm.files.10min), "__"), "[[", 1)))))

# subset cities metadata data frame to match final files, same for all sets
cities.meta.file <- "/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final/cities_meta.RData"
load(cities.meta.file)
cities.meta <- subset(cities.meta, Location %in% city.names)
save(cities.meta, file=cities.meta.file) # okay to save over original file
rm(cities.meta.file, topDir)

# @knitr save_metadata
save.image(file.path("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final/meta.RData"))
