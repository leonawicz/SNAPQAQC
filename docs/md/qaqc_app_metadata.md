# QAQC Shiny App Metadata Finalization



The `qaqc_app_metadata.R` script loads the `meta.RData` metadata workspace file used by the `cmip3_cmip5` QAQC app.
Final preparations of region and city data file path locations are completed.
This new workspace is saved over the original.

This script current handles subsetting of city files as best it can,
but ultimately time is needed to investigate why GCM and CRU extractions of the same set of cities (specific grid cells) are yielding different, overlapping sets of cities.

## R code

File path information for GCM and CRU files are finalized.


```r
#### This is the final script run after all app data and metadata files are
#### compiled.  This is for final formatting related to regions and cities.

topDir <- "/workspace/Shared/Users/mfleonawicz/AR4_AR5_Comparisons"

region.gcm.stats.path <- list.files(file.path(topDir, "region_files_GCM/stats"), 
    full = T)
region.gcm.stats.files <- lapply(region.gcm.stats.path, list.files, full = T)
region.cru.stats.path <- list.files(file.path(topDir, "region_files_CRU/stats"), 
    full = T)
region.cru.stats.files <- lapply(region.cru.stats.path, list.files, full = T)

region.gcm.samples.path <- list.files(file.path(topDir, "region_files_GCM/samples"), 
    full = T)
region.gcm.samples.files <- lapply(region.gcm.samples.path, list.files, full = T)
region.cru.samples.path <- list.files(file.path(topDir, "region_files_CRU/samples"), 
    full = T)
region.cru.samples.files <- lapply(region.cru.samples.path, list.files, full = T)
names(region.cru.stats.files) <- names(region.gcm.stats.files) <- names(region.cru.samples.files) <- names(region.gcm.samples.files) <- basename(region.gcm.stats.path)
```

File path information for GCM and CRU cities, and associated metadata, are fortified as best as currently possible.
See code comments for details regarding pending investigation of discrepancies between data sets.


```r
# City inclusion is not perfect, some cities are dropped out CRU has fewer
# city files, no time to investigate yet
city.cru.files.path <- file.path(topDir, "city_files_CRU")
city.cru.files <- list.files(city.cru.files.path, full = T)

city.gcm.files.path <- file.path(topDir, "city_files_GCM")
# restrict GCM files to subset with corresponding CRU files
city.gcm.files <- gsub("city_files_CRU", "city_files_GCM", city.cru.files)  # list.files(city.gcm.files.path, full=T)

# city names to appear in app menu
city.names <- gsub("APOS", "\\'", gsub("--", ", ", sapply(strsplit(basename(city.gcm.files), 
    "__"), "[[", 1)))

# subset cities metadata data frame to match CRU files
cities.meta.file <- "/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/cities_meta_akcan2km.RData"
load(cities.meta.file)
cities.meta.akcan2km <- subset(cities.meta.akcan2km, Location %in% city.names)
save(cities.meta.akcan2km, file = cities.meta.file)  # okay to save over original file
rm(cities.meta.file, topDir)

# reverse for consistency, extra precaution
city.names <- city.names[city.names %in% cities.meta.akcan2km$Location]
```

The `meta.RData` workspace used by the master QAQC Shiny app is loaded and re-saved.


```r
load("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/meta.RData")
save.image(file.path("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/meta.RData"))
```
