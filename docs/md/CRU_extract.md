# CRU 3.x Data Extraction Code



## Introduction

The `CRU_extract.R` script extracts CRU 3.x temperature and precipitation data from regional shapefiles and the full Alaska-Canada 2-km extent.
Optionally it also extracts data for specific point locations via the raster grid cell in which given spatial coordinates fall.
Data are saved to **R** workspaces (.RData files) for analysis and graphing by subsequent **R** code.

### Motivation
The primary motivation for this code is not just to extract regional and point data from a large number of high-resolution geotiffs,
but to limit the routine recurrence and redundancy of such extractions.
It allows for storing commonly required data in a more compact format that can be quickly shared and digested by other projects.

### Details

#### Capabilities
This script is set up to use parallel processing via the base `parallel` package and the function `mclapply` similar to its companion script, `AR4_AR5_extract.R`.
However, there are some subtle differences. First, it processes precipitation and temperature back to back, serially.
Between this and the fact that CRU 3.x is a single data set, compared to the multiple climate models and scenarios, there is nothing remaining to parallelize across.
The `mclapply` call only runs the extraction function on a single CPU core.
But in the future this code can easily be extended in parallel across more than just CRU 3.x, should other comparable data products become available.

The code ran be run strictly to obtain community data, regional data, or both. This allows for redoing only one type of extraction, for instance if new shapefiles or additional city locations are desired.
It can be passed an arbitrary list of cities and/or regions.

#### Limitations
This code is typically called by a SLURM script. It requires no command line arguments.
Although it uses only a single CPU core, it is still run as its own SLURM job, whereas it might be better to fold `CRU_extract.slurm` in with `AR4_AR5_extract.slurm`.
It is not set up to parallelize across climate variables.

Cities coordinates are assumed to be WGS84 Lat/Lon and are projected to NAD 83 Alaska Albers prior to use in extraction.
Regional information are expected in the form of a named nested list object of a specific type.
Non-list type sub-elements are vectors of cell indices previously extracted from 2-km Alaska-Canada extent template rasters using named groups of various shapefiles.
These objects are created in advance. There are multiple versions, as indexing by shapefile is based on the origin, resolution and extent of the rasters and geotiffs of intended use,
as well as subsampling and/or NA-removal, if desired.
The specific **R** workspace intended for use in data extraction by this script is hard coded and matches that used in `AR4_AR5_extract.R`.

### Files and data
The input files include 2-km Alaska-Canada extent data from 10 combined CMIP3 and CMIP5 downscaled GCMs
as well as community locations and raster layer cell indices pertaining to various regional shapefile polygons.
`CRU_extract.R` is called via SLURM script, `CRU_extract.slurm`.

The `CRU_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

Each of these groups of outputs is handled separately by additional **R** scripts.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R` which performs similar operations on downscaled GCM data.

## **R** Code

### Setup

Setup consists of loading required **R** packages and additional files, preparing any command line arguments for use, and defining functions.

#### Required packages


```r
library(raster)
library(maptools)
library(parallel)
library(plyr)
library(data.table)
```

#### Command line arguments


```r
comArgs <- commandArgs(TRUE)
if (length(comArgs)) for (i in 1:length(comArgs)) eval(parse(text = comArgs[[i]]))

if (!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if (!exists("regions")) regions <- FALSE
if (!exists("cities")) cities <- FALSE
if (!(regions | cities)) stop("regions and cities both FALSE. Nothing to process.")
if (!exists("cru")) cru <- "32"
```

#### Sourced code


```r
if (domain == "akcan2km") {
    # For regions and/or cities
    if (as.numeric(cru) == 31) 
        topDir <- file.path("/Data/Base_Data/Climate/AK_CAN_2km", paste0("historical/singleBand/CRU/cru_TS", 
            cru, "/historical"))
    if (as.numeric(cru) == 32) 
        topDir <- file.path("/Data/Base_Data/Climate/AK_CAN_2km", paste0("historical/CRU/CRU_TS", 
            cru))
    if (regions) {
        load("/workspace/UA/mfleonawicz/leonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN2km_5pct.RData")
    } else cells_shp_list_5pct <- region.names.out <- n.shp <- NULL
} else if (domain == "world10min") {
    # Currently for cities only
    topDir <- file.path("/Data/Base_Data/Climate/World/World_10min", paste0("historical/CRU/CRU_TS", 
        cru))  # files are not read, but metadata parsed from filenames list
    cells_shp_list_5pct <- region.names.out <- n.shp <- NULL
    #### Need to insert a load() command analogous to that above for regions
}

locs <- read.csv("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/locs.csv")
```

#### Define **R** objects


```r
if (exists("years")) yr1 <- years[1] else yr1 <- 1901
if (exists("years")) yr2 <- tail(years, 1) else if (as.numeric(cru) == 32) yr2 <- 2013 else yr2 <- 2009
years <- yr1:yr2
varid <- c("tas", "pr")

if (cities) {
    if (domain != "world10min") 
        locs <- subset(locs, region != "NWT")
    # locs <- locs[is.na(locs$pop) | locs$pop > 10,]
    l <- paste(locs$region, locs$loc)
    lu <- unique(l)
    dup <- which(duplicated(l))
    l.dup.u <- unique(l[dup])
    drp <- c()
    for (i in 1:length(l.dup.u)) {
        ind <- which(l == l.dup.u[i])
        ind <- ind[-which.max(locs$pop[ind])[1]]
        drp <- c(drp, ind)
    }
    cities <- locs[-drp, ]
    if (exists("cities.batch")) {
        batch.bounds <- round(seq(1, nrow(cities) + 1, length = 11))[c(cities.batch, 
            cities.batch + 1)] - c(0, 1)
        cities <- cities[batch.bounds[1]:batch.bounds[2], ]
    } else cities.batch <- 1
    d.cities <- cities[-c(which(names(locs) %in% c("lon_albers", "lat_albers")))]
    cities <- if (domain == "akcan2km") 
        cbind(cities$lon_albers, cities$lat_albers) else if (domain == "world10min") 
        cbind(cities$lon + 360, cities$lat)  # +360 for PC lat lon rasters
} else cities <- NULL

seasons <- "annual"
n.samples <- 100
n2 <- 2 * n.samples
agg.stat.colnames <- c("Mean", "SD", paste0("Pct_", c("05", 10, 25, 50, 75, 
    90, 95)))
agg.stat.names <- c("Mean", "Std Dev", paste0(c(5, 10, 25, 50, 75, 90, 95), 
    "th percentile"))
agg.stat.names[agg.stat.names == "50th percentile"] <- "Median"
```

#### Support functions


```r
# Density estimation
denFun <- function(x, n, adj = 0.25, variable) {
    x <- x[!is.na(x)]
    dif <- diff(range(x))
    z <- density(x, adjust = adj, n = n, from = min(x) - 0.05 * dif, to = max(x) + 
        0.05 * dif)
    if (variable == "pr" && any(z$x < 0)) 
        z <- density(x, adjust = adj, n = n, from = 0, to = max(x) + 0.05 * 
            dif)
    as.numeric(c(z$x, z$y))
}
```

----


```r
# Processing function
getData <- function(i, varid, cells.list = NULL, shp.names = NULL, n.shp = NULL, 
    seed = 232, regions = TRUE, n.samples = 512, cities = NULL, start.year = NULL, 
    end.year = NULL) {
    print(i)
    for (p in 1:length(varid)) {
        path <- file.path(topDir, varid[p])
        files <- list.files(path, pattern = ".tif$", full = T)
        if (!is.null(start.year)) 
            files <- files[substr(files, nchar(files) - 7, nchar(files) - 4) >= 
                start.year]
        if (!is.null(end.year)) 
            files <- files[substr(files, nchar(files) - 7, nchar(files) - 4) <= 
                end.year]
        
        if (!is.null(cities)) {
            r <- readAll(raster(files[1]))  # template done
            cells_cities <- extract(r, cities, cellnumbers = T)[, 1]
            print("Raster cell indices for point locations obtained.")
        }
        
        mo.tmp <- substr(files, nchar(files) - 10, nchar(files) - 9)
        yr.tmp <- substr(files, nchar(files) - 7, nchar(files) - 4)
        yr.mo.tmp <- paste0(yr.tmp, mo.tmp)
        files <- files[order(yr.mo.tmp)]
        yr.tmp <- substr(files, nchar(files) - 7, nchar(files) - 4)
        
        if (regions) {
            seq.q <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)  #seq(0, 1, by=0.05) # Hard coded
            m <- matrix(NA, nrow = length(files), ncol = length(seq.q) + 2)  # +2 for mean and SD below
            
            mden <- matrix(NA, nrow = 2 * n.samples, ncol = length(files))
            samples.list <- rapply(cells.list, f = function(x, m) m, classes = "integer", 
                how = "replace", m = mden)
            names(samples.list) <- names(shp.names)
            for (l1 in 1:length(shp.names)) names(samples.list[[l1]]) <- shp.names[[l1]]
        }
        
        if (regions) 
            for (l1 in 1:length(shp.names)) for (l2 in 1:length(shp.names[[l1]])) assign(paste("m", 
                names(shp.names)[l1], shp.names[[l1]][l2], sep = "__"), m)
        if (regions) 
            rm(m)
        gc()
        if (!is.null(cities)) 
            m.cities <- matrix(NA, nrow = length(cells_cities), ncol = length(files))
        for (b in 1:length(unique(yr.tmp))) {
            pat <- gsub("expression", "", paste(bquote(expression(".", .(unique(yr.tmp)[b]), 
                ".tif$")), collapse = ""))
            files.sub <- list.files(path, pattern = pat, full = T)
            mat <- getValues(stack(files.sub, quick = T))
            n <- ncol(mat)
            if (regions) {
                for (l1 in 1:length(shp.names)) {
                  for (l2 in 1:length(shp.names[[l1]])) {
                    cells.tmp <- cells.list[[names(shp.names)[l1]]][[shp.names[[l1]][l2]]]
                    m.tmp <- get(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], 
                      sep = "__"))
                    m.tmp[(1:n) + (12 * (b - 1)), 1] <- as.numeric(colMeans(mat[cells.tmp, 
                      ], na.rm = T))
                    m.tmp[(1:n) + (12 * (b - 1)), 2] <- as.numeric(apply(mat[cells.tmp, 
                      ], 2, sd, na.rm = T))
                    m.tmp[(1:n) + (12 * (b - 1)), 2 + (1:length(seq.q))] <- t(apply(mat[cells.tmp, 
                      ], 2, quantile, probs = seq.q, na.rm = T))
                    assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], 
                      sep = "__"), m.tmp)
                    samples.tmp <- apply(mat[cells.tmp, ], 2, denFun, n = n.samples, 
                      variable = varid[p])  # Time across columns
                    samples.list[[names(shp.names)[l1]]][[shp.names[[l1]][l2]]][, 
                      (1:n) + (12 * (b - 1))] <- samples.tmp
                  }
                }
                gc()
            }
            if (!is.null(cities)) 
                m.cities[, (1:n) + (12 * (b - 1))] <- mat[cells_cities, ]
            print(paste0("Process", i, ": ", unique(yr.tmp)[b]))
        }
        if (p == 1) {
            if (regions) {
                for (l1 in 1:length(shp.names)) for (l2 in 1:length(shp.names[[l1]])) assign(paste("m", 
                  names(shp.names)[l1], shp.names[[l1]][l2], "hold", sep = "__"), 
                  get(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], 
                    sep = "__")))
                samples.list.hold <- samples.list
            }
            if (!is.null(cities)) 
                m.cities.hold <- m.cities
            yr.hold <- as.numeric(yr.tmp)
        } else {
            if (regions) {
                for (l1 in 1:length(shp.names)) {
                  for (l2 in 1:length(shp.names[[l1]])) {
                    assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], 
                      "hold", sep = "__"), rbind(get(paste("m", names(shp.names)[l1], 
                      shp.names[[l1]][l2], "hold", sep = "__")), get(paste("m", 
                      names(shp.names)[l1], shp.names[[l1]][l2], sep = "__"))))
                  }
                }
                for (l1 in 1:length(shp.names)) {
                  for (l2 in 1:length(shp.names[[l1]])) {
                    samples.list.hold[[names(shp.names)[l1]]][[shp.names[[l1]][l2]]] <- cbind(samples.list.hold[[names(shp.names)[l1]]][[shp.names[[l1]][l2]]], 
                      samples.list[[names(shp.names)[l1]]][[shp.names[[l1]][l2]]])  # Time across columns
                  }
                }
            }
            if (!is.null(cities)) 
                m.cities.hold <- cbind(m.cities.hold, m.cities)
            yr.hold <- c(yr.hold, as.numeric(yr.tmp))
        }
    }
    if (regions) {
        output.names <- ls(pattern = "^m__.*.__hold$")
        m.list <- mget(output.names)
        regions.list <- c(m.list, list(yr = yr.hold), list(samples.list.hold))
        names(regions.list) <- c(sapply(strsplit(output.names, "__"), "[[", 
            3), "yr", "samples")
    }
    if (regions & is.null(cities)) 
        l <- list(regions.list, list("No cities"))
    if (regions & !is.null(cities)) 
        l <- list(regions.list, list(m.cities = m.cities.hold))
    if (!regions & !is.null(cities)) 
        l <- list(c(list(yr = yr.hold), list("empty")), list(m.cities = m.cities.hold))
    l
}
```

### Processing

Compile annual average time series data for each variable serially.
`mclapply` is called here in keeping with the method applied in the corresponding GCM extraction code found in `AR4_AR5_extract.R`, but defaults to serial processing.


```r
results <- mclapply(X = 1, FUN = getData, varid = varid, cells.list = cells_shp_list_5pct, 
    shp.names = region.names.out, n.shp = n.shp, regions = regions, n.samples = n.samples, 
    cities = cities, start.year = yr1, end.year = yr2, mc.cores = 1)
```

### Results


```r
# Organize and save results
for (k in 1:2) {
    if (!regions & k == 1) 
        next
    stats <- lapply(results, "[[", k)
    if (is.character(stats[[1]][[1]])) 
        next
    len <- ifelse(k == 1, 2, 0)
    if (k == 1) {
        # regional stats
        for (i in 1:length(stats)) {
            results.years <- results[[i]][[1]][[length(results[[i]][[1]]) - 
                1]]
            for (j in 1:(length(stats[[i]]) - len)) stats[[i]][[j]] <- stats[[i]][[j]][which(results.years %in% 
                years), ]  # matrix (regional stats)
        }
        stats.out.names <- names(results[[1]][[1]])
        stats.out.names <- stats.out.names[-which(stats.out.names %in% c("yr", 
            "samples"))]
        stats.out <- list()
        for (i in 1:length(stats.out.names)) stats.out[[i]] <- do.call(rbind, 
            lapply(stats, "[[", i))
        names(stats.out) <- stats.out.names
        time.seq <- rep(seq(years[1], length.out = nrow(stats.out[[1]])/(2 * 
            12)), each = 12)
        for (i in 1:length(stats.out)) {
            colnames(stats.out[[i]]) <- agg.stat.names
            rownames(stats.out[[i]]) <- rep(time.seq, 2)
            print(length(stats.out) - i)
        }
        save(stats.out, results.years, region.names.out, agg.stat.names, agg.stat.colnames, 
            file = paste0("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/regional/stats/", 
                "CRU", cru, "_annual_regions_stats.RData"))
        
        # regional samples
        samples <- lapply(results, function(x) x[[1]][["samples"]])
        samples <- samples[[1]]
        list2df <- function(x, yr_mo) {
            x <- ldply(x, data.frame)
            names(x) <- c("Location", yr_mo)
            x
        }
        for (l1 in 1:length(samples)) samples[[l1]] <- list2df(samples[[l1]], 
            yr_mo = paste(results.years, month.abb, sep = "_"))
        samples <- rbindlist(samples)
        samples.names <- unique(samples$Location)
        mid <- (ncol(samples) - 1)/2 + 1
        names.hold <- names(samples)[1:mid]
        samples <- data.frame(Location = samples$Location, rbind(as.matrix(samples[, 
            2:mid]), as.matrix(samples[, (mid + 1):ncol(samples)])))
        samples$Location <- rep(rep(samples.names, each = n2))
        samples.out <- list()
        for (i in 1:length(samples.names)) {
            samples.out[[i]] <- subset(samples, Location == samples.names[i])
            rownames(samples.out[[i]]) <- NULL
        }
        names(samples.out) <- samples.names
        save(samples.out, samples.names, region.names.out, n.samples, file = paste0("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/regional/samples/", 
            "CRU", cru, "_annual_regions_samples.RData"))
        
    }
    if (k == 2 & is.matrix(cities)) {
        for (i in 1:length(stats)) {
            results.years <- results[[i]][[1]][[length(results[[i]][[1]]) - 
                1]]
            for (j in 1:(length(stats[[i]]) - len)) stats[[i]][[j]] <- stats[[i]][[j]][, 
                which(results.years %in% years)]  # matrix (multiple cities)
        }
        stats.out <- list(cities = sapply(stats, function(x) t(x[[1]])))
        for (i in 1:length(stats.out)) rownames(stats.out[[i]]) <- rep(rep(rep(seq(yr1, 
            yr2), each = 12), nrow(cities)), length(varid))
        d <- stats.out[[1]]
        save(d, d.cities, results.years, file = paste0("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/cities/", 
            "CRU", cru, "_annual_cities_batch", cities.batch, "_", domain, ".RData"))
    }
}
```
