# AR4/CMIP3 and AR5/CMIP5 Data Extraction Code



## Introduction

The `AR4_AR5_extract.R` script extracts AR4/CMIP3 and AR5/CMIP5 temperature and precipitation data from regional shapefiles and the full Alaska-Canada 2-km extent.
Optionally it also extracts data for specific point locations via the raster grid cell in which given spatial coordinates fall.
Data are saved to **R** workspaces (.RData files) for analysis and graphing by subsequent **R** code.

### Motivation
The primary motivation for this code is not just to extract regional and point data from a large number of high-resolution geotiffs,
but to limit the routine recurrence and redundancy of such extractions.
It allows for storing commonly required data in a more compact format that can be quickly shared and digested by other projects.

### Details

#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
It uses 12 CPU cores to process all 12 combinations of climate variable (precipitation and temperature) and emissions scenario (three per CMIP phase).
These 12 processes are applied to one GCM at a time from each of the two phases.
This is why they are listed in the code in pairs of models, one from each phase.
Processing is serial across model pairs.

The code ran be run strictly to obtain community data, regional data, or both. This allows for redoing only one type of extraction, for instance if new shapefiles or additional city locations are desired.
It can be passed an arbitrary list of cities and/or regions.

#### Limitations
This code is intended to be run on all ten models hard coded in the list of pairs below.
It can be run on fewer and in fact is typically called via SLURM script with a command line argument referring to a specific batch pair.
Outputs for each pair of two models are saved independently.
However, subsequent **R** scripts assume that all ten GCM data extractions have been completed and that all necessary output files are available for further processing.
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.

Cities coordinates are assumed to be WGS84 Lat/Lon and are projected to NAD 83 Alaska Albers prior to use in extraction.
Regional information are expected in the form of a named nested list object of a specific type.
Non-list type sub-elements are vectors of cell indices previously extracted from 2-km Alaska-Canada extent template rasters using named groups of various shapefiles.
These objects are created in advance. There are multiple versions, as indexing by shapefile is based on the origin, resolution and extent of the rasters and geotiffs of intended use,
as well as subsampling and/or NA-removal, if desired.
The specific **R** workspace intended for use in data extraction by this script is hard coded and matches that used in `AR4_AR5_extract.R`.

### Files and data
The input files include 2-km Alaska-Canada extent data from 10 combined CMIP3 and CMIP5 downscaled GCMs
as well as community locations and raster layer cell indices pertaining to various regional shapefile polygons.
`AR4_AR5_extract.R` is called via SLURM script, `AR4_AR5_extract.slurm`.

The `AR4_AR5_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

Each of these groups of outputs is handled separately by additional **R** scripts.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R` which performs similar operations on downscaled CRU 3.x data.

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

if (!exists("batch")) stop("Model pair batch integer not supplied.")
if (!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if (!exists("regions")) regions <- FALSE
if (!exists("cities")) cities <- FALSE
if (!(regions | cities)) stop("regions and cities both FALSE. Nothing to process.")
```

#### Sourced code


```r
if (domain == "akcan2km") {
    # For regions and/or cities
    topDir <- file.path("/Data/Base_Data/Climate/AK_CAN_2km", c("historical", 
        "projected"))  # files are not read, but metadata parsed from filenames list
    if (regions) {
        load("/workspace/UA/mfleonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN2km_5pct.RData")
    } else cells_shp_list_5pct <- region.names.out <- n.shp <- NULL
} else if (domain == "world10min") {
    # Currently for cities only
    topDir <- file.path("/Data/Base_Data/Climate/World/World_10min", c("historical", 
        "projected"))  # files are not read, but metadata parsed from filenames list
    cells_shp_list_5pct <- region.names.out <- n.shp <- NULL
    #### Need to insert a load() command analogous to that above for regions
}

locs <- read.csv("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/locs.csv")
```

#### Define **R** objects


```r
varid <- rep(c("tas", "pr"), each = 6)
scenid <- if (length(topDir) == 2) {
    c(rep("", 12), rep(c("sresb1", "sresa1b", "sresa2", "rcp45", "rcp60", "rcp85"), 
        2))
} else if (basename(topDir) == "projected") {
    rep(c("sresb1", "sresa1b", "sresa2", "rcp45", "rcp60", "rcp85"), 2)
} else if (basename(topDir) == "historical") {
    rep("", 12)
}
arid <- rep(rep(c("AR4_CMIP3_models", "AR5_CMIP5_models"), each = 3), 2)

model.pairs <- list(c("cccma-cgcm3-1-t47", "CCSM4"), c("gfdl-cm2-1", "GFDL-CM3"), 
    c("miroc3-2-medres", "GISS-E2-R"), c("mpi-echam5", "IPSL-CM5A-LR"), c("ukmo-hadcm3", 
        "MRI-CGCM3"))

if (exists("years")) yr1 <- years[1] else yr1 <- 1870
if (exists("years")) yr2 <- tail(years, 1) else yr2 <- 2099
years <- yr1:yr2

if (cities) {
    if (domain != "world10min") 
        locs <- subset(locs, region != "Northwest Territories")
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
    if (domain == "world10min") {
        # hack to deal with specific NT locations
        nt.na <- which(cities$loc %in% c("Paulatuk", "Sachs Harbour"))
        cities$lon[nt.na] <- cities$lon[nt.na] + 0.1666667
    }
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
getData <- function(i, model, cells.list = NULL, shp.names = NULL, n.shp = NULL, 
    seed = 232, regions = TRUE, n.samples = 512, cities = NULL, start.year = NULL, 
    end.year = NULL) {
    print(i)
    if (length(model) == 1) 
        model <- rep(model, i)
    for (p in 1:length(topDir)) {
        if (p == 2) 
            print(paste("Process", i, "beginning projected data extraction."))
        if (p == 1) 
            scenid.tmp <- scenid[1:12] else if (p == 2) 
            scenid.tmp <- scenid[13:24]
        path <- gsub("/remove", "", gsub("//", "/", file.path(topDir[p], arid[i], 
            scenid.tmp[i], model[i], varid[i])))
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
        yr.tmp <- substr(files, nchar(files) - 7, nchar(files) - 4)
        if (regions) {
            seq.q <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)  #seq(0, 1, by=0.05) #Hard coded
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
        if ((i %in% c(1, 4, 7, 10)) | p == 2) {
            # no need to read historical files more than once for each CMIP phase,
            # scenario, and variable.
            for (b in 1:length(unique(yr.tmp))) {
                pat <- gsub("expression", "", paste(bquote(expression(".", .(unique(yr.tmp)[b]), 
                  ".tif$")), collapse = ""))
                files.sub <- list.files(path, pattern = pat, full = T)
                mat <- getValues(stack(files.sub))
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
                        variable = varid[i])  # Time across columns
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
        } else print(paste("Process", i, "skipping historical files to avoid redundancy."))
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

Compile annual average time series data for each model, scenario/rcp, variable, from CMIP3 and CMIP5, two models at a time (one from each phase)
**batch** is passed as a command line argument


```r
models <- rep(rep(model.pairs[[batch]], each = 3), 2)
results <- mclapply(X = 1:12, FUN = getData, model = models, cells.list = cells_shp_list_5pct, 
    shp.names = region.names.out, n.shp = n.shp, regions = regions, n.samples = n.samples, 
    cities = cities, start.year = yr1, end.year = yr2, mc.cores = 12)
```

### Results


```r
# Organize and save results
for (k in 1:2) {
    if (!regions & k == 1) 
        next
    stats <- lapply(results, "[[", k)
    if (is.character(stats[[1]])) 
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
        time.seq <- rep(seq(years[1], length.out = nrow(stats.out[[1]])/(12 * 
            12)), each = 12)
        for (i in 1:length(stats.out)) {
            colnames(stats.out[[i]]) <- agg.stat.colnames
            rownames(stats.out[[i]]) <- rep(time.seq, 12)
            for (block in c(2, 5, 8, 11)) {
                block.rows <- 1:length(time.seq) + length(time.seq) * (block - 
                  1)
                na.rows <- which(is.na(stats.out[[i]][block.rows, 1]))
                stats.out[[i]][block.rows[na.rows], ] <- stats.out[[i]][block.rows[na.rows] + 
                  length(time.seq), ] <- stats.out[[i]][block.rows[na.rows] - 
                  length(time.seq), ]
            }
            print(length(stats.out) - i)
        }
        save(stats.out, results.years, region.names.out, agg.stat.names, agg.stat.colnames, 
            file = paste0("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/regional/stats/", 
                models[1], "_", models[4], "_annual_regions_stats.RData"))
        
        # regional samples
        samples <- lapply(results, function(x) x[[1]][["samples"]])
        list2df <- function(x, yr_mo) {
            x <- ldply(x, data.frame)
            names(x) <- c("Location", yr_mo)
            x
        }
        for (z in 1:length(arid)) {
            # 12 combinations
            lab <- paste(substr(arid[z], 1, 3), scenid[scenid != ""][z], models[z], 
                varid[z], sep = "_")
            samples.sub <- samples[[z]]
            for (l1 in 1:length(samples.sub)) samples.sub[[l1]] <- list2df(samples.sub[[l1]], 
                yr_mo = paste(results.years, month.abb, sep = "_"))
            samples.sub <- rbindlist(samples.sub)
            samples.names <- unique(samples.sub$Location)
            samples.out <- list()
            for (i in 1:length(samples.names)) {
                samples.out[[i]] <- subset(samples.sub, Location == samples.names[i])
                rownames(samples.out[[i]]) <- NULL
            }
            names(samples.out) <- samples.names
            save(samples.out, samples.names, region.names.out, n.samples, file = paste0("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/regional/samples/", 
                lab, "_annual_regions_samples.RData"))
        }
    }
    if (k == 2 & !is.null(cities)) {
        for (i in 1:length(stats)) {
            results.years <- results[[i]][[1]][[length(results[[i]][[1]]) - 
                1]]
            for (j in 1:(length(stats[[i]]) - len)) stats[[i]][[j]] <- stats[[i]][[j]][, 
                which(results.years %in% years)]  # matrix (multiple cities)
        }
        stats.out <- list(cities = sapply(stats, function(x) t(x[[1]])))
        for (i in 1:length(stats.out)) {
            rownames(stats.out[[i]]) <- rep(rep(seq(yr1, yr2), each = 12), nrow(cities))
            na.rows1 <- which(is.na(stats.out[[i]][, 2]))
            na.rows2 <- which(is.na(stats.out[[i]][, 5]))
            stats.out[[i]][na.rows1, 2] <- stats.out[[i]][na.rows1, 3] <- stats.out[[i]][na.rows1, 
                1]
            stats.out[[i]][na.rows2, 5] <- stats.out[[i]][na.rows2, 6] <- stats.out[[i]][na.rows2, 
                4]
            stats.out[[i]][na.rows1, 8] <- stats.out[[i]][na.rows1, 9] <- stats.out[[i]][na.rows1, 
                7]
            stats.out[[i]][na.rows2, 11] <- stats.out[[i]][na.rows2, 12] <- stats.out[[i]][na.rows2, 
                10]
        }
        d <- stats.out[[1]]
        save(d, d.cities, results.years, file = paste0("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/cities/", 
            models[1], "_", models[4], "_annual_cities_batch", cities.batch, 
            "_", domain, ".RData"))
    }
}
```
