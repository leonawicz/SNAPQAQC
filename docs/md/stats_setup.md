# GCM Statistics Setup Code




## Introduction
The `stats_setup.R` script prepares extracted AR4/CMIP3 and AR5/CMIP5 regional temperature and precipitation data for direct use.

### Motivation
This script transforms and organizes output data from a parent script into a format that is easier to work with in certain frameworks.

### Details
The data are extracted based on regional shapefiles and the full Alaska-Canada 2-km extent for each of SNAP's five models, from both CMIP3 and CMIP5, and three emissions scenarios and RCPs, respectively.
Data manipulation is performed to organize the inputs into more usefully formatted outputs.


#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
It uses as many CPU cores as are available since it is parallelized across regions, which exceed the number of cores.
It becomes serial beyond the number of cores simultaneously in use.

#### Limitations
This code expects all ten models from CMIP3 and CMIP5 combined have been processed by the parent script.
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.

## Related items

### Files and Data
The input files are produced by the **R** script, `AR4_AR5_extract.R`. `stats_setup.R` assumes a complete and successful run of this precursory code.

The `AR4_AR5_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

For this reason, `stats_setup.R` has two companion **R** scripts which share the parent script, `AR4_AR5_extract.R`. These are `cities_setup.R` and `samples_setup.R`.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R`.
As a result, this script and its two companion scripts mentioned above are related to a similar set of three **R** scripts which perform similar intermediary data manipulation on extracted CRU 3.1 data.

## **R** Code

### Setup

Setup is minimal. Set working directory. List GCM files while making sure to avoid CRU 3.1 files.
SNAP's ten combined CMIP3 and CMIP5 models are in a hardcoded list as in the parent script.


```r
setwd("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/regional/stats")

library(data.table)

files <- list.files(pattern = "regions_stats.RData$")
files <- files[substr(files, 1, 3) != "CRU"]

models <- list(c("CCCMAcgcm31", "CCSM4"), c("GFDLcm21", "GFDLcm3"), c("MIROC32m", 
    "GISSe2-r"), c("MPIecham5", "IPSLcm5a-lr"), c("ukmoHADcm3", "MRIcgcm3"))
```

#### Load and compile data


```r
dlist <- vector("list", length(files))
for (i in 1:length(files)) {
    load(files[i])
    m <- do.call(rbind, stats.out)
    n <- nrow(m)
    p <- length(stats.out)
    dlist[[i]] <- data.frame(Phase = rep(c("CMIP3", "CMIP5"), each = n/(4 * 
        p)), Scenario = rep(c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", 
        "RCP 8.5"), each = n/(12 * p)), Model = rep(models[[i]], each = n/(4 * 
        p)), Var = rep(c("Temperature", "Precipitation"), each = n/(2 * p)), 
        Location = rep(names(stats.out), each = n/p), m, stringsAsFactors = F)
    print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(stats.out, n, m, i, p, files, models, dlist)
```

#### Organize data frame and support objects


```r
d[d$Var == "Temperature", 6:(ncol(d) - 1)] <- round(d[d$Var == "Temperature", 
    6:(ncol(d) - 1)], 1)
d[d$Var == "Precipitation", 6:(ncol(d) - 1)] <- round(d[d$Var == "Precipitation", 
    6:(ncol(d) - 1)])
d$Year <- results.years
rm(results.years)
d$Month <- month.abb
d$Month <- factor(d$Month, levels = month.abb)
d$Scenario <- factor(d$Scenario, levels = unique(d$Scenario))

phases <- unique(d$Phase)
years <- unique(d$Year)
decades <- years[years%%10 == 0]
d$Decade <- paste0(substr(d$Year, 1, 3), 0)
varnames <- unique(d$Var)
modnames <- lapply(1:length(phases), function(i, x, p) unique(x$Model[x$Phase == 
    p[i]]), x = d, p = phases)
scennames <- lapply(1:length(phases), function(i, x, p) unique(as.character(x$Scenario)[x$Phase == 
    p[i]]), x = d, p = phases)
stats.columns <- seq(which(names(d) == "Mean"), length.out = length(agg.stat.names))
stats.multiplier <- 1  # Multiply by 10 to shrink file size (doesn't seem to help with files this small)
# save.image('../../final/region_stats_data.RData')
```

#### Save output files


```r
library(parallel)
f <- function(i, index, multiplier) {
    name.tmp <- as.character(unlist(region.names.out))[i]
    region.dat <- subset(d, Location == name.tmp)
    region.dat <- multiplier * unlist(region.dat[, index])
    names(region.dat) <- NULL
    grp <- rep(names(region.names.out), times = sapply(region.names.out, length))[i]
    dir.create(outDir <- file.path("../../final/region_files_GCM/stats", grp, 
        name.tmp), recursive = T, showWarnings = F)
    save(region.dat, file = file.path(outDir, "stats_climate.RData"))
    print(i)
}

mclapply(1:sum(sapply(region.names.out, length)), f, index = stats.columns, 
    multiplier = stats.multiplier, mc.cores = 32)

gcm.stats.df <- subset(d, Location == as.character(unlist(region.names.out))[1])
gcm.stats.df[, stats.columns] <- NA
gcm.stats.df$Location <- NA
stats.colnames <- names(gcm.stats.df)[stats.columns] <- agg.stat.colnames
rm(d, f)
save.image("../../final/meta.RData")
```
