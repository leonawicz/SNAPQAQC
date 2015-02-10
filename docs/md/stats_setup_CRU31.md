# CRU 3.1 Statistics Setup Code




## Introduction
The `stats_setup_CRU31.R` script prepares extracted CRU 3.1 regional temperature and precipitation data for direct use.

### Motivation
This script transforms and organizes output data from a parent script into a format that is easier to work with in certain frameworks.

### Details
The data are extracted based on regional shapefiles and the full Alaska-Canada 2-km extent for SNAP's CRU 3.1 data set.
Data manipulation is performed to organize the inputs into more usefully formatted outputs.


#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
It uses as many CPU cores as are available since it is parallelized across regions, which exceed the number of cores.
It becomes serial beyond the number of cores simultaneously in use.

#### Limitations
This code expects CRU 3.1 data have been processed by the parent script.
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.

## Related items

### Files and Data
The input files are produced by the **R** script, `CRU_extract.R`. `stats_setup_CRU31.R` assumes a complete and successful run of this precursory code.

The `CRU_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

For this reason, `stats_setup_CRU31.R` has two companion **R** scripts which share the parent script, `CRU_extract.R`. These are `cities_setup_CRU31.R` and `samples_setup_CRU31.R`.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R`.
As a result, this script and its two companion scripts mentioned above are related to a similar set of three **R** scripts which perform similar intermediary data manipulation on extracted GCM data.

## **R** Code

### Setup

Setup is minimal. Set working directory. List CRU 3.1 files.

#### Load and compile data


```r
dlist <- vector("list", length(files))
for (i in 1:length(files)) {
    load(files[i])
    m <- do.call(rbind, stats.out)
    n <- nrow(m)
    p <- length(stats.out)
    dlist[[i]] <- data.frame(Var = rep(c("Temperature", "Precipitation"), each = n/(2 * 
        p)), Location = rep(names(stats.out), each = n/p), m, stringsAsFactors = F)
    print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(stats.out, n, m, i, p, files, models, dlist)
```

#### Organize data frame and support objects


```r
d[d$Var == "Temperature", 3:ncol(d)] <- round(d[d$Var == "Temperature", 3:ncol(d)], 
    1)
d[d$Var == "Precipitation", 3:ncol(d)] <- round(d[d$Var == "Precipitation", 
    3:ncol(d)])
d$Year <- results.years
rm(results.years)
d$Month <- month.abb
d$Month <- factor(d$Month, levels = month.abb)
d$Decade <- paste0(substr(d$Year, 1, 3), 0)

# agg.stat.names <- c('Mean', 'Std. Dev.', '5th percentile','10th
# percentile', '25th percentile', '50th (Median)', '75th percentile', '90th
# percentile', '95th percentile')
stats.columns.cru <- seq(which(names(d) == "Mean"), length.out = length(agg.stat.names))
# save.image('../../final/CRU31_region_stats_data.RData')
```

#### Save output files


```r
library(parallel)

f <- function(i) {
    name.tmp <- as.character(unlist(region.names.out))[i]
    # assign(name.tmp, subset(d.cities, Location==cities.meta$Location[i]))
    # save(list=c(name.tmp), file=paste0('../final/city_files_GCM/', gsub(', ',
    # '--', name.tmp), '.RData'))
    region.cru.dat <- subset(d, Location == name.tmp)
    names(region.cru.dat)[stats.columns.cru] <- agg.stat.colnames
    grp <- rep(names(region.names.out), times = sapply(region.names.out, length))[i]
    dir.create(outDir <- file.path("../../final/region_files_CRU/stats", grp, 
        name.tmp), recursive = T, showWarnings = F)
    save(region.cru.dat, file = file.path(outDir, "stats_climate.RData"))
    print(i)
}

mclapply(1:length(unique(d$Location)), f, mc.cores = 32)
```

NEW_CODE_CHUNKS

```r
setwd("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/regional/stats")

library(data.table)

files <- list.files(pattern = "^CRU31.*.regions_stats.RData$")

models <- "CRU31"
```
