# GCM Cities Setup Code



## Introduction
The `cities_setup.R` script prepares extracted AR4/CMIP3 and AR5/CMIP5 community temperature and precipitation data for direct use.

### Motivation
This script transforms and organizes output data from a parent script into a format that is easier to work with in certain frameworks.

### Details
The data are extracted based on coordinates of communities.
Using the 2-km Alaska_Canada extent, single 2-km raster grid cells in which coordinate pairs for each community fall are used as point location estimates.
This is done for each of SNAP's five models, from both CMIP3 and CMIP5, and three emissions scenarios and RCPs, respectively.
Data manipulation is performed to organize the inputs into more usefully formatted outputs.

#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
It uses as many CPU cores as are available since it is parallelized across cities, which exceed the number of cores.
It becomes serial beyond the number of cores simultaneously in use.

#### Limitations
This code expects all ten models from CMIP3 and CMIP5 combined have been processed by the parent script.
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.
This script is called by `cities_setup.slurm`.
It is suggested that for large numbers of cities, batch processing be used, passing one batch of input files to this script via `cities_setup.slurm` at a time at the command line.
In such a case, there will be sets of input .RData files covering all the models and scenarios, for each batch, where a batch includes a subset of the cities.
It may be too much to include several thousand cities in a single batch, depending on other parameters of the upstream extractions such as the time period.

## Related items

### Files and data
The input files are produced by the **R** script, `AR4_AR5_extract.R`. `cities_setup.R` assumes a complete and successful run of this precursory code.

The `AR4_AR5_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

For this reason, `cities_setup.R` has two companion **R** scripts which share the parent script, `AR4_AR5_extract.R`. These are `stats_setup.R` and `samples_setup.R`.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R`.
As a result, this script and its two companion scripts mentioned above are related to a similar set of three **R** scripts which perform similar intermediary data manipulation on extracted CRU 3.x data.

## **R** Code

### Setup

Setup is minimal. Set working directory. List GCM files while making sure to avoid CRU 3.x files.
SNAP's ten combined CMIP3 and CMIP5 models are in a hardcoded list as in the parent script.
Cities may be processed in batches via command line argument.


```r
comArgs <- commandArgs(TRUE)
if (length(comArgs)) for (i in 1:length(comArgs)) eval(parse(text = comArgs[[i]]))
if (!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if (!exists("cities.batch")) cities.batch <- 1

setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/cities")

library(data.table)

files <- list.files(pattern = paste0("batch", cities.batch, "_", domain, ".RData$"))
files <- files[substr(files, 1, 3) != "CRU"]

models <- list(c("CCCMAcgcm31", "CCSM4"), c("GFDLcm21", "GFDLcm3"), c("MIROC32m", 
    "GISSe2-r"), c("MPIecham5", "IPSLcm5a-lr"), c("ukmoHADcm3", "MRIcgcm3"))
```

#### Load and compile data


```r
dlist <- vector("list", length(files))
for (i in 1:length(files)) {
    load(files[i])
    cities.meta <- d.cities
    m <- as.numeric(d)
    n <- length(m)
    dlist[[i]] <- data.frame(Phase = rep(c("AR4", "AR5"), each = n/4), Scenario = rep(c("SRES B1", 
        "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5"), each = n/12), 
        Model = rep(models[[i]], each = n/4), Var = rep(c("Temperature", "Precipitation"), 
            each = n/2), Location = rep(paste0(cities.meta$loc, ", ", cities.meta$region), 
            each = n/(12 * nrow(cities.meta))), Val = m, stringsAsFactors = F)
    gc()
    print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(n, m, i, files, models, dlist)
gc()
```

#### Organize data frame and support objects


```r
d$Val[d$Var == "Temperature"] <- round(d$Val[d$Var == "Temperature"], 1)
d$Val[d$Var == "Precipitation"] <- round(d$Val[d$Var == "Precipitation"])
d$Month <- month.abb
d$Year <- results.years  #rep(1870:2099,each=12)
d$Month <- factor(d$Month, levels = d$Month[1:12])
d$Scenario <- factor(d$Scenario, levels = unique(d$Scenario))
cities.meta$loc <- paste0(cities.meta$loc, ", ", cities.meta$region)
cities.meta <- cities.meta[c(1, 3:6)]
names(cities.meta) <- c("Country", "Location", "Population", "Lat", "Lon")
d$Decade <- as.character(10 * (d$Year%/%10))
gc()
d.cities <- d
rm(d, results.years)
gc()
# save.image('../final/data_cities.RData')
```

#### Save output files


```r
# Save individual R workspace files for each city for use in master QAQC
# Shiny app
library(parallel)

f <- function(i, overwrite = FALSE) {
    name.tmp <- gsub("\\.", "PER", gsub("/", "FSLASH", gsub("`", "", gsub("~", 
        "", gsub("?", "", gsub("\\'", "APOS", cities.meta$Location[i]))))))
    filename <- paste0("../final/city_files_GCM/", domain, "/", gsub(", ", "--", 
        name.tmp), "__", domain, ".RData")
    if (overwrite | !file.exists(filename)) {
        city.dat <- subset(d.cities, Location == cities.meta$Location[i])
        save(city.dat, file = filename)
    }
    print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores = 32)
```
