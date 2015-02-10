# CRU 3.1 Cities Setup Code
Matthew Leonawicz  



## Introduction
The `cities_setup_CRU31.R` script prepares extracted CRU 3.1 community temperature and precipitation data for direct use.

### Motivation
This script transforms and organizes output data from a parent script into a format that is easier to work with in certain frameworks.

### Details
The data are extracted based on coordinates of communities.
Using the 2-km Alaska_Canada extent, single 2-km raster grid cells in which coordinate pairs for each community fall are used as point location estimates.
This is done for SNAP's CRU 3.1 data set.
Data manipulation is performed to organize the inputs into more usefully formatted outputs.

#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
It uses as many CPU cores as are available since it is parallelized across cities, which exceed the number of cores.
It becomes serial beyond the number of cores simultaneously in use.

#### Limitations
This code expects CRU 3.1 data have been processed by the parent script.
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.

## Related items

### Files and Data
The input files are produced by the **R** script, `CRU_extract.R`. `cities_setup_CRU31.R` assumes a complete and successful run of this precursory code.

The `CRU_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

For this reason, `cities_setup_CRU31.R` has two companion **R** scripts which share the parent script, `CRU_extract.R`. These are `stats_setup_CRU31.R` and `samples_setup_CRU31.R`.

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
    cities.meta <- d.cities
    m <- as.numeric(d)
    n <- length(m)
    dlist[[i]] <- data.frame(Var = rep(c("Temperature", "Precipitation"), each = n/(2 * 
        nrow(cities.meta))), Location = rep(paste0(cities.meta$loc, ", ", cities.meta$region), 
        each = n/nrow(cities.meta)), Val = m, stringsAsFactors = F)
    gc()
    print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(n, m, i, files, models, dlist, cities.batch)
```

#### Organize data frame and support objects


```r
d$Val[d$Var == "Temperature"] <- round(d$Val[d$Var == "Temperature"], 1)
d$Val[d$Var == "Precipitation"] <- round(d$Val[d$Var == "Precipitation"])
d$Month <- month.abb
d$Year <- results.years  #rep(1901:2009,each=12)
d$Month <- factor(d$Month, levels = d$Month[1:12])
cities.meta$loc <- paste0(cities.meta$loc, ", ", cities.meta$region)
cities.meta <- cities.meta[c(1, 3:6)]
names(cities.meta) <- c("Country", "Location", "Population", "Lat", "Lon")
d$Decade <- paste0(substr(d$Year, 1, 3), 0)
gc()
d.cities.cru31 <- d
rm(d, results.years)
gc()
# save.image('../final/data_cities_CRU31.RData')
```

#### Save output files


```r
library(parallel)
f <- function(i) {
    name.tmp <- gsub("`", "", gsub("~", "", gsub("?", "", gsub("\\'", "APOS", 
        cities.meta$Location[i]))))
    city.cru.dat <- subset(d.cities.cru31, Location == cities.meta$Location[i])
    save(city.cru.dat, file = paste0("../final/city_files_CRU/", gsub(", ", 
        "--", name.tmp), "__", domain, ".RData"))
    print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores = 32)

assign(paste0("cities.meta.", domain), cities.meta)

rm(d.cities.cru31, f, cities.meta, domain)
save(cities.meta.akcan2km, file = "../final/cities_meta_akcan2km.RData")
# load('../final/meta.RData') save.image('../final/meta.RData')
```

NEW_CODE_CHUNKS

```r
comArgs <- commandArgs(TRUE)
if (length(comArgs)) for (i in 1:length(comArgs)) eval(parse(text = comArgs[[i]]))
if (!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if (!exists("cities.batch")) cities.batch <- ""

setwd("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/cities")

library(data.table)

files <- list.files(pattern = paste0("^CRU31.*.batch", cities.batch, ".*.", 
    domain, ".RData$"))

models <- "CRU31"
```
