# CRU 3.x Samples Setup Code



## Introduction
The `samples_setup_CRU.R` script prepares extracted CRU 3.x regional temperature and precipitation data for direct use.

### Motivation
This script transforms and organizes output data from a parent script into a format that is easier to work with in certain frameworks.

### Details
The data are extracted based on regional shapefiles and the full Alaska-Canada 2-km extent for SNAP's CRU 3.x data set.
Data manipulation is performed to organize the inputs into more usefully formatted outputs.


#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
It uses as many CPU cores as are available since it is parallelized across regions, which exceed the number of cores.
It becomes serial beyond the number of cores simultaneously in use.

#### Limitations
This code expects CRU 3.x data have been processed by the parent script.
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.

## Related items

### Files and data
The input files are produced by the **R** script, `CRU_extract.R`. `samples_setup_CRU.R` assumes a complete and successful run of this precursory code.

The `CRU_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

For this reason, `samples_setup_CRU.R` has two companion **R** scripts which share the parent script, `CRU_extract.R`. These are `cities_setup_CRU.R` and `stats_setup_CRU.R`.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R`.
As a result, this script and its two companion scripts mentioned above are related to a similar set of three **R** scripts which perform similar intermediary data manipulation on extracted GCM data.

### Code flow

The collection of project code of which this script is a part exists within the context of a broader hierarchy of **R** code and data sets (shown below) spanning a multitude of other projects (not shown).
Similarly, outputs from this script and its companion scripts are used across a range of other projects.






## **R** Code

### Setup

Setup is minimal. Set working directory. List CRU 3.x files.


```r
setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/regional/samples")

library(parallel)
library(data.table)
library(reshape2)

cru <- "32"
files <- list.files(pattern = paste0("^CRU", cru, ".*.regions_samples.RData$"))

models <- paste0("CRU", cru)
```

#### Load and compile data


```r
dlist <- vector("list", length(files))
for (i in 1:length(files)) {
    load(files[i])
    m <- do.call(rbind, samples.out)
    n <- nrow(m)
    dlist[[i]] <- data.frame(Var = rep(c("Temperature", "Precipitation"), each = n/(2 * 
        length(samples.out))), m, stringsAsFactors = F)
    print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(samples.out, n, m, i, files, models, dlist)
names(d)[3:ncol(d)] <- gsub("X", "", names(d)[3:ncol(d)])
```

#### Organize data frame and support objects


```r
f <- function(i, n, index, cru, multiplier = c(1, 1), decimals = NULL) {
    name.tmp <- samples.names[i]
    rsd <- subset(d, Location == name.tmp)
    rsd <- melt(rsd, id.vars = names(rsd)[1:2], measure.vars = names(rsd)[-c(1:2)], 
        variable.name = "Time", value.name = "Vals_Probs")
    rsd$vals.ind <- rep(rep(c(T, F), each = n), length = nrow(rsd))
    rsd$dcastIsDumb <- rep(1:n, length = nrow(rsd))
    rsd <- dcast(rsd, Var + Location + Time + dcastIsDumb ~ vals.ind, value.var = "Vals_Probs")
    names(rsd)[6:5] <- c("Val", "Prob")
    rsd$Year <- substr(as.character(rsd$Time), 1, 4)
    rsd$Month <- substr(as.character(rsd$Time), 6, 8)
    rsd$Month <- factor(rsd$Month, levels = month.abb)
    rsd$Decade <- paste0(substr(rsd$Year, 1, 3), 0)
    rownames(rsd) <- NULL
    rsd <- rsd[c(1:2, 6, 5, 7:9)]
    grp <- rep(names(region.names.out), times = sapply(region.names.out, length))[i]
    dir.create(outDir <- paste0("../../final/region_files_CRU", cru, "/samples/", 
        grp, "/", name.tmp), recursive = T, showWarnings = F)
    if (!is.null(decimals)) 
        rsd$Val <- round(rsd$Val, decimals[1]) * multiplier[1]
    if (!is.null(decimals)) 
        rsd$Prob <- round(rsd$Prob, decimals[2]) * multiplier[2]
    if (i == 1) 
        x <- subset(rsd, Var == "Precipitation") else x <- NULL
    gc()
    rsd.cru <- unlist(subset(rsd, Var == "Precipitation")[, index])
    names(rsd.cru) <- NULL
    save(rsd.cru, file = paste0(outDir, "/", "precipitation.RData"))
    gc()
    rsd.cru <- unlist(subset(rsd, Var == "Temperature")[, index])
    names(rsd.cru) <- NULL
    save(rsd.cru, file = paste0(outDir, "/", "temperature.RData"))
    gc()
    print(i)
    x
}
```

#### Save output files


```r
if (file.exists("../../final/meta.RData")) load("../../final/meta.RData")
samples.columns.cru <- 3:4
samples.multipliers.cru <- c(10, 1e+16)  # see function arguments 'multiplier' and 'decimals'

out <- mclapply(1:length(samples.names), f, n = n.samples, index = samples.columns.cru, 
    cru = cru, multiplier = samples.multipliers.cru, decimals = c(1, 16), mc.cores = 32)
cru.samples.df <- out[[1]]
cru.samples.df[, samples.columns.cru] <- NA
cru.samples.df$Var <- NA
cru.samples.df$Location <- NA
if (as.numeric(cru) == 31) cru31.samples.df <- cru.samples.df else if (as.numeric(cru) == 
    32) cru32.samples.df <- cru.samples.df

rm(cru, cru.samples.df, d, f, out)
save.image("../../final/meta.RData")
```
