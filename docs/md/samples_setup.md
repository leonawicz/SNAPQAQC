# GCM Samples Setup Code



## Introduction
The `samples_setup.R` script prepares extracted AR4/CMIP3 and AR5/CMIP5 regional temperature and precipitation data for direct use.

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

### Files and data
The input files are produced by the **R** script, `AR4_AR5_extract.R`. `samples_setup.R` assumes a complete and successful run of this precursory code.

The `AR4_AR5_extract.R` script produces three main sets of output data:
* Data for variables at the scale and location of individual communities
* Spatially aggregated summary statistics of this data over larger regions
* Regional probability distribution estimates for these same regions.

For this reason, `samples_setup.R` has two companion **R** scripts which share the parent script, `AR4_AR5_extract.R`. These are `cities_setup.R` and `stats_setup.R`.

In the hierarchy of these extraction and data preparation scripts, `AR4_AR5_extract.R` exists alongside `CRU_extract.R`.
As a result, this script and its two companion scripts mentioned above are related to a similar set of three **R** scripts which perform similar intermediary data manipulation on extracted CRU 3.x data.

## **R** Code

### Setup

Setup is minimal. Set working directory. List GCM files while making sure to avoid CRU 3.x files.
SNAP's ten combined CMIP3 and CMIP5 models are in a hardcoded list as in the parent script.


```r
setwd("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/regional/samples")

library(parallel)
library(reshape2)
library(data.table)

files <- list.files(pattern = "regions_samples.RData$")
files <- files[substr(files, 1, 3) != "CRU"]
lab <- do.call(rbind, strsplit(files, "_"))[, 1:4]
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")

swapVar <- function(x) switch(x, pr = "Precipitation", tas = "Temperature")
swapScen <- function(x) switch(x, sresb1 = "SRES B1", sresa1b = "SRES A1B", 
    sresa2 = "SRES A2", rcp45 = "RCP 4.5", rcp60 = "RCP 6.0", rcp85 = "RCP 8.5")
swapMod <- function(x) {
    switch(x, `cccma-cgcm3-1-t47` = "CCCMAcgcm31", `gfdl-cm2-1` = "GFDLcm21", 
        `miroc3-2-medres` = "MIROC32m", `mpi-echam5` = "MPIecham5", `ukmo-hadcm3` = "ukmoHADcm3", 
        CCSM4 = "CCSM4", `GFDL-CM3` = "GFDLcm3", `GISS-E2-R` = "GISSe2-r", `IPSL-CM5A-LR` = "IPSLcm5a-lr", 
        `MRI-CGCM3` = "MRIcgcm3")
}

lab[, 2:4] <- cbind(sapply(lab[, 2], swapScen), sapply(lab[, 3], swapMod), sapply(lab[, 
    4], swapVar))
rm(swapVar, swapScen, swapMod)
```

#### Load and compile data



#### Organize data frame and support objects



#### Save output files


