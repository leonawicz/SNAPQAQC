# ALFRESCO Statistics and Densities



## Introduction

The `getAlfStatsAndDensities.R` script loads, combines, and organizes extracted ALFRESCO output statistics and distributional information from intermediary .RData workspace files produced by upstream **R** scripts in the processing chain.
Curated outputs are saved to **R** workspaces for analysis and graphing by subsequent **R** code.

### Motivation
The primary motivation for this code is not just to extract regional and point data from a large number of high-resolution geotiffs,
but to limit the routine recurrence and redundancy of such extractions.
It allows for storing commonly required data in a more compact format that can be quickly shared and digested by other projects.

### Details

#### Capabilities
This script uses parallel processing via the base `parallel` package and the function `mclapply`.
Parallel processing is across regions.
Processing is serial across model-scenario pairs.
`getAlfStatsAndDensities.R` represents the stage in the processing chain at which the 200 simulation replicates from each ALFRESCO run are finally aggregated together, greatly reducing the overall amount of data carry-through,
despite that additional distributional estimation leads to a smaller uptick in total information.

#### Limitations
The code does not make use of `Rmpi` or similar options so it cannot take advantage of multi-node cluster parallel processing.
`getAlfStatsAndDensities.R` is not currently called via slurm script, but should be.

### Files and data
Input files include .RData workspace files storing fire, vegetation, and age data unique to each model, scenario, simulation replicate and spatial region.

The `getAlfStatsAndDensities.R` script produces the following curated outputs:
* Spatially aggregated summary statistics of burn area and fire frequency by region
* Spatially aggregated summary statistics of vegetation cover by vegetation class and region
* Regional probability distribution estimates for these same regions for burn area, fire frequency, and vegetation cover, and vegetation age.

Things to note:
* Age probability distributions have a joint support involving the spatially explicit distribution of ages on the landscape and the variation across the simulation replicates.
* As spatially aggregated statistics by definition, burn area, fire frequency, and vegetation cover probability distributions have a univariate support based strictly on variation across simulation replicates.
* Outputs, particularly the distributional information, are stored compactly, e.g., as a numeric vector in an **R** workspace, while a template data frame is stored in a metadata workspace used across various projects.

## R code

### Setup
Setup consists of loading required **R** packages and defining file paths and other **R** objects including scenario names, vegetation classes, numbers of processing cores, and sampling factor for distributional estimation.


```r
comargs <- (commandArgs(TRUE))
if (!length(comargs)) q("no") else for (z in 1:length(comargs)) eval(parse(text = comargs[[z]]))

library(parallel)
library(reshape2)
library(data.table)
library(dplyr)

if (!exists("mainDir")) mainDir <- "/big_scratch/mfleonawicz/Rmpi/outputs"
if (!exists("variable")) stop("Must provide 'variable' argument in escaped quotes. Options are 'age' (vegetation age), 'veg' (vegetated area), 'abfc' (area burned and fire counts), or 'fsv' (fire sizes by vegetation class).")
ageDirs <- list.files(file.path(mainDir, "vegDir"), full = T)
empty <- which(lapply(ageDirs, function(x) length(list.files(x))) == 0)
if (length(empty)) ageDirs <- ageDirs[-empty]
fsvDir <- file.path(mainDir, "fsv")

n.samples <- 1000  # new density estimation
n.samples.in <- 1000  # known, based on upstream settings for vegetation age
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", 
    "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
datnames <- c("Phase", "Scenario", "Model", "Location", "Var", "Vegetation", 
    "Year", "Val")

# All regions for which stage-1 outputs currently exist on disk.  Checking
# abfc files, which have maximum regions, but only regions common to fsv,
# abfc, veg, and age outputs are processed.  It is assumed stage one
# processing has been done on a common set of regions for all four variable
# sets.
regions <- unique(sapply(strsplit(list.files(file.path(fsvDir)), "__"), "[", 
    2))
n.regions <- length(regions)
n.cores <- min(n.regions, 32)
```

### Support functions
Define support functions which assist in assembling curated data frames.


```r
# Support functions include later: # CCSM4='CCSM4', GFDLcm3='GFDLcm3',
# GISSe2-r='GISSe2-r', IPSLcm5a-lr='IPSLcm5a-lr', MRIcgcm3='MRIcgcm3'
swapModelName <- function(x) {
    switch(x, cccma_cgcm3_1 = "CCCMAcgcm31", gfdl_cm2_1 = "GFDLcm21", miroc3_2_medres = "MIROC32m", 
        mpi_echam5 = "MPIecham5", ukmo_hadcm3 = "ukmoHADcm3")
}

# include later: # rcp45='RCP 4.5', rcp60='RCP 6.0', rcp85='RCP 8.5'
swapScenarioName <- function(x) {
    switch(x, sresb1 = "SRES B1", sresa1b = "SRES A1B", sresa2 = "SRES A2")
}

# include later: # rcp45='AR5', rcp60='AR5', rcp85='AR5'
getPhase <- function(x) {
    switch(x, sresb1 = "AR4", sresa1b = "AR4", sresa2 = "AR4")
}
```

Define support functions for density estimation, boostrap resampling from estimated densities, and extraction of mean, standard deviation, and quantile statistics.


```r
# Density estimation
denFun <- function(x, n = 1000, adj = 0.1, min.zero = TRUE, diversify = FALSE, 
    missing.veg.NA = TRUE, fire = FALSE) {
    if (all(is.na(x))) 
        return(rep(NA, 2 * n))
    x <- x[!is.na(x)]
    lx <- length(x)
    if (sum(x) == 0 & missing.veg.NA & !fire) 
        return(rep(NA, 2 * n))
    if (diversify && length(unique(x)) == 1) 
        x <- rnorm(max(10, lx), mean = x[1])  # diversify constant values
    if (lx == 1) 
        x <- x + c(-1:1)  #single pixel of veg type, add and subtract one age year to make procedure possible
    dif <- diff(range(x))
    z <- density(x, adjust = adj, n = n, from = min(x) - max(1, 0.05 * dif), 
        to = max(x) + max(1, 0.05 * dif))
    if (min.zero && any(z$x < 0)) 
        z <- density(x, adjust = adj, n = n, from = 0, to = max(x) + max(1, 
            0.05 * dif))
    as.numeric(c(z$x, z$y))
}

# Bootstrapping
btfun_dt <- function(p1, p2, n.samples = length(p)/2, n.boot = 10000, interp = FALSE, 
    n.interp = 1000, ...) {
    p <- c(p1, p2)  # concatenate separate columns from data table
    if (!length(p)) 
        return(NULL)
    if (all(is.na(p))) 
        return(rep(NA, n.boot))
    p <- list(x = p[1:n.samples], y = p[(n.samples + 1):(2 * n.samples)])
    if (interp && length(unique(p[1:n.samples])) > 1) 
        p <- approx(p$x, p$y, n = n.interp)
    p <- round(sample(p$x, n.boot, prob = p$y, rep = T), ...)
    p
}
```

### Processing functions
#### get_AgeDensities



#### get_FireStats_FireDensities



#### get_VegStats_VegDensities



### Processing
Prepare for processing.



Process fire statistics and distributions.


```r
if ("abfc" %in% variable) out.fire <- mclapply(1:n.regions, get_fireStats_fireDensities, 
    inDir = abfcDir, n.samples = n.samples, datnames = datnames, scen.levels = scen.levels, 
    mc.cores = n.cores)
```

Process vegetation class statistics and distributions.


```r
if ("veg" %in% variable) out.veg <- mclapply(1:n.regions, get_vegStats_vegDensities, 
    inDir = vegDir, n.samples = n.samples, datnames = datnames, veg.labels = veg.labels, 
    scen.levels = scen.levels, mc.cores = n.cores)
```

Process vegetation age distributions.


```r
if ("age" %in% variable) out.age <- mclapply(1:n.regions, get_ageStats_ageDensities, 
    dirs = ageDirs, n.samples = n.samples, n.samples.in = n.samples.in, datnames = datnames, 
    veg.labels = veg.labels, scen.levels = scen.levels, mc.cores = n.cores)
```

Save metadata into `meta.RData`, which is used by the master QA/QC Shiny app.
Currently, this is saved in a backup file as a specific version of metadata since the integration of ALFRESCO output statistics into the app is in an early developmental stage and not used by most repo branches.


