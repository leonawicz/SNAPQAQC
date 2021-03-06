# ALFRESCO Data Prep



## Introduction

The `alfPrep.R` script loads, combines, and organizes extracted ALFRESCO random variable distributional information from intermediary .RData workspace files produced by upstream **R** scripts in the processing chain.
Curated outputs are saved to **R** workspaces for analysis and graphing by subsequent **R** code.

### Motivation
The primary motivation for this code is not just to extract regional and point data from a large number of high-resolution geotiffs,
but to limit the routine recurrence and redundancy of such extractions.
It allows for storing commonly required data in a more compact format that can be quickly shared and digested by other projects.

### Details

This script uses parallel processing via the base `parallel` package and the function `mclapply`.
Parallel processing is across regions.
Processing is serial across model-scenario pairs.
The code does not make use of `Rmpi` or similar options so it does not take advantage of multi-node cluster parallel processing.
`alfPrep.R` is called via slurm script, `alfPrep.slurm`.

### Files and data
Input files include .RData workspace files storing fire, vegetation, and age data unique to each model, scenario, simulation replicate and spatial region.

The `alfPrep.R` script produces the following curated outputs. Spatially aggregated summary statistics and modeled distributions by region and vegetation class for:

* burn area
* fire frequency
* fire size
* vegetation age
* vegetated cover area

Things to note:
* Age probability distributions have a joint support involving the spatially explicit distribution of ages on the landscape and the variation across the simulation replicates.
* As spatially aggregated statistics by definition, burn area, fire frequency, and vegetation cover probability distributions have a univariate support based strictly on variation across simulation replicates.

## R code

### Setup


```r
comargs <- (commandArgs(TRUE))
if (!length(comargs)) q("no") else for (z in 1:length(comargs)) eval(parse(text = comargs[[z]]))

library(parallel)
library(reshape2)
library(data.table)
library(dplyr)

if (!exists("mainDir")) mainDir <- "/big_scratch/mfleonawicz/Rmpi/outputs"
if (!exists("variable")) stop("Must provide 'variable' argument in escaped quotes. Options are 'age' (vegetation age), 'veg' (vegetated area), or 'fsv' (fire sizes by vegetation class).")
stopifnot(length(variable) == 1 && variable %in% c("age", "veg", "fsv"))
inDir <- file.path(mainDir, variable)
outDir <- "/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final/alfresco"

# All regions for which stage-1 outputs currently exist on disk.  Checking
# fsv files, but it is assumed stage one processing has been done on a
# common set of regions for all variables.
regions <- unique(sapply(strsplit(list.files(file.path(inDir)), "__"), "[", 
    2))
n.regions <- length(regions)
n.cores <- min(n.regions, 32)
```

### Support functions
Define support functions which assist in assembling curated data frames.


```r
# Support functions Density estimation
dtDen <- function(x, n = 1000, adj = 0.1, out = "vector", min.zero = TRUE, diversify = FALSE) {
    x <- x[!is.na(x)]
    lx <- length(x)
    if (diversify && length(unique(x)) == 1) 
        x <- rnorm(max(10, lx), mean = x[1])  # diversify constant values
    if (lx == 1) 
        x <- x + c(-1:1)  #single pixel of veg type, add and subtract one age year to make procedure possible
    b <- max(1, 0.05 * diff(range(x)))
    z <- density(x, adjust = adj, n = n, from = min(x) - b, to = max(x) + b)
    if (min.zero && any(z$x < 0)) 
        z <- density(x, adjust = adj, n = n, from = 0, to = max(x) + b)
    if (out == "vector") 
        return(as.numeric(c(z$x, z$y))) else if (out == "list") 
        return(z)
}

# Bootstrapping
dtBoot <- function(p, p2 = NULL, n.boot = 10000, interp = TRUE, n.interp = 1e+05, 
    round.samples = FALSE) {
    stopifnot(is.logical(round.samples) || is.na(as.integer(round.samples)))
    lp <- length(p)
    p <- if (is.null(p2)) 
        list(x = p[1:(lp/2)], y = p[(lp/2 + 1):lp]) else list(x = p, y = p2)
    if (interp) 
        p <- approx(p$x, p$y, n = n.interp)
    p <- sample(p$x, n.boot, prob = p$y, rep = T)
    if (round.samples == FALSE) 
        return(p) else if (round.samples == TRUE) 
        return(round(p)) else return(p, round.samples)
}
```

### Processing function
#### prep_data


```r
# Primary processing functions
prep_data <- function(j, inDir, outDir, n.samples = 1000, n.boot = 10000, ...) {
    id <- basename(inDir)
    exact <- list(...)$exact
    if (is.null(exact) || !is.logical(exact) || id != "veg") 
        exact <- FALSE
    files <- list.files(inDir, full = T, pattern = paste0("^", id, "__.*.RData$"))
    files.locs <- sapply(strsplit(files, "__"), "[", 2)
    locs <- unique(files.locs)
    if (j > length(locs)) 
        return(NULL)
    loc <- locs[j]
    files <- files[which(files.locs %in% loc)]
    dat <- stat <- vector("list", length(files))
    for (i in 1:length(files)) {
        load(files[i], envir = environment())
        d <- get(ls(pattern = "^d\\."))
        rm(list = ls(pattern = "^d\\."))
        loc.grp <- d$LocGroup[1]
        loc <- d$Location[1]
        if (id == "fsv") {
            d2 <- group_by(d, Phase, Scenario, Model, Location, Var, Year, Replicate, 
                FID) %>% summarise(Val = sum(Val)) %>% mutate(Vegetation = "All")  # agg-veg FS
            d <- data.table(bind_rows(d, d2)) %>% mutate(Vegetation = factor(Vegetation, 
                levels = unique(Vegetation))) %>% group_by(Phase, Scenario, 
                Model, Location, Var, Vegetation, Year)  # individual and aggregate-veg fire sizes
            d2 <- group_by(d, Replicate, add = T) %>% summarise(BA = sum(Val), 
                FC = length(Val))  # burn area and fire frequency
            d <- summarise(d, Val = dtDen(Val, n = n.samples, out = "list")$x, 
                Prob = dtDen(Val, n = n.samples, out = "list")$y)
            d2.ba <- summarise(d2, Val = dtDen(BA, n = n.samples, out = "list")$x, 
                Prob = dtDen(BA, n = n.samples, out = "list")$y) %>% mutate(Var = "Burn Area")
            d2.fc <- summarise(d2, Val = dtDen(FC, n = n.samples, out = "list")$x, 
                Prob = dtDen(FC, n = n.samples, out = "list")$y) %>% mutate(Var = "Fire Count")
            rm(d2)
            d <- data.table(bind_rows(d, d2.ba, d2.fc))
            rm(d2.ba, d2.fc)
            gc()
        }
        d <- group_by(d, Phase, Scenario, Model, Location, Var, Vegetation, 
            Year)
        if (id == "age") 
            d <- summarise(d, Val = dtDen(sample(Age, n.boot, T, Freq), n = n.samples, 
                out = "list")$x, Prob = dtDen(sample(Age, n.boot, T, Freq), 
                n = n.samples, out = "list")$y) %>% group_by(Year, add = T)
        if (!exact & id == "veg") 
            d <- summarise(d, Val = dtDen(Val, n = n.samples, out = "list")$x, 
                Prob = dtDen(Val, n = n.samples, out = "list")$y) %>% group_by(Year, 
                add = T)
        get_stats <- function(data, exact = FALSE) {
            if (!exact) 
                data <- summarise(data, Val = dtBoot(Val, Prob, n.boot = n.boot)) %>% 
                  group_by(Year, add = T)
            summarise(data, Mean = round(mean(Val)), SD = round(sd(Val), 1), 
                Min = round(min(Val)), Pct_05 = round(quantile(Val, 0.05)), 
                Pct_10 = round(quantile(Val, 0.1)), Pct_25 = round(quantile(Val, 
                  0.25)), Pct_50 = round(quantile(Val, 0.5)), Pct_75 = round(quantile(Val, 
                  0.75)), Pct_90 = round(quantile(Val, 0.9)), Pct_95 = round(quantile(Val, 
                  0.95)), Max = round(max(Val))) %>% group_by(Year, add = T)
        }
        s <- get_stats(d, exact = exact)
        dat[[i]] <- d
        stat[[i]] <- s
        print(i)
    }
    dir.create(statsDir <- file.path(outDir, "stats", loc.grp, loc), recursive = T, 
        showWarnings = F)
    dir.create(samplesDir <- file.path(outDir, "samples", loc.grp, loc), recursive = T, 
        showWarnings = F)
    if (id == "fsv") {
        dat <- rbindlist(dat)
        d.alf.fs <- filter(dat, Var == "Fire Size")
        d.alf.ba <- filter(dat, Var == "Burn Area")
        d.alf.fc <- filter(dat, Var == "Fire Count")
        stats.alf.fire <- rbindlist(stat)
        save(d.alf.fs, file = file.path(samplesDir, "fsByVeg.RData"))
        save(d.alf.ba, file = file.path(samplesDir, "baByVeg.RData"))
        save(d.alf.fc, file = file.path(samplesDir, "fcByVeg.RData"))
        save(stats.alf.fire, file = file.path(statsDir, "stats_fsbafcByVeg.RData"))
    } else {
        data.obj.name <- switch(id, age = "d.alf.vegage", veg = "d.alf.vegarea")
        stats.obj.name <- switch(id, age = "stats.alf.vegage", veg = "stats.alf.vegarea")
        assign(data.obj.name, rbindlist(dat))
        assign(stats.obj.name, rbindlist(stat))
        filename <- paste0(tail(strsplit(data.obj.name, "\\.")[[1]], 1), ".RData")
        save(list = data.obj.name, file = file.path(samplesDir, filename))
        save(list = stats.obj.name, file = file.path(statsDir, filename))
    }
    return()
}
```

### Processing


```r
mclapply(1:n.regions, prep_data, inDir, outDir, mc.cores = n.cores)
```
