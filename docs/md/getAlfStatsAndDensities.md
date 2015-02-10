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

### Files and Data
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
library(parallel)
library(reshape2)
n.regions <- 53  # known, fixed
n.cores <- 27  # of 32 available
n.samples <- 20  # known, based on upstream settings for vegetation age
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", 
    "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
ageDirs <- list.files("/big_scratch/mfleonawicz/Rmpi/outputs/ageDensities", 
    full = T)
fireDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire"
vegDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/veg"
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

# include later: # rcp45='CMIP5', rcp60='CMIP5', rcp85='CMIP5'
getPhase <- function(x) {
    switch(x, sresb1 = "CMIP3", sresa1b = "CMIP3", sresa2 = "CMIP3")
}
```

Define support functions for density estimation, boostrap resampling from estimated densities, and extraction of mean, standard deviation, and quantile statistics.


```r
# Density estimation
denFun <- function(x, n = 20, min.zero = TRUE, diversify = FALSE, missing.veg.NA = TRUE, 
    fire = FALSE) {
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
    z <- density(x, adjust = 2, n = n, from = min(x) - max(1, 0.05 * dif), to = max(x) + 
        max(1, 0.05 * dif))
    if (min.zero && any(z$x < 0)) 
        z <- density(x, adjust = 2, n = n, from = 0, to = max(x) + max(1, 0.05 * 
            dif))
    as.numeric(c(z$x, z$y))
}

# Bootstrapping
btfun <- function(p, n.samples = length(p)/2, n.boot = 10000, interp = FALSE, 
    n.interp = 1000, ...) {
    if (!length(p)) 
        return(p)
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


```r
# Primary processing functions
get_AgeDensities <- function(j, dirs, n.samples = 20, n.samples.in, samples.index, 
    multiplier, veg.labels, scen.levels) {
    pat <- paste0("^veg.age.loc", j, ".rep.*.RData$")
    pat2 <- substr(pat, 1, nchar(pat) - 9)
    files.list <- lapply(1:length(dirs), function(i, dirs, ...) list.files(dirs[i], 
        ...), dirs = dirs, full = T, pattern = pat)
    for (i in 1:length(dirs)) {
        modname <- basename(dirs[i])
        # Load all reps for one location, all model and scenarios, into local
        # location-specific environment
        print(system.time(for (p in 1:length(files.list[[i]])) {
            load(files.list[[i]][p], envir = environment())
            print(paste("Local environment: loading workspace", p, "of", length(files.list[[i]])))
        }))
        d.names <- ls(pattern = pat2, envir = environment())
        d.list <- mget(d.names, envir = environment())
        print(paste("EVIRONMENT", paste0("e", j), "OBJECTS INCLUDE:", paste(ls(pattern = pat2, 
            envir = environment()), collapse = ", "), collapse = " "))
        rm(list = d.names, envir = environment())
        gc()
        all.years <- unique(d.list[[1]]$Year)
        loc.grp <- d.list[[1]]$LocGroup[1]
        loc <- d.list[[1]]$Location[1]
        nr <- sapply(d.list, nrow)
        year.veg.list <- lapply(d.list, function(x) as.numeric(unique(paste0(x$Year, 
            ".", x$VegID))))
        year.veg.uni <- sort(unique(unlist(year.veg.list)))
        a <- sapply(year.veg.list, function(x) all(paste0(all.years, ".", 1:8) %in% 
            x))
        dl <- length(d.list)
        for (z in 1:dl) {
            x <- as.matrix(d.list[[z]][, 3:5])
            x <- cbind(x[, 2], x[, 3] + x[, 1]/10)
            if (!a[z]) {
                x.ids <- unique(x[, 2])
                veg.missing <- paste0(rep(all.years, each = 8), ".", 1:8)
                veg.missing <- as.numeric(veg.missing[!(veg.missing %in% x.ids)])
                x <- rbind(x, cbind(NA, rep(veg.missing, each = 2 * n.samples)))
                x <- x[order(x[, 2]), ]
                rownames(x) <- NULL
            }
            d.list[[z]] <- split(x[, 1], x[, 2])
            print(paste("Obtain annual vegetation samples: replicate", z, "of", 
                dl))
        }
        d.list <- rapply(d.list, btfun, classes = "numeric", how = "replace", 
            n.samples = n.samples.in, n.boot = 100, interp = TRUE, n.interp = 1000)
        d.list <- unlist(d.list, recursive = FALSE)
        p <- length(d.list)
        m <- p/length(d.names)
        system.time(for (k in 1:m) {
            d.list[[k]] <- denFun(unlist(d.list[seq(1, p, by = m) + k - 1]), 
                n = n.samples)
            print(paste("Empirical density estimation: sample", k, "of", m))
        })
        d.list <- d.list[1:m]
        gc()
        d.tmp <- data.frame(Var = "Age", Location = loc, VegID = rep(1:8, each = 2 * 
            n.samples), Vals_Probs = unlist(d.list), Year = rep(all.years, each = 8 * 
            2 * n.samples), stringsAsFactors = FALSE)
        rm(d.list)
        gc()
        ind.val <- as.numeric(mapply("+", seq(1, nrow(d.tmp), by = 2 * n.samples), 
            MoreArgs = list(0:(n.samples - 1))))
        d.tmp$Vals_Probs[ind.val] <- round(d.tmp$Vals_Probs[ind.val])
        mod.scen <- unlist(strsplit(basename(dirs[i]), "\\."))
        d.tmp$Model <- swapModelName(mod.scen[1])
        d.tmp$Scenario <- swapScenarioName(mod.scen[2])
        d.tmp$Phase <- getPhase(mod.scen[2])
        d.tmp <- d.tmp[c(8:6, 1:5)]
        rownames(d.tmp) <- NULL
        if (i == 1) 
            d.alf.age <- d.tmp else d.alf.age <- rbind(d.alf.age, d.tmp)
        print(i)
    }
    
    rm(d.tmp)
    d.alf.age$vals.ind <- rep(rep(c(T, F), each = n.samples), length = nrow(d.alf.age))
    d.alf.age$dcastIsDumb <- rep(1:n.samples, length = nrow(d.alf.age))
    d.alf.age <- dcast(d.alf.age, Phase + Scenario + Model + Var + Location + 
        VegID + Year + dcastIsDumb ~ vals.ind, value.var = "Vals_Probs")
    names(d.alf.age)[10:9] <- c("Val", "Prob")
    d.alf.age$Scenario <- factor(d.alf.age$Scenario, levels = scen.levels)
    d.alf.age$Decade <- paste0(substr(d.alf.age$Year, 1, 3), 0)
    d.alf.age <- d.alf.age[order(paste(d.alf.age$Year, d.alf.age$VegID)), c(1:5, 
        10:9, 6:7, 11)]
    rownames(d.alf.age) <- NULL
    for (p in 1:length(veg.labels)) d.alf.age$VegID[d.alf.age$VegID == p] <- veg.labels[p]
    names(d.alf.age)[which(names(d.alf.age) == "VegID")] <- "Vegetation"
    # dir.create(statsDir <-
    # file.path('/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/region_files_GCM/stats',
    # loc.grp, loc), recursive=T, showWarnings=F)
    dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/region_files_GCM/samples", 
        loc.grp, loc), recursive = T, showWarnings = F)
    d.alf.age$Val <- round(d.alf.age$Val)
    d.alf.age$Prob <- round(d.alf.age$Prob, 8) * multiplier[2]  # round to eight seems decent for now
    r.alf.age <- unlist(d.alf.age[, samples.index])
    names(r.alf.age) <- NULL
    if (j == 1) 
        x <- d.alf.age else x <- NULL
    # if(j==1) y <- region.dat else y <- NULL
    save(r.alf.age, file = paste0(samplesDir, "/", "vegetationAge.RData"))
    gc()
    # region.dat <- unlist(region.dat[,stats.index]) names(region.dat) <- NULL
    # save(region.dat, file=file.path(statsDir, 'stats_age.RData'))
    return(list(alf.ageSamples.df = x))  #, alf.fireStats.df=y))
}
```

#### get_FireStats_FireDensities


```r
get_FireStats_FireDensities <- function(j, inDir, n.samples = 20, stats.index, 
    samples.index, multiplier, scen.levels) {
    pat <- paste0("^abfc.loc", j, "\\..*.RData$")
    files <- list.files(inDir, full = T, pattern = pat)
    for (i in 1:length(files)) {
        load(files[i], envir = environment())
        d.tmp <- f.dat
        rm(f.dat)
        gc()
        loc.grp <- d.tmp$LocGroup[1]
        loc <- d.tmp$Location[1]
        ind <- with(d.tmp, paste(Var, Year))
        v <- unlist(tapply(d.tmp$Val, ind, denFun, n = n.samples, fire = TRUE))
        d.tmp <- d.tmp[d.tmp$Replicate == 1, c(1:3, 6, 5, 7, 8)]
        stats.tmp <- d.tmp
        d.tmp <- data.frame(lapply(d.tmp, rep, n.samples), stringsAsFactors = FALSE)
        d.tmp$Val <- round(v[rep(c(T, F), each = n.samples)])
        d.tmp$Prob <- v[rep(c(F, T), each = n.samples)]
        d.tmp <- d.tmp[c(1:6, 8, 7)]
        v <- tapply(v, rep(1:(length(v)/(2 * n.samples)), each = 2 * n.samples), 
            btfun, n.samples = n.samples, n.boot = 1000, interp = TRUE)
        qtiles <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
        n.stats <- 2 + length(qtiles)
        v <- do.call(rbind, lapply(v, getStats, seq.q = qtiles))
        stats.tmp <- data.frame(stats.tmp, v)
        stats.tmp <- stats.tmp[c(1:5, 7 + 1:ncol(v), 7)]
        if (i == 1) 
            d.alf.fire <- d.tmp else d.alf.fire <- rbind(d.alf.fire, d.tmp)
        if (i == 1) 
            region.dat <- stats.tmp else region.dat <- rbind(region.dat, stats.tmp)
    }
    rm(d.tmp, stats.tmp)
    gc()
    d.alf.fire$Scenario <- factor(d.alf.fire$Scenario, levels = scen.levels)
    d.alf.fire$Decade <- paste0(substr(d.alf.fire$Year, 1, 3), 0)
    region.dat$Scenario <- factor(region.dat$Scenario, levels = scen.levels)
    region.dat$Decade <- paste0(substr(region.dat$Year, 1, 3), 0)
    dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/region_files_GCM/stats", 
        loc.grp, loc), recursive = T, showWarnings = F)
    dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/region_files_GCM/samples", 
        loc.grp, loc), recursive = T, showWarnings = F)
    d.alf.fire$Prob <- round(d.alf.fire$Prob, 8) * multiplier[2]  # round to eight seems decent for now
    r.alf.fire <- unlist(subset(d.alf.fire, Var == "Burn Area")[, samples.index])
    names(r.alf.fire) <- NULL
    if (j == 1) 
        x <- subset(d.alf.fire, Var == "Burn Area") else x <- NULL
    # if(j==1) y <- subset(region.dat, Var=='Burn Area') else y <- NULL if(j==1)
    # x <- d.alf.fire else x <- NULL
    if (j == 1) 
        y <- region.dat else y <- NULL
    save(r.alf.fire, file = paste0(samplesDir, "/", "burnArea.RData"))
    gc()
    r.alf.fire <- unlist(subset(d.alf.fire, Var == "Fire Count")[, samples.index])
    names(r.alf.fire) <- NULL
    save(r.alf.fire, file = file.path(samplesDir, "fireCount.RData"))
    gc()
    region.dat <- unlist(region.dat[, stats.index])
    names(region.dat) <- NULL
    save(region.dat, file = file.path(statsDir, "stats_fire.RData"))
    return(list(alf.fireSamples.df = x, alf.fireStats.df = y))
}
```

#### get_VegStats_VegDensities


```r
get_VegStats_VegDensities <- function(j, inDir, n.samples = 20, stats.index, 
    samples.index, multiplier, veg.labels, scen.levels) {
    pat <- paste0("^v.loc", j, "\\..*.RData$")
    files <- list.files(inDir, full = T, pattern = pat)
    for (i in 1:length(files)) {
        load(files[i], envir = environment())
        # locs <- unique(v.dat$Location) d.tmp <- v.dat[v.dat$Location==locs[j],]
        d.tmp <- data.frame(v.dat[1, 1:6], Val = 0, VegID = 1:8, Year = rep(unique(v.dat$Year), 
            each = 8), Replicate = rep(unique(v.dat$Replicate), each = 8 * length(unique(v.dat$Year))))
        obs.ind <- paste(d.tmp$Year, d.tmp$VegID, d.tmp$Replicate) %in% paste(v.dat$Year, 
            v.dat$VegID, v.dat$Replicate)
        d.tmp$Val[obs.ind] <- v.dat$Val
        rm(v.dat)
        gc()
        loc.grp <- d.tmp$LocGroup[1]
        loc <- d.tmp$Location[1]
        ind <- with(d.tmp, paste(Year, VegID))
        v <- unlist(tapply(d.tmp$Val, ind, denFun, n = n.samples))
        d.tmp <- d.tmp[d.tmp$Replicate == 1, c(1:3, 6, 5, 7, 8, 9)]
        stats.tmp <- d.tmp
        d.tmp <- data.frame(lapply(d.tmp, rep, n.samples), stringsAsFactors = FALSE)
        d.tmp <- d.tmp[order(paste(d.tmp$Year, d.tmp$VegID)), ]
        for (p in 1:length(veg.labels)) stats.tmp$VegID[stats.tmp$VegID == p] <- veg.labels[p]
        for (p in 1:length(veg.labels)) d.tmp$VegID[d.tmp$VegID == p] <- veg.labels[p]
        names(stats.tmp)[7] <- names(d.tmp)[7] <- "Vegetation"
        d.tmp$Val <- round(v[rep(c(T, F), each = n.samples)])
        d.tmp$Prob <- v[rep(c(F, T), each = n.samples)]
        d.tmp <- d.tmp[c(1:6, 9, 7, 8)]
        v <- tapply(v, rep(1:(length(v)/(2 * n.samples)), each = 2 * n.samples), 
            btfun, n.samples = n.samples, n.boot = 1000, interp = TRUE)
        qtiles <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
        n.stats <- 2 + length(qtiles)
        v <- do.call(rbind, lapply(v, getStats, seq.q = qtiles))
        stats.tmp <- data.frame(stats.tmp, v)
        stats.tmp <- stats.tmp[c(1:5, 8 + 1:ncol(v), 7, 8)]
        if (i == 1) 
            d.alf.veg <- d.tmp else d.alf.veg <- rbind(d.alf.veg, d.tmp)
        if (i == 1) 
            region.dat <- stats.tmp else region.dat <- rbind(region.dat, stats.tmp)
    }
    rm(d.tmp, stats.tmp)
    gc()
    d.alf.veg$Scenario <- factor(d.alf.veg$Scenario, levels = scen.levels)
    d.alf.veg$Decade <- paste0(substr(d.alf.veg$Year, 1, 3), 0)
    region.dat$Scenario <- factor(region.dat$Scenario, levels = scen.levels)
    region.dat$Decade <- paste0(substr(region.dat$Year, 1, 3), 0)
    dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/region_files_GCM/stats", 
        loc.grp, loc), recursive = T, showWarnings = F)
    dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/region_files_GCM/samples", 
        loc.grp, loc), recursive = T, showWarnings = F)
    d.alf.veg$Prob <- round(d.alf.veg$Prob, 8) * multiplier[2]  # round to eight seems decent for now
    r.alf.veg <- unlist(d.alf.veg[, samples.index])
    names(r.alf.veg) <- NULL
    if (j == 1) 
        x <- d.alf.veg else x <- NULL
    if (j == 1) 
        y <- region.dat else y <- NULL
    save(r.alf.veg, file = file.path(samplesDir, "vegetatedArea.RData"))
    gc()
    region.dat <- unlist(region.dat[, stats.index])
    names(region.dat) <- NULL
    save(region.dat, file = file.path(statsDir, "stats_veg.RData"))
    return(list(alf.vegSamples.df = x, alf.vegStats.df = y))
}
```

### Processing
Prepare for processing.


```r
# Processing
agg.stat.colnames <- c("Mean", "SD", paste0("Pct_", c("05", 10, 25, 50, 75, 
    90, 95)))
stats.columns <- 5 + 1:length(agg.stat.colnames)
samples.columns <- 6:7
samples.multipliers.alf <- c(10, 1e+08)  # First not used
```

Process fire statistics and distributions.


```r
system.time(out.fire <- mclapply(1:n.regions, get_FireStats_FireDensities, inDir = fireDir, 
    n.samples = n.samples, stats.index = stats.columns, samples.index = samples.columns, 
    multiplier = samples.multipliers.alf, mc.cores = n.cores))
alf.fireStats.df <- out.fire[[1]]$alf.fireStats.df
alf.fireSamples.df <- out.fire[[1]]$alf.fireSamples.df
rm(out.fire)
gc()

alf.fireStats.df[, stats.columns] <- NA
alf.fireStats.df$Location <- NA
names(alf.fireStats.df)[stats.columns] <- agg.stat.colnames

alf.fireSamples.df[, samples.columns] <- NA
alf.fireSamples.df$Var <- NA
alf.fireSamples.df$Location <- NA
```

Process vegetation class statistics and distributions.


```r
system.time(out.veg <- mclapply(1:n.regions, get_VegStats_VegDensities, inDir = vegDir, 
    n.samples = n.samples, stats.index = stats.columns, samples.index = samples.columns, 
    multiplier = samples.multipliers.alf, veg.labels = veg.labels, mc.cores = n.cores))
alf.vegStats.df <- out.veg[[1]]$alf.vegStats.df
alf.vegSamples.df <- out.veg[[1]]$alf.vegSamples.df
rm(out.veg)
gc()

alf.vegStats.df[, stats.columns] <- NA
alf.vegStats.df$Location <- NA
names(alf.vegStats.df)[stats.columns] <- agg.stat.colnames

alf.vegSamples.df[, samples.columns] <- NA
alf.vegSamples.df$Var <- NA
alf.vegSamples.df$Location <- NA
```

Process vegetation age distributions.


```r
system.time(out.age <- mclapply(1:n.regions, get_AgeDensities, dirs = ageDirs, 
    n.samples = n.samples, n.samples.in = 20, samples.index = samples.columns, 
    multiplier = samples.multipliers.alf, veg.labels = veg.labels, mc.cores = n.cores))
# alf.ageStats.df <- out.age[[1]]$alf.ageStats.df
alf.ageSamples.df <- out.age[[1]]$alf.ageSamples.df
rm(out.age)
gc()

# alf.ageStats.df[,stats.columns] <- NA alf.ageStats.df$Location <- NA
# names(alf.ageStats.df)[stats.columns] <- agg.stat.colnames

alf.ageSamples.df[, samples.columns] <- NA
alf.ageSamples.df$Var <- NA
alf.ageSamples.df$Location <- NA
```

Save metadata into `meta.RData`, which is used by the master QAQC Shiny app.
Currently, this is saved in a backup file as a specific version of metadata since the integration of ALFRESCO output statistics into the app is in an early developmental stage and not used by most repo branches.


```r
# Remove unwanted objects, load metadata workspace, save along with age
# metadata
rm(ageDirs, agg.stat.colnames, btfun, denFun, fireDir, get_AgeDensities, get_FireStats_FireDensities, 
    get_VegStats_VegDensities, getPhase, getSamples, getStats, lapply_C, n.cores, 
    n.regions, n.samples, samples.columns, scen.levels, stats.columns, swapModelName, 
    swapScenarioName, vegDir)
# Create a backup of the meta.RData file, load and save over that.
# Inclusion of ALFRESCO output in QAQC R Shiny app is under early
# development.  All but one shiny-apps development repo branch should not
# include this version of meta.RData in the cmip3_cmip5 app
load("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/meta_backup.RData")
save.image("/workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/data/final/meta_backup.RData")
```
