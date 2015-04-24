# ALFRESCO Statistics: Rmpi Code



## Introduction

The `alfStatsByRep_Rmpi.R` script uses `Rmpi` to efficiently extract fire and vegetation statistics from calibrated, finalized ALFRESCO run output over the full Alaska-Canada 1-km resolution ALFRESCO extent.
Data are saved to **R** workspaces (.RData files) for analysis and graphing by subsequent **R** code.

### Motivation
The primary motivation for this code is not just to extract regional data from a large number of high-resolution geotiffs,
but to limit the routine recurrence and redundancy of such extractions.
It allows for storing commonly required data in a more compact format that can be quickly shared and digested by other projects.

### Details

#### Capabilities
Extractions occur across multiple climate model and scenario-driven runs for multiple simulation replicates.
Statistics are extracted from geotiff output map layers for the entire spatial extent as well as an arbitrary number of spatial subregions.

#### Limitations
This script is intended to run on a single climate model-driven set of ALFRESCO output files at a time; one model-scenario pair.
This is because the outputs for each model-scenario pair consist of 200 simulation replicates.
This is plenty to process at one time.
Scripts further downstream are written with some expectation of this upstream process compiling output statistics based on all replicates,
regardless of how many models or scenarios are processed at once.

### Files and data
The input files include 1-km Alaska-Canada extent data from ALFRESCO simulations
as well as raster layer cell indices pertaining to various regional shapefile polygons.
`alfStatsByRep_Rpmi.R` is called via slurm script, `alfStatsByRep_Rmpi.slurm` and is itself a wrapper around `alfStatsByRep.R` which performs data extractions on multiple CPU cores across multiple compute nodes.

The `alfStatsByRep_Rpmi.R` gathers on a head node the data extracted from ALFRESCO outputs on the CPU cores across all the processing nodes and writes these data to intermediary .RData workspace files:
* Burn area by region
* Fire frequency by region
* Vegetated area by region and vegetation class

Other data are written to files directly by the compute nodes, avoiding returning too much data to the head node.
This consists of distributional information, in contrast to specific computed statistics.

Intermediate outputs are handled by additional **R** scripts.

## R code

### Setup
ADD_TEXT_HERE: EXAMPLE
Setup consists of loading required **R** packages and additional files, preparing any command line arguments for use, and defining functions and other **R** objects.


```r
# Command line arguments
args = (commandArgs(TRUE))
if (!length(args)) q("no") else for (i in 1:length(args)) eval(parse(text = args[[i]]))

if (!exists("model.index")) stop("Must provide model(s) as integer(s) 1:15, e.g., model.index=1:2")
if (!exists("reps")) stop("Must provide replicates as integer(s) 1:200, e.g., reps=1:25")
if (!exists("years")) years <- NULL  # Assume all years if none specified
if (!exists("Rmpi")) Rmpi <- TRUE
if (!exists("doFire")) doFire <- TRUE
if (!exists("doAgeVeg")) doAgeVeg <- TRUE
if (exists("repSample") && is.numeric(repSample)) {
    set.seed(47)
    reps <- sort(sample(reps, min(repSample, length(reps))))
}

library(data.table)

# Rmpi setup
if (Rmpi) {
    library(Rmpi)
    mpi.spawn.Rslaves(needlog = TRUE)
    mpi.bcast.cmd(id <- mpi.comm.rank())
    mpi.bcast.cmd(np <- mpi.comm.size())
    mpi.bcast.cmd(host <- mpi.get.processor.name())
} else {
    library(parallel)
    n.cores <- 32
}

# Load nested list of cell indices defining groups of shapefile polygons
load("/workspace/UA/mfleonawicz/leonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN1km.RData")
load("/workspace/UA/mfleonawicz/leonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN1km_rmNA.RData")
if (exists("locgroup")) {
    cells_shp_list <- cells_shp_list[locgroup]
    cells_shp_list_rmNA <- cells_shp_list_rmNA[locgroup]
    region.names.out <- region.names.out[locgroup]
    n.shp <- sum(sapply(region.names.out, length))
}

dirs <- list.files("/big_scratch/apbennett/Calibration/FinalCalib", pattern = ".*.sres.*.", 
    full = T)
mainDirs <- rep(paste0(dirs, "/Maps")[model.index], each = length(reps))
modnames <- basename(dirname(mainDirs))  # Should only be one name at a time due to $ScenMod column naming below
dir.create(ageDenDir <- file.path("/big_scratch/mfleonawicz/Rmpi/outputs/ageDensities", 
    modnames[1]), recursive = T, showWarnings = F)  # also assumes one model at a time
repid <- rep(reps, length(model.index))
# breaks <- c(-1,5,25,50,100,200,10000) # Age class breaks # -1 inclusive of
# zero, not using cut(), using .bincode() n.brks <- length(breaks)
# age.labels <- c(paste(breaks[-c(n.brks-1, n.brks)] + 1, breaks[-c(1,
# n.brks)], sep='-'), paste0(breaks[n.brks-1], '+')) veg.labels <- c('Black
# Spruce', 'White Spruce', 'Deciduous', 'Shrub Tundra', 'Graminoid Tundra',
# 'Wetland Tundra', 'Barren lichen-moss', 'Temperate Rainforest')
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
mod.scen <- unlist(strsplit(modnames[1], "\\."))  # assume one model at a time
```

### Functions
Define support functions.


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
paste("Remaining support objects created. Now pushing objects to slaves.")
```

### Send information from head node to CPUs on compute nodes
Send **R** objects from master to slaves.


```r
# Export objects to slaves
if (Rmpi) {
    mpi.bcast.Robj2slave(veg.labels)
    mpi.bcast.Robj2slave(cells_shp_list)
    mpi.bcast.Robj2slave(cells_shp_list_rmNA)
    mpi.bcast.Robj2slave(region.names.out)
    mpi.bcast.Robj2slave(n.shp)
    mpi.bcast.Robj2slave(years)
    mpi.bcast.Robj2slave(mainDirs)
    mpi.bcast.Robj2slave(modnames)
    mpi.bcast.Robj2slave(ageDenDir)
    mpi.bcast.Robj2slave(repid)
    print("mpi.bcast.Robj2slave calls completed.")
}
```

Send **R** commands from master to slaves.


```r
# Issue commands to slaves
if (Rmpi) {
    mpi.bcast.cmd(mainDir <- mainDirs[id])
    mpi.bcast.cmd(source("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/code/alfStatsByRep.R"))
    mpi.bcast.cmd(dir.create(tmpDir <- paste0("/big_scratch/mfleonawicz/tmp/proc", 
        id), showWarnings = F))
    mpi.bcast.cmd(rasterOptions(chunksize = 1e+11, maxmemory = 1e+12, tmpdir = tmpDir))
    print("mpi.bcast.cmd calls completed. Now running mpi.remote.exec...")
} else {
    mainDir <- mainDirs[1]
    source("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/code/alfStatsByRep.R")
    tmpDir <- paste0("/big_scratch/mfleonawicz/tmp/procX")
    rasterOptions(chunksize = 1e+11, maxmemory = 1e+12, tmpdir = tmpDir)
}
```

### Gather data
Extract and compile fire statistics.


```r
# Compile fire statistics
if (doFire) {
    if (Rmpi) 
        abfc.fsv.dat <- mpi.remote.exec(getFireStats(i = repid[id], mainDir = mainDir, 
            years = years, cells.list = cells_shp_list, shp.names.list = region.names.out, 
            n = n.shp))
    if (!Rmpi) 
        abfc.fsv.dat <- mclapply(repid, getFireStats, mainDir = mainDir, years = years, 
            cells.list = cells_shp_list, shp.names.list = region.names.out, 
            n = n.shp, mc.cores = n.cores)
    print("Fire size by vegetation class completed.")
    
    abfc.dat <- rbindlist(lapply(abfc.fsv.dat, "[[", 1))
    fsv.dat <- rbindlist(lapply(abfc.fsv.dat, "[[", 2))
    
    rownames(abfc.dat) <- NULL
    abfc.dat$Model <- swapModelName(mod.scen[1])
    abfc.dat$Scenario <- swapScenarioName(mod.scen[2])
    abfc.dat$Phase <- getPhase(mod.scen[2])
    abfc.dat <- data.frame(abfc.dat)[, c(9:7, 1:6)]
    
    rownames(fsv.dat) <- NULL
    fsv.dat$Model <- swapModelName(mod.scen[1])
    fsv.dat$Scenario <- swapScenarioName(mod.scen[2])
    fsv.dat$Phase <- getPhase(mod.scen[2])
    fsv.dat <- data.frame(fsv.dat)[, c(9:7, 1:6)]
    
    print("Converted list to area burned and fire frequency data frame and fire size by vegetation class data frame.")
    print("Saving area burned and fire frequency data frames by location to .RData file.")
    
    locs <- unique(abfc.dat$Location)
    dir.create(abfcDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire/abfc", 
        recursive = TRUE, showWarnings = FALSE)
    
    for (j in 1:length(locs)) {
        filename.tmp <- paste0("abfc__", locs[j], "__", modnames[1])  # assume one model
        d.abfc <- abfc.dat[abfc.dat$Location == locs[j], ]
        save(d.abfc, locs, file = paste0(abfcDir, "/", filename.tmp, ".RData"))
        print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
    }
    rm(abfc.dat, d.abfc)
    gc()
    
    print("Saving fire size by vegetation class data frames by location to .RData file.")
    
    locs <- unique(fsv.dat$Location)
    dir.create(fsvDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire/fsv", recursive = TRUE, 
        showWarnings = FALSE)
    
    for (j in 1:length(locs)) {
        filename.tmp <- paste0("fsv__", locs[j], "__", modnames[1])  # assume one model
        d.fsv <- fsv.dat[fsv.dat$Location == locs[j], ]
        save(d.fsv, locs, file = paste0(fsvDir, "/", filename.tmp, ".RData"))
        print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
    }
    rm(fsv.dat, d.fsv)
    gc()
}
```

Extract and compile vegetation class and age statistics.


```r
# Compile vegetation class and age statistics
if (doAgeVeg) {
    if (Rmpi) 
        va.dat <- mpi.remote.exec(getAgeVegStats(i = repid[id], mainDir = mainDir, 
            denDir = ageDenDir, years = years, cells.list = cells_shp_list_rmNA, 
            shp.names.list = region.names.out, n = n.shp, veg.lab = veg.labels, 
            n.samples = 1000))
    if (!Rmpi) 
        va.dat <- mclapply(repid, getAgeVegStats, mainDir = mainDir, denDir = ageDenDir, 
            years = years, cells.list = cells_shp_list_rmNA, shp.names.list = region.names.out, 
            n = n.shp, n.samples = 1000, mc.cores = n.cores)
    print("Veg area data frame list returned from slaves.")
    
    va.dat <- rbindlist(va.dat)
    rownames(va.dat) <- NULL
    va.dat$Model <- swapModelName(mod.scen[1])
    va.dat$Scenario <- swapScenarioName(mod.scen[2])
    va.dat$Phase <- getPhase(mod.scen[2])
    va.dat <- data.frame(va.dat)[, c(10:8, 1:7)]
    print("Converted list to single veg area data frame.")
    print("Saving veg area data frames by location to .RData files.")
    
    locs <- unique(va.dat$Location)
    dir.create(vegDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/veg", showWarnings = F)
    
    for (j in 1:length(locs)) {
        filename.tmp <- paste0("veg__", locs[j], "__", modnames[1])  # assume one model
        d.veg <- va.dat[va.dat$Location == locs[j], ]
        save(d.veg, locs, file = paste0(vegDir, "/", filename.tmp, ".RData"))
        print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
    }
    rm(va.dat, d.veg)
    gc()
}

# All done
if (Rmpi) {
    mpi.close.Rslaves(dellog = FALSE)
    mpi.exit()
}
```
