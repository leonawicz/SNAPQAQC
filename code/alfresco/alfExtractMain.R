# @knitr setup
# Command line arguments
args=(commandArgs(TRUE))
if(!length(args)) q("no") else for(i in 1:length(args)) eval(parse(text=args[[i]]))

if(!exists("model.index")) stop("Must provide a model.index 1 to 15, e.g., model.index=1")
stopifnot(length(model.index)==1)
if(!exists("reps")) stop("Must provide replicates as integer(s) 1:200, e.g., reps=1:25")
if(!exists("years")) years <- 2008:2100 # Assume if not specified
if(!exists("Rmpi")) Rmpi <- TRUE
if(Rmpi) {
    if(!exists("mpiBy")) mpiBy <- "year"
    itervar <- if(mpiBy=="rep") reps else if(mpiBy=="year") 1:length(years) else stop("mpiBy must be 'rep' or 'year'.")
    loopBy <- if(mpiBy=="rep") "year" else "rep"
}
if(!exists("doFire")) doFire <- TRUE
if(!exists("doAgeVeg")) doAgeVeg <- TRUE
if(exists("repSample") && is.numeric(repSample)){
	set.seed(47)
	reps <- sort(sample(reps, min(repSample, length(reps))))
	cat("Sampled replicates:\n", reps, "\n")
}

library(raster)
library(data.table)
library(dplyr)

# Rmpi setup
if(Rmpi){
	library(Rmpi)
	mpi.spawn.Rslaves(needlog = TRUE)
	mpi.bcast.cmd( id <- mpi.comm.rank() )
	mpi.bcast.cmd( np <- mpi.comm.size() )
	mpi.bcast.cmd( host <- mpi.get.processor.name() )
} else {
	library(parallel)
	n.cores <- 32
}

# Load data table of cell indices defining groups of shapefile polygons
load("/workspace/UA/mfleonawicz/projects/DataExtraction/workspaces/shapes2cells_akcan1km2km.RData")
cells <- filter(cells, Source=="akcan1km") %>% group_by %>% select(-Source) %>% group_by(LocGroup, Location)
if(exists("locgroup")){
    cat("locgroup = "); cat(locgroup); cat("\n")
	cells <- filter(cells, LocGroup %in% unique(cells$LocGroup)[locgroup])
    print(unique(cells$LocGroup))
    stopifnot(nrow(cells) > 0)
}

dirs <- list.files("/big_scratch/apbennett/Calibration/FinalCalib", pattern=".*.sres.*.", full=T)
mainDirs <- rep(paste0(dirs,"/Maps")[model.index], each=length(itervar))
modname <- unique(basename(dirname(mainDirs)))
if(mpiBy=="rep") dir.create(ageDir <- file.path("/big_scratch/mfleonawicz/Rmpi/outputs/veg", modname), recursive=T, showWarnings=F) else ageDir <- NULL

#veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
mod.scen <- unlist(strsplit(modname, "\\."))

# @knitr functions
# Support functions
# include later: # CCSM4="CCSM4", GFDLcm3="GFDLcm3", GISSe2-r="GISSe2-r", IPSLcm5a-lr="IPSLcm5a-lr", MRIcgcm3="MRIcgcm3"
swapModelName <- function(x){
	switch(x,
		cccma_cgcm3_1="CCCMAcgcm31", gfdl_cm2_1="GFDLcm21", miroc3_2_medres="MIROC32m", mpi_echam5="MPIecham5", ukmo_hadcm3="ukmoHADcm3"
	)
}

# include later: # rcp45="RCP 4.5", rcp60="RCP 6.0", rcp85="RCP 8.5"
swapScenarioName <- function(x){
	switch(x,
		sresb1="SRES B1", sresa1b="SRES A1B", sresa2="SRES A2"
	)
}

# include later: # rcp45="AR5", rcp60="AR5", rcp85="AR5"
getPhase <- function(x){
	switch(x,
		sresb1="AR4",	sresa1b="AR4", sresa2="AR4"
	)
}
paste("Remaining support objects created. Now pushing objects to slaves.")

# @knitr obj2slaves
# Export objects to slaves
if(Rmpi){
	mpi.bcast.Robj2slave(cells)
    mpi.bcast.Robj2slave(reps)
	mpi.bcast.Robj2slave(years)
	mpi.bcast.Robj2slave(mainDirs)
	mpi.bcast.Robj2slave(modname)
	mpi.bcast.Robj2slave(ageDir)
	mpi.bcast.Robj2slave(itervar)
    mpi.bcast.Robj2slave(loopBy)
	print("mpi.bcast.Robj2slave calls completed.")
}

# @knitr commands2slaves
# Issue commands to slaves
if(Rmpi){
	mpi.bcast.cmd( mainDir <- mainDirs[id] )
	mpi.bcast.cmd( source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfExtract.R") )
	mpi.bcast.cmd( dir.create(tmpDir <- paste0("/big_scratch/mfleonawicz/tmp/proc",id), showWarnings=F) )
	mpi.bcast.cmd( rasterOptions(chunksize=10e10, maxmemory=10e11, tmpdir=tmpDir) )
	print("mpi.bcast.cmd calls completed. Now running mpi.remote.exec...")
} else {
	mainDir <- mainDirs[1]
	source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfExtract.R")
	tmpDir <- paste0("/big_scratch/mfleonawicz/tmp/procX")
	rasterOptions(chunksize=10e10, maxmemory=10e11, tmpdir=tmpDir)
}

# @knitr fire_stats
# Compile fire statistics
if(doFire){
    print("#### Compiling fire statistics... ####")
	if(Rmpi){
        fsv.dat <- mpi.remote.exec( extract_data(i=itervar[id], type="fsv", loopBy=loopBy, mainDir=mainDir, reps=reps, years=years, cells=select(cells, -Cell_rmNA)) )
        fsv.dat <- rbindlist(fsv.dat)
    } else {
        len <- length(itervar)
        if(len <= n.cores){
            fsv.dat <- mclapply(itervar, extract_data, type="fsv", loopBy=loopBy, mainDir=mainDir, reps=reps, years=years, cells=select(cells, -Cell_rmNA), mc.cores=n.cores)
            fsv.dat <- rbindlist(fsv.dat)
        } else {
            serial.iters <- ceiling(len/n.cores)
            n.cores2 <- which(len/(1:n.cores) < serial.iters)[1]
            fsv.dat <- vector("list", serial.iters)
            for(j in 1:serial.iters){
                itervar.tmp <- 1:n.cores2 + (j-1)*n.cores2
                itervar.tmp <- itervar.tmp[itervar.tmp <= max(itervar)]
                fsv.tmp <- mclapply(itervar.tmp, extract_data, type="fsv", loopBy=loopBy, mainDir=mainDir, reps=reps, years=years, cells=select(cells, -Cell_rmNA), mc.cores=n.cores)
                fsv.dat[[j]] <- rbindlist(fsv.tmp)
                rm(fsv.tmp)
                gc()
                print(paste("Replicate batch", j, "of", serial.iters, "complete."))
            }
            fsv.dat <- rbindlist(fsv.dat)
        }
    }
	
    fsv.dat.names.ini <- copy(names(fsv.dat))
	fsv.dat[, Model := swapModelName(mod.scen[1])]
	fsv.dat[, Scenario := swapScenarioName(mod.scen[2])]
	fsv.dat[, Scenario := factor(Scenario, levels=scen.levels)]
	fsv.dat[, Phase := getPhase(mod.scen[2])]
	fsv.dat <- setcolorder(fsv.dat, c("Phase", "Scenario", "Model", fsv.dat.names.ini))
	setkey(fsv.dat, Location)
	
	print("Fire size by vegetation class completed.")
	print("Saving fire size by vegetation class data frames by location to .RData file.")
	locs <- unique(fsv.dat$Location)
	dir.create(fsvDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fsv", recursive=TRUE, showWarnings=FALSE)
	for(j in 1:length(locs)){
		filename.tmp <- paste0("fsv__", locs[j], "__", modname)
		d.fsv <- fsv.dat[locs[j]]
		save(d.fsv, file=paste0(fsvDir, "/", filename.tmp, ".RData"))
		print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	}
    print(tables())
	rm(fsv.dat, d.fsv)
	gc()
}

# @knitr age_veg_stats
# Compile vegetation class and age statistics
if(doAgeVeg){
    print("#### Compiling vegetation class and age statistics... ####")
	if(Rmpi){
        va.dat <- mpi.remote.exec( extract_data(i=itervar[id], type="av", loopBy=loopBy, mainDir=mainDir, ageDir=ageDir, reps=reps, years=years, cells=select(cells, -Cell)) )
        d.area <- rbindlist(lapply(va.dat, function(x) x$d.area))
        if(mpiBy=="year") d.age <- rbindlist(lapply(va.dat, function(x) x$d.age))
	} else {
        len <- length(itervar)
        if(len <= n.cores){
            va.dat <- mclapply(itervar, extract_data, type="av", loopBy=loopBy, mainDir=mainDir, ageDir=ageDir, reps=reps, years=years, cells=select(cells, -Cell), mc.cores=n.cores)
            d.area <- rbindlist(lapply(va.dat, function(x) x$d.area))
            if(mpiBy=="year") d.age <- rbindlist(lapply(va.dat, function(x) x$d.age))
        } else {
            serial.iters <- ceiling(len/n.cores)
            n.cores2 <- which(len/(1:n.cores) < serial.iters)[1]
            d.age <- d.area <- vector("list", serial.iters)
            for(j in 1:serial.iters){
                itervar.tmp <- 1:n.cores2 + (j-1)*n.cores2
                itervar.tmp <- itervar.tmp[itervar.tmp <= max(itervar)]
                va.dat <- mclapply(itervar.tmp, extract_data, type="av", loopBy=loopBy, mainDir=mainDir, ageDir=ageDir, reps=reps, years=years, cells=select(cells, -Cell), mc.cores=n.cores)
                d.area[[j]] <- rbindlist(lapply(va.dat, function(x) x$d.area))
                if(mpiBy=="year") d.age[[j]] <- rbindlist(lapply(va.dat, function(x) x$d.age))
                rm(va.dat)
                gc()
                print(paste("Replicate batch", j, "of", serial.iters, "complete."))
            }
            d.area <- rbindlist(d.area)
            if(mpiBy=="year") d.age <- rbindlist(d.age)
        }
    }

    d.area.names.ini <- copy(names(d.area))
    d.area[, Model := swapModelName(mod.scen[1])]
	d.area[, Scenario := swapScenarioName(mod.scen[2])]
	d.area[, Scenario := factor(Scenario, levels=scen.levels)]
	d.area[, Phase := getPhase(mod.scen[2])]
	d.area <- setcolorder(d.area, c("Phase", "Scenario", "Model", d.area.names.ini))
	setkey(d.area, Location)
    print("Vegetation area completed.")
    print("Saving vegetation area data tables by location to .RData files.")
	locs <- unique(d.area$Location)
	dir.create(vegDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/veg", showWarnings=F)

	for(j in 1:length(locs)){
		filename.tmp <- paste0("veg__", locs[j], "__", modname)
		d.vegarea <- d.area[locs[j]]
		save(d.vegarea, file=paste0(vegDir, "/", filename.tmp, ".RData"))
		print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	}
	print(tables())
	rm(d.area, d.vegarea)
	gc()
    
    if(mpiBy=="year"){
    d.age.names.ini <- copy(names(d.age))
    d.age[, Model := swapModelName(mod.scen[1])]
    d.age[, Scenario := swapScenarioName(mod.scen[2])]
    d.age[, Scenario := factor(Scenario, levels=scen.levels)]
    d.age[, Phase := getPhase(mod.scen[2])]
    d.age <- setcolorder(d.age, c("Phase", "Scenario", "Model", d.age.names.ini))
    setkey(d.age, Location)
    print("Vegetation area by age completed.")
    print("Saving vegetation area by age data tables by location to .RData files.")
	locs <- unique(d.age$Location)
	dir.create(ageDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/age", showWarnings=F)

	for(j in 1:length(locs)){
		filename.tmp <- paste0("age__", locs[j], "__", modname)
		d.vegage <- d.age[locs[j]]
		save(d.vegage, file=paste0(ageDir, "/", filename.tmp, ".RData"))
		print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	}
	print(tables())
	rm(d.age, d.vegage)
	gc()
    }
}

# All done
if(Rmpi){
	mpi.close.Rslaves(dellog = FALSE)
	mpi.exit()
}
