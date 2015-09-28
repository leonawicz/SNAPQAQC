# @knitr setup
# Command line arguments
args=(commandArgs(TRUE))
if(!length(args)) q("no") else for(i in 1:length(args)) eval(parse(text=args[[i]]))

if(!exists("model.index")) stop("Must provide a model.index 1 to 15, e.g., model.index=1")
stopifnot(length(model.index)==1)
if(!exists("reps")) stop("Must provide replicates as integer(s) 1:200, e.g., reps=1:25")
if(!exists("years")) years <- NULL # Assume all years if none specified
if(!exists("Rmpi")) Rmpi <- TRUE
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
	cells <- filter(cells, LocGroup %in% locgroup)
}

dirs <- list.files("/big_scratch/apbennett/Calibration/FinalCalib", pattern=".*.sres.*.", full=T)
mainDirs <- rep(paste0(dirs,"/Maps")[model.index], each=length(reps))
modname <- unique(basename(dirname(mainDirs)))
dir.create(vegDir <- file.path("/big_scratch/mfleonawicz/Rmpi/outputs/veg", modname), recursive=T, showWarnings=F)
repid <- rep(reps,length(model.index))
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
	mpi.bcast.Robj2slave(years)
	mpi.bcast.Robj2slave(mainDirs)
	mpi.bcast.Robj2slave(modname)
	mpi.bcast.Robj2slave(vegDir)
	mpi.bcast.Robj2slave(repid)
	print("mpi.bcast.Robj2slave calls completed.")
}

# @knitr commands2slaves
# Issue commands to slaves
if(Rmpi){
	mpi.bcast.cmd( mainDir <- mainDirs[id] )
	mpi.bcast.cmd( source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfStatsByRep.R") )
	mpi.bcast.cmd( dir.create(tmpDir <- paste0("/big_scratch/mfleonawicz/tmp/proc",id), showWarnings=F) )
	mpi.bcast.cmd( rasterOptions(chunksize=10e10, maxmemory=10e11, tmpdir=tmpDir) )
	print("mpi.bcast.cmd calls completed. Now running mpi.remote.exec...")
} else {
	mainDir <- mainDirs[1]
	source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfStatsByRep.R")
	tmpDir <- paste0("/big_scratch/mfleonawicz/tmp/procX")
	rasterOptions(chunksize=10e10, maxmemory=10e11, tmpdir=tmpDir)
}

# @knitr fire_stats
# Compile fire statistics
if(doFire){
    print("#### Compiling fire statistics... ####")
	if(Rmpi){
        fsv.dat <- mpi.remote.exec( getFireStats(i=repid[id], mainDir=mainDir, years=years, cells=select(cells, -Cell_rmNA)) )
        fsv.dat <- rbindlist(fsv.dat)
    } else {
        len <- length(repid)
        if(len <= n.cores){
            fsv.dat <- mclapply(repid, getFireStats, mainDir=mainDir, years=years, cells=select(cells, -Cell_rmNA), mc.cores=n.cores)
            fsv.dat <- rbindlist(fsv.dat)
        } else {
            serial.iters <- ceiling(len/n.cores)
            n.cores2 <- which(len/(1:n.cores) < serial.iters)[1]
            fsv.dat <- vector("list", serial.iters)
            for(j in 1:serial.iters){
                repid.tmp <- 1:n.cores2 + (j-1)*n.cores2
                repid.tmp <- repid.tmp[repid.tmp <= max(repid)]
                fsv.tmp <- mclapply(repid.tmp, getFireStats, mainDir=mainDir, years=years, cells=select(cells, -Cell_rmNA), mc.cores=n.cores)
                fsv.dat[[j]] <- rbindlist(fsv.tmp)
                rm(fsv.tmp)
                gc()
                print(paste("Replicate batch", j, "of", serial.iters, "complete."))
            }
            fsv.dat <- rbindlist(fsv.dat)
        }
    }
	
	fsv.dat[, Model := swapModelName(mod.scen[1])]
	fsv.dat[, Scenario := swapScenarioName(mod.scen[2])]
	fsv.dat[, Scenario := factor(Scenario, levels=scen.levels)]
	fsv.dat[, Phase := getPhase(mod.scen[2])]
	fsv.dat <- setcolorder(fsv.dat, c("Phase", "Scenario", "Model", "LocGroup", "Location", "Vegetation", "FS", "FID", "Year", "Replicate"))
	setkey(fsv.dat, Location)
	
	print("Fire size by vegetation class completed.")
	print("Saving fire size by vegetation class data frames by location to .RData file.")
	locs <- unique(fsv.dat$Location)
	dir.create(fsvDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fsv", recursive=TRUE, showWarnings=FALSE)
	for(j in 1:length(locs)){
		filename.tmp <- paste0("fsv__", locs[j], "__", modname)
		d.fsv <- fsv.dat[locs[j]]
		save(d.fsv, locs, file=paste0(fsvDir, "/", filename.tmp, ".RData"))
		print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	}
	rm(fsv.dat, d.fsv)
	gc()
}

# @knitr age_veg_stats
# Compile vegetation class and age statistics
if(doAgeVeg){
    print("#### Compiling vegetation class and age statistics... ####")
	if(Rmpi){
        va.dat <- mpi.remote.exec( getAgeVegStats(i=repid[id], mainDir=mainDir, vegDir=vegDir, years=years, cells=select(cells, -Cell)) )
	} else {
        len <- length(repid)
        if(len <= n.cores){
            va.dat <- mclapply(repid, getAgeVegStats, mainDir=mainDir, vegDir=vegDir, years=years, cells=select(cells, -Cell), mc.cores=n.cores)
        } else {
            serial.iters <- ceiling(len/n.cores)
            n.cores2 <- which(len/(1:n.cores) < serial.iters)[1]
            va.dat <- vector("list", serial.iters)
            for(j in 1:serial.iters){
                repid.tmp <- 1:n.cores2 + (j-1)*n.cores2
                repid.tmp <- repid.tmp[repid.tmp <= max(repid)]
                va.dat.tmp <- mclapply(repid.tmp, getAgeVegStats, mainDir=mainDir, vegDir=vegDir, years=years, cells=select(cells, -Cell), mc.cores=n.cores)
                #print(sapply(va.dat.tmp, class))
                #err.ind <- which(sapply(va.dat.tmp, class)=="try-error")
                #if(length(err.ind)) for(e in 1:length(err.ind)) print(va.dat.tmp[[e]])
                #print(va.dat.tmp[[1]])
                #va.dat[[j]] <- rbindlist(va.dat.tmp)
                rm(va.dat.tmp)
                gc()
                print(paste("Replicate batch", j, "of", serial.iters, "complete."))
            }
        }
    }

	#va.dat <- rbindlist(va.dat)
	#va.dat[, Model := swapModelName(mod.scen[1])]
	#va.dat[, Scenario := swapScenarioName(mod.scen[2])]
	#va.dat[, Scenario := factor(Scenario, levels=scen.levels)]
	#va.dat[, Phase := getPhase(mod.scen[2])]
	#va.dat <- setcolorder(va.dat, c("Phase", "Scenario", "Model", "LocGroup", "Location", "Vegetation", "Age", "Area", "Year", "Replicate"))
	#setkey(va.dat, Location)
	print("Vegetation area by age completed.")
	#print("Saving veg area by age data tables by location to .RData files.")

	#locs <- unique(va.dat$Location)
	#dir.create(vegAreaDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/vegArea", showWarnings=F)

	#for(j in 1:length(locs)){
	#	filename.tmp <- paste0("veg__", locs[j], "__", modname)
	#	d.veg <- va.dat[locs[j]]
	#	save(d.veg, locs, file=paste0(vegDir, "/", filename.tmp, ".RData"))
	#	print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	#}
	#print(tables())
	#rm(va.dat, d.veg)
	#gc()
}

# All done
if(Rmpi){
	mpi.close.Rslaves(dellog = FALSE)
	mpi.exit()
}
