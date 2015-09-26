# @knitr setup
# Command line arguments
args=(commandArgs(TRUE))
if(!length(args)) q("no") else for(i in 1:length(args)) eval(parse(text=args[[i]]))

if(!exists("model.index")) stop("Must provide model(s) as integer(s) 1:15, e.g., model.index=1:2")
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

library(data.table)

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

# Load nested list of cell indices defining groups of shapefile polygons
load("/workspace/UA/mfleonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN1km.RData")
load("/workspace/UA/mfleonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN1km_rmNA.RData")
if(exists("locgroup")){
    print(paste("locgroup =", locgroup))
	cells_shp_list <- cells_shp_list[locgroup]
	cells_shp_list_rmNA <- cells_shp_list_rmNA[locgroup]
	region.names.out <- region.names.out[locgroup]
	n.shp <- sum(sapply(region.names.out, length))
}

dirs <- list.files("/big_scratch/apbennett/Calibration/FinalCalib/b1_rerun", pattern=".*.sres.*.", full=T)
mainDirs <- rep(paste0(dirs,"/Maps")[model.index], each=length(reps))
modnames <- basename(dirname(mainDirs)) # Should only be one name at a time due to $ScenMod column naming below
dir.create(ageDenDir <- file.path("/big_scratch/mfleonawicz/Rmpi/outputs/ageDensities", modnames[1]), recursive=T, showWarnings=F) # also assumes one model at a time
repid <- rep(reps,length(model.index))
#breaks <- c(-1,5,25,50,100,200,10000) # Age class breaks # -1 inclusive of zero, not using cut(), using .bincode()
#n.brks <- length(breaks)
#age.labels <- c(paste(breaks[-c(n.brks-1, n.brks)] + 1, breaks[-c(1, n.brks)], sep="-"), paste0(breaks[n.brks-1], "+"))
#veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
mod.scen <- unlist(strsplit(modnames[1], "\\.")) # assume one model at a time

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
        abfc.fsv.dat <- mpi.remote.exec( getFireStats(i=repid[id], mainDir=mainDir, years=years, cells.list=cells_shp_list, shp.names.list=region.names.out, n=n.shp) )
        abfc.dat <- rbindlist(lapply(abfc.fsv.dat, "[[", 1))
        fsv.dat <- rbindlist(lapply(abfc.fsv.dat, "[[", 2))
    } else {
        len <- length(repid)
        if(len <= n.cores){
            abfc.fsv.dat <- mclapply(repid, getFireStats, mainDir=mainDir, years=years, cells.list=cells_shp_list, shp.names.list=region.names.out, n=n.shp, mc.cores=n.cores)
            abfc.dat <- rbindlist(lapply(abfc.fsv.dat, "[[", 1))
            fsv.dat <- rbindlist(lapply(abfc.fsv.dat, "[[", 2))
        } else {
            serial.iters <- ceiling(len/n.cores)
            n.cores2 <- which(len/(1:n.cores) < serial.iters)[1]
            abfc.dat <- fsv.dat <- vector("list", serial.iters)
            for(j in 1:serial.iters){
                repid.tmp <- 1:n.cores2 + (j-1)*n.cores2
                repid.tmp <- repid.tmp[repid.tmp <= max(repid)]
                abfc.fsv <- mclapply(repid.tmp, getFireStats, mainDir=mainDir, years=years, cells.list=cells_shp_list, shp.names.list=region.names.out, n=n.shp, mc.cores=n.cores)
                abfc.dat[[j]] <- rbindlist(lapply(abfc.fsv, "[[", 1))
                fsv.dat[[j]] <- rbindlist(lapply(abfc.fsv, "[[", 2))
                rm(abfc.fsv)
                gc()
                print(paste("Replicate batch", j, "of", serial.iters, "complete."))
            }
            abfc.dat <- rbindlist(abfc.dat)
            fsv.dat <- rbindlist(fsv.dat)
        }
    }
	print("Fire size by vegetation class completed.")
    
	abfc.dat[, Model := swapModelName(mod.scen[1])]
	abfc.dat[, Scenario := swapScenarioName(mod.scen[2])]
	abfc.dat[, Scenario := factor(Scenario, levels=scen.levels)]
	abfc.dat[, Phase := getPhase(mod.scen[2])]
	abfc.dat <- setcolorder(abfc.dat, c("Phase", "Scenario", "Model", "LocGroup", "Location", "Var", "Val", "Year", "Replicate"))
	setkey(abfc.dat, Location)
	
	fsv.dat[, Model := swapModelName(mod.scen[1])]
	fsv.dat[, Scenario := swapScenarioName(mod.scen[2])]
	fsv.dat[, Scenario := factor(Scenario, levels=scen.levels)]
	fsv.dat[, Phase := getPhase(mod.scen[2])]
	fsv.dat <- setcolorder(fsv.dat, c("Phase", "Scenario", "Model", "LocGroup", "Location", "VegID", "Val", "FID", "Year", "Replicate"))
	setkey(fsv.dat, Location)
	
	print("Converted list to area burned and fire frequency data table and fire size by vegetation class data table.")
	print("Saving area burned and fire frequency data tables by location to .RData files.")

	locs <- unique(abfc.dat$Location)
	dir.create(abfcDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire/abfc", recursive=TRUE, showWarnings=FALSE)

	for(j in 1:length(locs)){
		filename.tmp <- paste0("abfc__", locs[j], "__", modnames[1]) # assume one model
		d.abfc <- abfc.dat[locs[j]]
		save(d.abfc, locs, file=paste0(abfcDir, "/", filename.tmp, ".RData"))
		print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	}
	print(tables())
	rm(abfc.dat, d.abfc)
	gc()

	print("Saving fire size by vegetation class data frames by location to .RData file.")

	locs <- unique(fsv.dat$Location)
	dir.create(fsvDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire/fsv", recursive=TRUE, showWarnings=FALSE)

	for(j in 1:length(locs)){
		filename.tmp <- paste0("fsv__", locs[j], "__", modnames[1]) # assume one model
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
        va.dat <- mpi.remote.exec( getAgeVegStats(i=repid[id], mainDir=mainDir, denDir=ageDenDir, years=years, cells.list=cells_shp_list_rmNA, shp.names.list=region.names.out, n=n.shp, n.samples=100) )
	} else {
        len <- length(repid)
        if(len <= n.cores){
            va.dat <- mclapply(repid, getAgeVegStats, mainDir=mainDir, denDir=ageDenDir, years=years, cells.list=cells_shp_list_rmNA, shp.names.list=region.names.out, n=n.shp, n.samples=100, mc.cores=n.cores)
        } else {
            serial.iters <- ceiling(len/n.cores)
            n.cores2 <- which(len/(1:n.cores) < serial.iters)[1]
            va.dat <- vector("list", serial.iters)
            for(j in 1:serial.iters){
                repid.tmp <- 1:n.cores2 + (j-1)*n.cores2
                repid.tmp <- repid.tmp[repid.tmp <= max(repid)]
                va.dat.tmp <- mclapply(repid.tmp, getAgeVegStats, mainDir=mainDir, denDir=ageDenDir, years=years, cells.list=cells_shp_list_rmNA, shp.names.list=region.names.out, n=n.shp, n.samples=100, mc.cores=n.cores)
                print(sapply(va.dat.tmp, class))
                err.ind <- which(sapply(va.dat.tmp, class)=="try-error")
                if(length(err.ind)) for(e in 1:length(err.ind)) print(va.dat.tmp[[e]])
                print(va.dat.tmp[[1]])
                va.dat[[j]] <- rbindlist(va.dat.tmp)
                rm(va.dat.tmp)
                gc()
                print(paste("Replicate batch", j, "of", serial.iters, "complete."))
            }
        }
    }
    print("Veg area data frame list returned from slaves.")

	va.dat <- rbindlist(va.dat)
	va.dat[, Model := swapModelName(mod.scen[1])]
	va.dat[, Scenario := swapScenarioName(mod.scen[2])]
	va.dat[, Scenario := factor(Scenario, levels=scen.levels)]
	va.dat[, Phase := getPhase(mod.scen[2])]
	va.dat <- setcolorder(va.dat, c("Phase", "Scenario", "Model", "LocGroup", "Location", "VegID", "Var", "Val", "Year", "Replicate"))
	setkey(va.dat, Location)
	print("Converted list to single veg area data table.")
	print("Saving veg area data tables by location to .RData files.")

	locs <- unique(va.dat$Location)
	dir.create(vegDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/veg", showWarnings=F)

	for(j in 1:length(locs)){
		filename.tmp <- paste0("veg__", locs[j], "__", modnames[1]) # assume one model
		d.veg <- va.dat[locs[j]]
		save(d.veg, locs, file=paste0(vegDir, "/", filename.tmp, ".RData"))
		print(paste(filename.tmp, "object", j, "of", length(locs), "saved."))
	}
	print(tables())
	rm(va.dat, d.veg)
	gc()
}

# All done
if(Rmpi){
	mpi.close.Rslaves(dellog = FALSE)
	mpi.exit()
}
