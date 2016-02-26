setwd("/atlas_scratch/mfleonawicz/projects/SNAPQAQC/workspaces")

library(parallel)
library(raster)
tmpDir <- paste0("/atlas_scratch/mfleonawicz/tmp")
rasterOptions(chunksize=10e10, maxmemory=10e11, tmpdir=tmpDir)

# setup
prism.yrs <- 1961:1990
#cru.yrs <- 1970:1999
gcm.yrs <- 2010:2039
topDir <- "/Data/Base_Data/Climate/AK_CAN_2km"
#topDir <- "/Data/Base_Data/Climate/World/World_10min" # for world grid 10-minute res data
vars <- c("tas", "pr")
scenarios <- c("sresb1", "sresa1b", "sresa2")
rcps <- c("rcp45", "rcp60", "rcp85")
ar4models <- c("cccma-cgcm3-1-t47", "gfdl-cm2-1", "miroc3-2-medres", "mpi-echam5", "ukmo-hadcm3")
ar5models <- c("CCSM4", "GFDL-CM3", "GISS-E2-R", "IPSL-CM5A-LR", "MRI-CGCM3")
n.s <- length(scenarios)
n.r <- length(rcps)
n.m1 <- length(ar4models)
n.m2 <- length(ar5models)
seasons <- list(Winter=1:3, Spring=4:6, Summer=7:9, Fall=10:12)
prism6190dirs <- file.path(topDir, "historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs", vars, "ak83albers")
#cru32dirs <- file.path(topDir, "historical/CRU/CRU_TS32", vars)
ar4dirs <- file.path(topDir, "projected/AR4_CMIP3_models", rep(scenarios, each=n.m1), ar4models, rep(vars, each=n.s*n.m1))
ar5dirs <- file.path(topDir, "projected/AR5_CMIP5_models", rep(rcps, each=n.m2), ar5models, rep(vars, each=n.r*n.m2))
outDir <- paste0("../data/delta_maps/climate_2km/base_climatology_PRISM_", paste(range(prism.yrs), collapse="_"))
#outDir <- paste0("../data/delta_maps/climate_10min/base_climatology_CRU_", paste(range(cru.yrs), collapse="_"))

# functions
# monthly or seasonal (DJF, MAM, JJA, SON) climatology layers
get_clim <- function(k=NULL, dir, years, seasons=FALSE, prev.december=TRUE, round.values=FALSE, internal.par=FALSE){
    if(!is.null(k)){ dir <- dir[k]; years <- years[[k]] }
    path.vec <- strsplit(dir, "/")[[1]]
    prism <- "prism" %in% path.vec
    rnd <- if("tas" %in% path.vec) 1 else 0
    p <- if(seasons) list(Winter=c(1,2,12), Spring=3:5, Summer=6:8, Fall=9:11) else
        if(prev.december) split(1:13, factor(c("PrevDec", month.abb), levels=c("PrevDec", month.abb))) else
            split(1:12, factor(month.abb, levels=month.abb))
    files <- list.files(dir, pattern="\\.tif$", full=T)
    mo.vec <- if(prism) c(paste0("0", 1:9), 10:12) else substr(files, nchar(files)-10, nchar(files)-9)
    files <- split(files, mo.vec)
    if(!seasons && prev.december) files <- c(files[12], files)
    if(!prism){
    files <- lapply(1:length(p), function(i, x, years, seasons, pre){
        yr.vec <- as.numeric(substr(x[[i]], nchar(x[[i]])-7, nchar(x[[i]])-4))
        if((seasons && pre && i==12) || (!seasons && pre && i==1)) years <- years-1
        x[[i]][yr.vec %in% years]
        }, x=files, years=years, seasons=seasons, pre=prev.december)
    }
    files <- lapply(p, function(x, files) unlist(files[x]), files=files)
    if(internal.par) b <- brick(stack(mclapply(files, function(x) calc(stack(x, quick=T), mean), mc.cores=length(files))))
    if(!internal.par) b <- brick(stack(lapply(files, function(x) calc(stack(x, quick=T), mean))))
    if(round.values) b <- round(b, rnd)
    names(b) <- names(p)
    b
}

# append seasons to monthly layers
append_seasons <- function(x, seasons){
    n <- min(32, length(x))
    y <- mclapply(x, function(x, p) brick(lapply(p, function(p, x) calc(subset(x, p), mean), x=x)), p=seasons, mc.cores=n)
    gc()
    y <- mclapply(1:length(x), function(i, x, y) stack(x[[i]], y[[i]]), x=x, y=y, mc.cores=n)
    names(y) <- names(x)
    gc()
    y
}

# calculate deltas from historical and future climatologies
get_deltas <- function(x, y, v=NULL, n.cores=32){
    id <- names(y)
    if(is.null(v)) v <- factor(basename(id), levels=unique(basename(id)))
    y <- split(y, v)
    gc()
    for(i in 1:length(y)){
        nam <- names(y)[i]
        ind <- match(nam, unique(v))
        x1 <- x[[ind]]
        gc()
        y[[i]] <- if(nam=="tas") mclapply(y[[i]], function(y, x) round(y-x, 3), x=x1, mc.cores=n.cores) else
            if(nam=="pr") mclapply(y[[i]], function(y, x) round(y/x, 3), x=x1, mc.cores=n.cores)
        gc()
    }
    y <- unlist(y)
    names(y) <- id
    gc()
    y
}

# save tifs
write_delta_tifs <- function(deltas, outDir, years.range, scenarios, models, vars, n.cores=10, ...){
    n <- length(deltas)
    dir.create(outDir, recur=T, showWarnings=F)
    mclapply(1:n, function(i, ...){
        deltas <- deltas[i]
        gc()
        a <- strsplit(names(deltas), "/")[[1]]
        p <- intersect(c("ar4", "ar5"), tolower(substr(a, 1, 3)))
        s <- intersect(scenarios, a)
        m <- intersect(models, a)
        v <- intersect(vars, a)
        yrs <- paste(years.range, collapse="-")
        out <- paste0(outDir, paste("/deltas", v, p, s, m, yrs, sep="_"), ".tif")
        writeRaster(deltas[[1]], out, datatype="FLT4S", overwrite=T)
        gc()
    }, mc.cores=n.cores)
}

# process
#n <- length(cru32dirs)
#cru.clim <- lapply(1:length(cru32dirs), get_clim, cru32dirs, rep(list(cru.yrs), n), internal.par=TRUE) # 2 min
#names(cru.clim) <- cru32dirs
####cru.clim <- lapply(1:n, get_clim, cru32dirs, rep(list(cru.yrs), n), seasons=TRUE, internal.par=TRUE) # 3 min, not run
#cru.clim <- append_seasons(cru.clim, seasons)

n <- length(prism6190dirs)
prism.clim <- lapply(1:length(prism6190dirs), get_clim, prism6190dirs, rep(list(prism.yrs), n), internal.par=TRUE) # 2 min
prism.clim <- mclapply(prism.clim, trim, mc.cores=2)
names(prism.clim) <- prism6190dirs
#prism.clim <- lapply(1:n, get_clim, prism6190dirs, rep(list(prism.yrs), n), seasons=TRUE, internal.par=TRUE) # 3 min, not run
prism.clim <- append_seasons(prism.clim, seasons)

n <- length(ar4dirs)
y <- mclapply(1:n, get_clim, ar4dirs, rep(list(gcm.yrs), n), mc.cores=32) # 16.5 min
names(y) <- ar4dirs
#y <- mclapply(1:n, get_clim, ar4dirs, rep(list(gcm.yrs), n), seasons=TRUE, mc.cores=32) # 10 min, not run
y <- append_seasons(y, seasons)
gc()
y <- get_deltas(x=prism.clim, y) # 1 min
gc()
write_delta_tifs(y, outDir, range(gcm.yrs), scenarios, ar4models, vars)

n <- length(ar5dirs)
y <- mclapply(1:n, get_clim, ar5dirs, rep(list(gcm.yrs), n), mc.cores=32) # 16.5 min
names(y) <- ar5dirs
#y <- mclapply(1:n, get_clim, ar5dirs, rep(list(gcm.yrs), n), seasons=TRUE, mc.cores=32) # 10 min, not run
y <- append_seasons(y, seasons)
gc()
y <- get_deltas(x=prism.clim, y) # 1 min
gc()
write_delta_tifs(y, outDir, range(gcm.yrs), scenarios, ar5models, vars)
