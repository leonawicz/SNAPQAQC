##############################################################################
#### Stage two compilation of extracted Alfresco outputs into data tables ####
##############################################################################

#### Script author:  Matthew Leonawicz ####
#### Maintainted by: Matthew Leonawicz ####
#### Last updated:   07/22/2015        ####

# @knitr setup
comargs <- (commandArgs(TRUE))
if(!length(comargs)) q("no") else for(z in 1:length(comargs)) eval(parse(text=comargs[[z]]))

library(parallel)
library(reshape2)
library(data.table)
library(dplyr)

if(!exists("mainDir")) mainDir <- "/big_scratch/mfleonawicz/Rmpi/outputs"
if(!exists("variable")) stop("Must provide 'variable' argument in escaped quotes. Options are 'age' (vegetation age), 'veg' (vegetated area), 'abfc' (area burned and fire counts), or 'fsv' (fire sizes by vegetation class).")
ageDirs <- list.files(file.path(mainDir, "ageDensities"), full=T)
abfcDir <- file.path(mainDir, "fire/abfc")
fsvDir <- file.path(mainDir, "fire/fsv")
vegDir <- file.path(mainDir, "veg")

n.samples <- 1000 # new density estimation
n.samples.in <- 100 # known, based on upstream settings for vegetation age
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
datnames <- c("Phase", "Scenario", "Model", "Location", "Var", "Vegetation", "Year", "Val")

# All regions for which stage-1 outputs currently exist on disk.
# Checking abfc files, which have maximum regions, but only regions common to fsv, abfc, veg, and age outputs are processed.
# It is assumed stage one processing has been done on a common set of regions for all four variable sets.
regions <- unique(sapply(strsplit(list.files(file.path(abfcDir)), "__"), "[", 2))
n.regions <- length(regions)
n.cores <- min(n.regions, 32)

# @knitr functions1
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
		sresb1="AR4", sresa1b="AR4", sresa2="AR4"
	)
}

# @knitr functions2
# Density estimation
denFun <- function(x, n=1000, adj=0.25, min.zero=TRUE, diversify=FALSE, missing.veg.NA=TRUE, fire=FALSE){
	if(all(is.na(x))) return(rep(NA, 2*n))
	x <- x[!is.na(x)]
	lx <- length(x)
	if(sum(x)==0 & missing.veg.NA & !fire) return(rep(NA, 2*n))
	if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
	if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
	dif <- diff(range(x))
	z <- density(x, adjust=adj, n=n, from=min(x)-max(1, 0.05*dif), to=max(x)+max(1, 0.05*dif))
	if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+max(1, 0.05*dif))
	as.numeric(c(z$x, z$y))
}

# Bootstrapping
btfun_dt <- function(p1, p2, n.samples=length(p)/2, n.boot=10000, interp=FALSE, n.interp=1000, ...){
    p <- c(p1, p2) # concatenate separate columns from data table
	if(!length(p)) return(NULL)
	if(all(is.na(p))) return(rep(NA, n.boot))
	p <- list(x=p[1:n.samples], y=p[(n.samples+1):(2*n.samples)])
	if(interp && length(unique(p[1:n.samples])) > 1) p <- approx(p$x, p$y, n=n.interp)
	p <- round(sample(p$x, n.boot, prob=p$y, rep=T), ...)
	p
}

# @knitr not_in_use
# additional support functions, not in use
#getSamples <- function(d, v) d[d[,2]==v, 1]

#getStats <- function(x, seq.q){
#	if(all(is.na(x))) return(rep(NA, 2 + length(seq.q)))
#	x1 <- round(mean(x, na.rm=TRUE))
#	x2 <- round(sd(x, na.rm=TRUE), 1)
#	x3 <- round(quantile(x, probs=seq.q, na.rm=TRUE))
#	c(x1, x2, x3)
#}

# @knitr get_ageStats_ageDensities
# Primary processing functions
get_ageStats_ageDensities <- function(j, dirs, n.samples=100, n.samples.in, n.boot=10000, datnames, veg.labels, scen.levels){
	pat <- paste0("^age__.*.rep.*.RData$")
	pat2 <- substr(pat, 1, nchar(pat) - 9)
	files.list <- lapply(1:length(dirs), function(i, dirs, ...) list.files(dirs[i], ...), dirs=dirs, full=T, pattern=pat)
	for(i in 1:length(files.list)){
		files.locs <- sapply(strsplit(files.list[[i]], "__"), "[", 2)
		locs <- unique(files.locs)
        if(j > length(locs)) return(NULL)
        loc <- locs[j]
        if(loc %in% c("Manitoba", "Saskatchewan")) return(NULL)
		files.list[[i]] <- files.list[[i]][which(files.locs %in% loc)]
	}
    dat.list <- stat.list <- vector("list", length(dirs))
	for(i in 1:length(dirs)){
		modname <- basename(dirs[i])
        mod.scen <- unlist(strsplit(modname, "\\."))
		# Load all reps for one location, all model and scenarios, into local location-specific environment
		print(system.time( for(p in 1:length(files.list[[i]])) { load(files.list[[i]][p], envir=environment()); print(paste("Local environment: loading workspace", p, "of", length(files.list[[i]]))) } ))
		d.names <- ls(pattern=pat2, envir=environment())
		d.list <- mget(d.names, envir=environment())
		print(paste("EVIRONMENT", paste0("e", j), "OBJECTS INCLUDE:", paste(ls(pattern=pat2, envir=environment()), collapse=", "), collapse=" "))
		rm(list=d.names, envir=environment())
		gc()
        loc.grp <- d.list[[1]]$LocGroup[1]
		loc <- d.list[[1]]$Location[1]
		nr <- sapply(d.list, nrow)
        d.tmp <- rbindlist(d.list)
        denFun2 <- function(..., n=n.samples) denFun(..., n=n)
		d.tmp %>% filter(VegID %in% 1:6) %>% mutate(Var="Age", Vegetation=VegID) %>%
            select(Location, Var, Vegetation, Year, Val) %>%
            group_by(Location, Var, Vegetation, Year) %>% summarise(Vals_Probs=denFun2(Val)) -> d.tmp
        veg.id <- sort(unique(d.tmp$Vegetation))
        d.tmp[, Model:=swapModelName(mod.scen[1])]
		d.tmp[, Scenario:=swapScenarioName(mod.scen[2])]
		d.tmp[, Phase:=getPhase(mod.scen[2])]
        d.tmp[, vals.ind:=rep(rep(c(T,F), each=n.samples), length=nrow(d.tmp))]
        d.tmp[, dcastIsDumb:=rep(1:n.samples, length=nrow(d.tmp))]
        d.tmp <- data.table(dcast(d.tmp, Phase + Scenario + Model + Location + Var + Vegetation + Year + dcastIsDumb ~ vals.ind, value.var="Vals_Probs"))
        d.tmp[, dcastIsDumb:=NULL]
        d.tmp.names <- names(d.tmp)
        d.tmp.names[match(c("TRUE", "FALSE"), d.tmp.names)] <- c("Val", "Prob")
        setnames(d.tmp, d.tmp.names)
        setcolorder(d.tmp, c(datnames, "Prob"))
        d.tmp[, Vegetation:=veg.labels[veg.id][d.tmp[, Vegetation]]]
        d.tmp[, Val:=as.integer(round(Val))]
        setorder(d.tmp, Year, Vegetation)
        btfun2 <- function(..., n=n.samples, n.b=n.boot) btfun_dt(..., n.samples=n, n.boot=n.b)
        d.tmp %>% select(Phase, Scenario, Model, Location, Var, Vegetation, Year, Val, Prob) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Val=btfun2(Val, Prob)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>%
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1), Min=round(min(Val)),
                Pct_05=round(quantile(Val, 0.05)), Pct_10=round(quantile(Val, 0.10)), Pct_25=round(quantile(Val, 0.25)), Pct_50=round(quantile(Val, 0.50)),
                Pct_75=round(quantile(Val, 0.75)), Pct_90=round(quantile(Val, 0.90)), Pct_95=round(quantile(Val, 0.95)), Max=round(max(Val))) -> stats.tmp
        dat.list[[i]] <- d.tmp
        stat.list[[i]] <- stats.tmp
		print(i)
	}
	d.alf.age <- rbindlist(dat.list)
    region.dat <- rbindlist(stat.list)
	rm(dat.list, stat.list, d.tmp, stats.tmp)
    gc()
    d.alf.age[, Scenario:=factor(Scenario, levels=scen.levels)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
    save(d.alf.age, file=paste0(samplesDir, "/", "vegetationAge.RData"))
	gc()
	save(region.dat, file=file.path(statsDir, "stats_age.RData"))
    return(NULL)
}

# @knitr get_fireStats_fireDensities
get_fireStats_fireDensities <- function(j, inDir, n.samples=1000, n.boot=10000, datnames, scen.levels){
	files <- list.files(inDir, full=T, pattern="^abfc__.*.RData$")
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
    if(loc %in% c("Manitoba", "Saskatchewan")) return(NULL)
	files <- files[which(files.locs %in% loc)]
    dat.list <- stat.list <- vector("list", length(files))
	for(i in 1:length(files)){
		load(files[i], envir=environment())
		d.tmp <- data.table(d.abfc, Vegetation="All")
		rm(d.abfc)
		gc()
		loc.grp <- d.tmp$LocGroup[1]
		loc <- d.tmp$Location[1]
        denFun2 <- function(..., n=n.samples, isfire=TRUE) denFun(..., n=n, fire=isfire)
		d.tmp %>% select(which(names(d.tmp) %in% datnames)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Vals_Probs=denFun2(Val)) -> d.tmp
        d.tmp[, vals.ind:=rep(rep(c(T,F), each=n.samples), length=nrow(d.tmp))]
        d.tmp[, dcastIsDumb:=rep(1:n.samples, length=nrow(d.tmp))]
        d.tmp <- data.table(dcast(d.tmp, Phase + Scenario + Model + Location + Var + Vegetation + Year + dcastIsDumb ~ vals.ind, value.var="Vals_Probs"))
        d.tmp[, dcastIsDumb:=NULL]
        d.tmp.names <- names(d.tmp)
        d.tmp.names[match(c("TRUE", "FALSE"), d.tmp.names)] <- c("Val", "Prob")
        setnames(d.tmp, d.tmp.names)
        setcolorder(d.tmp, c(datnames, "Prob"))
        setorder(d.tmp, Var, Year)
        btfun2 <- function(..., n=n.samples, n.b=n.boot) btfun_dt(..., n.samples=n, n.boot=n.b)
        d.tmp %>% select(Phase, Scenario, Model, Location, Var, Vegetation, Year, Val, Prob) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Val=btfun2(Val, Prob)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>%
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1), Min=round(min(Val)),
                Pct_05=round(quantile(Val, 0.05)), Pct_10=round(quantile(Val, 0.10)), Pct_25=round(quantile(Val, 0.25)), Pct_50=round(quantile(Val, 0.50)),
                Pct_75=round(quantile(Val, 0.75)), Pct_90=round(quantile(Val, 0.90)), Pct_95=round(quantile(Val, 0.95)), Max=round(max(Val))) -> stats.tmp
        dat.list[[i]] <- d.tmp
        stat.list[[i]] <- stats.tmp
        print(i)
	}
	d.alf.fire <- rbindlist(dat.list)
    region.dat <- rbindlist(stat.list)
	rm(dat.list, stat.list, d.tmp, stats.tmp)
    gc()
	d.alf.fire[, Scenario:=factor(Scenario, levels=scen.levels)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
    d.alf.ba <- d.alf.fire[Var=="Burn Area",]
    d.alf.fc <- d.alf.fire[Var=="Fire Count",]
    rm(d.alf.fire)
    gc()
    save(d.alf.ba, file=paste0(samplesDir, "/", "burnArea.RData"))
    save(d.alf.fc, file=paste0(samplesDir, "/", "fireCount.RData"))
	save(region.dat, file=file.path(statsDir, "stats_fire.RData"))
    return(NULL)
}

# @knitr get_vegStats_vegDensities
get_vegStats_vegDensities <- function(j, inDir, n.samples=1000, n.boot=10000, datnames, veg.labels, scen.levels){
	files <- list.files(inDir, full=T, pattern="^veg__.*.RData$")
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
    if(loc %in% c("Manitoba", "Saskatchewan")) return(NULL)
	files <- files[which(files.locs %in% loc)]
    dat.list <- stat.list <- vector("list", length(files))
	for(i in 1:length(files)){
		load(files[i], envir=environment())
        d.veg %>% mutate(Var="Veg Area", Vegetation=VegID) -> d.tmp
        rm(d.veg)
		gc()
		loc.grp <- d.tmp$LocGroup[1]
		loc <- d.tmp$Location[1]
        veg.id <- sort(unique(d.tmp$Vegetation))
        denFun2 <- function(..., n=n.samples, isfire=TRUE) denFun(..., n=n, fire=isfire)
		d.tmp %>% select(which(names(d.tmp) %in% datnames)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Vals_Probs=denFun2(Val)) -> d.tmp
        d.tmp[, vals.ind:=rep(rep(c(T,F), each=n.samples), length=nrow(d.tmp))]
        d.tmp[, dcastIsDumb:=rep(1:n.samples, length=nrow(d.tmp))]
        d.tmp <- data.table(dcast(d.tmp, Phase + Scenario + Model + Location + Var + Vegetation + Year + dcastIsDumb ~ vals.ind, value.var="Vals_Probs"))
        d.tmp[, dcastIsDumb:=NULL]
        d.tmp.names <- names(d.tmp)
        d.tmp.names[match(c("TRUE", "FALSE"), d.tmp.names)] <- c("Val", "Prob")
        setnames(d.tmp, d.tmp.names)
        setcolorder(d.tmp, c(datnames, "Prob"))
        d.tmp[, Vegetation:=veg.labels[veg.id][d.tmp[, Vegetation]]]
		d.tmp[, Val:=as.integer(round(Val))]
        setorder(d.tmp, Var, Year)
        btfun2 <- function(..., n=n.samples, n.b=n.boot) btfun_dt(..., n.samples=n, n.boot=n.b)
        d.tmp %>% select(Phase, Scenario, Model, Location, Var, Vegetation, Year, Val, Prob) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Val=btfun2(Val, Prob)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>%
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1), Min=round(min(Val)),
                Pct_05=round(quantile(Val, 0.05)), Pct_10=round(quantile(Val, 0.10)), Pct_25=round(quantile(Val, 0.25)), Pct_50=round(quantile(Val, 0.50)),
                Pct_75=round(quantile(Val, 0.75)), Pct_90=round(quantile(Val, 0.90)), Pct_95=round(quantile(Val, 0.95)), Max=round(max(Val))) -> stats.tmp
        dat.list[[i]] <- d.tmp
        stat.list[[i]] <- stats.tmp
        print(i)
	}
	d.alf.veg <- rbindlist(dat.list)
    region.dat <- rbindlist(stat.list)
	rm(dat.list, stat.list, d.tmp, stats.tmp)
    gc()
	d.alf.veg[, Scenario:=factor(Scenario, levels=scen.levels)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
    save(d.alf.veg, file=file.path(samplesDir, "vegetatedArea.RData"))
	gc()
	save(region.dat, file=file.path(statsDir, "stats_veg.RData"))
    return(NULL)
}

# @knitr get_fsvStats_fsvDensities
get_fsvStats_fsvDensities <- function(j, inDir, n.samples=1000, n.boot=10000, datnames, veg.labels, scen.levels){
	files <- list.files(inDir, full=T, pattern="^fsv__.*.RData$")
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
    if(loc %in% c("Manitoba", "Saskatchewan")) return(NULL)
	files <- files[which(files.locs %in% loc)]
    dat.list <- stat.list <- vector("list", length(files))
	for(i in 1:length(files)){
		load(files[i], envir=environment())
        loc.grp <- d.fsv$LocGroup[1]
		loc <- d.fsv$Location[1]
		d.fsv %>% filter(VegID %in% 1:6) %>% mutate(Var="Fire Size", Vegetation=VegID) -> d.tmp
		rm(d.fsv)
		gc()
        veg.id <- sort(unique(d.tmp$Vegetation))
        denFun2 <- function(..., n=n.samples, isfire=TRUE) denFun(..., n=n, fire=isfire)
		d.tmp %>% select(which(names(d.tmp) %in% datnames)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Vals_Probs=denFun2(Val)) -> d.tmp
        d.tmp[, vals.ind:=rep(rep(c(T,F), each=n.samples), length=nrow(d.tmp))]
        d.tmp[, dcastIsDumb:=rep(1:n.samples, length=nrow(d.tmp))]
        d.tmp <- data.table(dcast(d.tmp, Phase + Scenario + Model + Location + Var + Vegetation + Year + dcastIsDumb ~ vals.ind, value.var="Vals_Probs"))
        d.tmp[, dcastIsDumb:=NULL]
        d.tmp.names <- names(d.tmp)
        d.tmp.names[match(c("TRUE", "FALSE"), d.tmp.names)] <- c("Val", "Prob")
        setnames(d.tmp, d.tmp.names)
        setcolorder(d.tmp, c(datnames, "Prob"))
        d.tmp[, Vegetation:=veg.labels[veg.id][d.tmp[, Vegetation]]]
		d.tmp[, Val:=as.integer(round(Val))]
        setorder(d.tmp, Var, Year)
        btfun2 <- function(..., n=n.samples, n.b=n.boot) btfun_dt(..., n.samples=n, n.boot=n.b)
        d.tmp %>% select(Phase, Scenario, Model, Location, Var, Vegetation, Year, Val, Prob) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% summarise(Val=btfun2(Val, Prob)) %>%
            group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>%
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1), Min=round(min(Val)),
                Pct_05=round(quantile(Val, 0.05)), Pct_10=round(quantile(Val, 0.10)), Pct_25=round(quantile(Val, 0.25)), Pct_50=round(quantile(Val, 0.50)),
                Pct_75=round(quantile(Val, 0.75)), Pct_90=round(quantile(Val, 0.90)), Pct_95=round(quantile(Val, 0.95)), Max=round(max(Val))) -> stats.tmp
        dat.list[[i]] <- d.tmp
        stat.list[[i]] <- stats.tmp
        print(i)
	}
	d.alf.fs <- rbindlist(dat.list)
    region.dat <- rbindlist(stat.list)
	rm(dat.list, stat.list, d.tmp, stats.tmp)
    gc()
	d.alf.fs$Scenario <- factor(d.alf.fs$Scenario, levels=scen.levels)
	region.dat$Scenario <- factor(region.dat$Scenario, levels=scen.levels)
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
    save(d.alf.fs, file=file.path(samplesDir, "fireSize.RData"))
	gc()
	save(region.dat, file=file.path(statsDir, "stats_fsv.RData"))
    return(NULL)
}

# Processing
# @knitr proc_fire
if("abfc" %in% variable) out.fire <- mclapply(1:n.regions, get_fireStats_fireDensities, inDir=abfcDir, n.samples=n.samples, datnames=datnames, scen.levels=scen.levels, mc.cores=n.cores)
# @knitr proc_veg
if("veg" %in% variable) out.veg <- mclapply(1:n.regions, get_vegStats_vegDensities, inDir=vegDir, n.samples=n.samples, datnames=datnames, veg.labels=veg.labels, scen.levels=scen.levels, mc.cores=n.cores)
# @knitr proc_fsv
if("fsv" %in% variable) out.fsv <- mclapply(1:n.regions, get_fsvStats_fsvDensities, inDir=fsvDir, n.samples=n.samples, datnames=datnames, veg.labels=veg.labels, scen.levels=scen.levels, mc.cores=n.cores)
# @knitr proc_age
if("age" %in% variable) out.age <- mclapply(1:n.regions, get_ageStats_ageDensities, dirs=ageDirs, n.samples=n.samples, n.samples.in=n.samples.in, datnames=datnames, veg.labels=veg.labels, scen.levels=scen.levels, mc.cores=n.cores)
