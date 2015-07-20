# @knitr setup
library(parallel)
library(reshape2)
library(data.table)
library(dplyr)
n.regions <- 12 # known, fixed
n.cores <- 12 # of 32 available
n.samples <- 1000 # new density estimation
n.samples.in <- 100 # known, based on upstream settings for vegetation age
scen.levels <- c("SRES B1", "SRES A1B", "SRES A2", "RCP 4.5", "RCP 6.0", "RCP 8.5")
veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
ageDirs <- list.files("/big_scratch/mfleonawicz/Rmpi/outputs/ageDensities", full=T)
abfcDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire/abfc"
fsvDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/fire/fsv"
vegDir <- "/big_scratch/mfleonawicz/Rmpi/outputs/veg"

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

# additional support functions
getSamples <- function(d, v) d[d[,2]==v, 1]

getStats <- function(x, seq.q){
	if(all(is.na(x))) return(rep(NA, 2 + length(seq.q)))
	x1 <- round(mean(x, na.rm=TRUE))
	x2 <- round(sd(x, na.rm=TRUE), 1)
	x3 <- round(quantile(x, probs=seq.q, na.rm=TRUE))
	c(x1, x2, x3)
}

getStatsList <- function(x, seq.q){
	if(all(is.na(x))) return(rep(NA, 2 + length(seq.q)))
	x1 <- as.list(round(mean(x, na.rm=TRUE)))
	x2 <- as.list(round(sd(x, na.rm=TRUE), 1))
	x3 <- as.list(round(quantile(x, probs=seq.q, na.rm=TRUE)))
	c(x1, x2, x3)
}

# @knitr get_AgeDensities
# Primary processing functions
get_AgeDensities <- function(j, dirs, n.samples=100, n.samples.in, n.boot=10000, datnames, multiplier, veg.labels, scen.levels){
    qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
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
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1),
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
	#d.alf.age[, Decade:=paste0(substr(Year,1,3),0)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	#region.dat[, Decade:=paste0(substr(Year,1,3),0)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
	#d.alf.age[, Prob:=round(Prob,22)*multiplier[2]]
	#r.alf.age <- unlist(d.alf.age[, samples.index, with=FALSE])
	#if(j==1) x <- d.alf.age else x <- NULL
	#if(j==1) y <- region.dat else y <- NULL
	#save(r.alf.age, file=paste0(samplesDir, "/", "vegetationAge.RData"))
    save(d.alf.age, file=paste0(samplesDir, "/", "vegetationAge.RData"))
	gc()
	#region.dat <- unlist(region.dat[,stats.index])
	#names(region.dat) <- NULL
	save(region.dat, file=file.path(statsDir, "stats_age.RData"))
	#return(list(alf.ageSamples.df=x))#, alf.fireStats.df=y))
    return(NULL)
}

# @knitr get_FireStats_FireDensities
get_FireStats_FireDensities <- function(j, inDir, n.samples=1000, n.boot=10000, datnames, multiplier, scen.levels){
    qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
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
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1),
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
	#d.alf.fire[, Decade:=paste0(substr(Year,1,3),0)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	#region.dat[, Decade:=paste0(substr(Year,1,3),0)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
	#d.alf.fire[, Prob:=round(Prob,22)*multiplier[2]]
    #r.alf.fire <- d.alf.fire[Var=="Burn Area", c(Val, Prob)]
	#if(j==1) x <- d.alf.fire[Var=="Burn Area",] else x <- NULL
	#if(j==1) y <- region.dat else y <- NULL
	#save(r.alf.fire, file=paste0(samplesDir, "/", "burnArea.RData"))
	#gc()
    #r.alf.fire <- d.alf.fire[Var=="Fire Count", c(Val, Prob)]
	#save(r.alf.fire, file=file.path(samplesDir, "fireCount.RData"))
    d.alf.ba <- d.alf.fire[Var=="Burn Area",]
    d.alf.fc <- d.alf.fire[Var=="Fire Count",]
    rm(d.alf.fire)
    gc()
    save(d.alf.ba, file=paste0(samplesDir, "/", "burnArea.RData"))
    save(d.alf.fc, file=paste0(samplesDir, "/", "fireCount.RData"))
    #region.dat <- unlist(region.dat[, which(names(region.dat) %in% statnames), with=FALSE])
	save(region.dat, file=file.path(statsDir, "stats_fire.RData"))
	#return(list(alf.fireSamples.df=x, alf.fireStats.df=y))
    return(NULL)
}

# @knitr get_VegStats_VegDensities
get_VegStats_VegDensities <- function(j, inDir, n.samples=1000, n.boot=10000, datnames, multiplier, veg.labels, scen.levels){
    qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
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
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1),
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
	#d.alf.veg[, Decade:=paste0(substr(Year,1,3),0)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	#region.dat[, Decade:=paste0(substr(Year,1,3),0)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
	#d.alf.veg[, Prob:=round(Prob,22)*multiplier[2]]
    #r.alf.veg <- d.alf.veg[, c(Val, Prob)]
	#if(j==1) x <- d.alf.veg else x <- NULL
	#if(j==1) y <- region.dat else y <- NULL
	#save(r.alf.veg, file=file.path(samplesDir, "vegetatedArea.RData"))
    save(d.alf.veg, file=file.path(samplesDir, "vegetatedArea.RData"))
	gc()
	#region.dat <- unlist(region.dat[, which(names(region.dat) %in% statnames), with=FALSE])
	save(region.dat, file=file.path(statsDir, "stats_veg.RData"))
	#return(list(alf.vegSamples.df=x, alf.vegStats.df=y))
    return(NULL)
}

# @knitr get_fsvStats_fsvDensities
get_fsvStats_fsvDensities <- function(j, inDir, n.samples=1000, n.boot=10000, datnames, multiplier, veg.labels, scen.levels){
    qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
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
            summarise(Mean=round(mean(Val)), SD=round(sd(Val),1),
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
	#d.alf.fs$Decade <- paste0(substr(d.alf.fs$Year,1,3),0)
	region.dat$Scenario <- factor(region.dat$Scenario, levels=scen.levels)
	#region.dat$Decade <- paste0(substr(region.dat$Year,1,3),0)
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alfresco/samples", loc.grp, loc), recursive=T, showWarnings=F)
	#d.alf.fs$Prob <- round(d.alf.fs$Prob,22)*multiplier[2]
	#r.alf.fs <- d.alf.fs[, c(Val, Prob)]
	#if(j==1) x <- d.alf.fs else x <- NULL
	#if(j==1) y <- region.dat else y <- NULL
	#save(r.alf.fs, file=file.path(samplesDir, "fireSize.RData"))
    save(d.alf.fs, file=file.path(samplesDir, "fireSize.RData"))
	gc()
	#region.dat <- unlist(region.dat[, which(names(region.dat) %in% statnames), with=FALSE])
	save(region.dat, file=file.path(statsDir, "stats_fsv.RData"))
	#return(list(alf.fsvSamples.df=x, alf.fsvStats.df=y))
    return(NULL)
}

# @knitr proc_setup
# Processing
datnames <- c("Phase", "Scenario", "Model", "Location", "Var", "Vegetation", "Year", "Val")
#agg.stat.colnames <- c("Mean", "SD", paste0("Pct_", c("05", 10, 25, 50, 75, 90, 95)), "Max")
#samples.multipliers.alf <- c(1e1, 1e22) # First not used
#sampnames <- c("Val", "Prob", "Var", "Location")

# @knitr proc_fire
system.time( out.fire <- mclapply(1:n.regions, get_FireStats_FireDensities, inDir=abfcDir, n.samples=n.samples,
    datnames=datnames, scen.levels=scen.levels, mc.cores=n.cores) )
#alf.fireStats.df <- out.fire[[1]]$alf.fireStats.df
#alf.fireSamples.df <- out.fire[[1]]$alf.fireSamples.df
#rm(out.fire)
#gc()

#nam <- names(alf.fireStats.df)
#alf.fireStats.df[, which(nam %in% c("Location", agg.stat.colnames)):=NA]
#nam <- names(alf.fireSamples.df)
#alf.fireSamples.df[, which(nam %in% sampnames):=NA]
#rm(nam)

# @knitr proc_veg
system.time( out.veg <- mclapply(1:n.regions, get_VegStats_VegDensities, inDir=vegDir, n.samples=n.samples,
    datnames=datnames, veg.labels=veg.labels, scen.levels=scen.levels, mc.cores=n.cores) )
#alf.vegStats.df <- out.veg[[1]]$alf.vegStats.df
#alf.vegSamples.df <- out.veg[[1]]$alf.vegSamples.df
#rm(out.veg)
#gc()

#nam <- names(alf.vegStats.df)
#alf.vegStats.df[, which(nam %in% c("Location", agg.stat.colnames)):=NA]
#nam <- names(alf.vegSamples.df)
#alf.vegSamples.df[, which(nam %in% sampnames):=NA]
#rm(nam)

# @knitr proc_fsv
system.time( out.fsv <- mclapply(1:n.regions, get_fsvStats_fsvDensities, inDir=fsvDir, n.samples=n.samples,
    datnames=datnames, veg.labels=veg.labels, scen.levels=scen.levels, mc.cores=n.cores) )
#alf.fsvStats.df <- out.fsv[[1]]$alf.fsvStats.df
#alf.fsvSamples.df <- out.fsv[[1]]$alf.fsvSamples.df
#rm(out.fsv)
#gc()

#nam <- names(alf.fsvStats.df)
#alf.fsvStats.df[, which(nam %in% c("Location", agg.stat.colnames)):=NA]
#nam <- names(alf.fsvSamples.df)
#alf.fsvSamples.df[, which(nam %in% sampnames):=NA]
#rm(nam)

# @knitr proc_age
system.time( out.age <- mclapply(1:n.regions, get_AgeDensities, dirs=ageDirs, n.samples=n.samples, n.samples.in=n.samples.in,
    datnames=datnames, veg.labels=veg.labels, scen.levels=scen.levels, mc.cores=n.cores) )
#alf.ageStats.df <- out.age[[1]]$alf.ageStats.df
#alf.ageSamples.df <- out.age[[1]]$alf.ageSamples.df
#rm(out.age)
#gc()

#nam <- names(alf.ageSamples.df)
#alf.ageSamples.df[, which(nam %in% sampnames):=NA]
#rm(nam)

# @knitr save_metadata
# Remove unwanted objects, load metadata workspace, save along with age metadata
#rm(ageDirs, agg.stat.colnames, btfun, denFun, abfcDir, fsvDir, get_AgeDensities, get_FireStats_FireDensities, get_fsvStats_fsvDensities, get_VegStats_VegDensities, getPhase,
#getSamples, getStats, n.cores, n.regions, n.samples, samples.columns, scen.levels, stats.columns, swapModelName, swapScenarioName, vegDir)
# Create a backup of the meta.RData file, load and save over that.
# Inclusion of ALFRESCO output in QAQC R Shiny app is under early development.
# All but one shiny-apps development repo branch should not include this version of meta.RData in the cmip3_cmip5 app
#load("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/meta_backup.RData")
#obj.keep <- c("alf.fireStats", "alf.fireSamples.df", "alf.vegStats.df", "alf.vegSamples", "alf.fsvStats.df", "alf.fsvSamples.df", "alf.ageSamples.df",
#    "datnames", "agg.stat.colnames", "samples.multipliers", "veg.labels", "scen.levels")
#save(list=obj.keep, file="/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/alf_meta.RData")
