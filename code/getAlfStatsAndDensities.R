# @knitr setup
library(parallel)
library(reshape2)
library(data.table)
library(dplyr)
n.regions <- 53 # known, fixed
n.cores <- 27 # of 32 available
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
btfun <- function(p, n.samples=length(p)/2, n.boot=10000, interp=FALSE, n.interp=1000, ...){
	if(!length(p)) return(p)
	if(all(is.na(p))) return(rep(NA, n.boot))
	p <- list(x=p[1:n.samples], y=p[(n.samples+1):(2*n.samples)])
	if(interp && length(unique(p[1:n.samples])) > 1) p <- approx(p$x, p$y, n=n.interp)
	p <- round(sample(p$x, n.boot, prob=p$y, rep=T), ...)
	p
}

# @knitr functions_not_in_use
#library(compiler)
getSamples <- function(d, v) d[d[,2]==v, 1]
#getSamples_cmp <- cmpfun(getSamples) # not in use

# not in use
library(Rcpp)
cppFunction(
'List lapply_C(List input, NumericMatrix data, Function f) {
  int n = input.size();
  List out(n);

  for(int i = 0; i < n; i++) {
    out[i] = f(_["v"] = input[i], _["d"] = data);
  }

  return out;
}'
)

getStats <- function(x, seq.q){
	if(all(is.na(x))) return(rep(NA, 2 + length(seq.q)))
	x1 <- round(mean(x, na.rm=TRUE))
	x2 <- round(sd(x, na.rm=TRUE), 1)
	x3 <- round(quantile(x, probs=seq.q, na.rm=TRUE))
	c(x1, x2, x3)
}

# @knitr get_AgeDensities
# Primary processing functions
get_AgeDensities <- function(j, dirs, n.samples=100, n.samples.in, samples.index, multiplier, veg.labels, scen.levels){
	pat <- paste0("^age__.*.rep.*.RData$")
	pat2 <- substr(pat, 1, nchar(pat) - 9)
	files.list <- lapply(1:length(dirs), function(i, dirs, ...) list.files(dirs[i], ...), dirs=dirs, full=T, pattern=pat)
	for(i in 1:length(files.list)){
		files.locs <- sapply(strsplit(files.list[[i]], "__"), "[", 2)
		locs <- unique(files.locs)
		loc <- locs[j]
		files.list[[i]] <- files.list[[i]][which(files.locs %in% loc)]
	}
	for(i in 1:length(dirs)){
		modname <- basename(dirs[i])
		# Load all reps for one location, all model and scenarios, into local location-specific environment
		print(system.time( for(p in 1:length(files.list[[i]])) { load(files.list[[i]][p], envir=environment()); print(paste("Local environment: loading workspace", p, "of", length(files.list[[i]]))) } ))
		d.names <- ls(pattern=pat2, envir=environment())
		d.list <- mget(d.names, envir=environment())
		print(paste("EVIRONMENT", paste0("e", j), "OBJECTS INCLUDE:", paste(ls(pattern=pat2, envir=environment()), collapse=", "), collapse=" "))
		rm(list=d.names, envir=environment())
		gc()
		all.years <- unique(d.list[[1]]$Year)
		loc.grp <- d.list[[1]]$LocGroup[1]
		loc <- d.list[[1]]$Location[1]
		nr <- sapply(d.list, nrow)
		year.veg.list <- lapply(d.list, function(x) as.numeric(unique(paste0(x$Year, ".", x$VegID))))
		year.veg.uni <- sort(unique(unlist(year.veg.list)))
		a <- sapply(year.veg.list, function(x) all(paste0(all.years, ".", 1:8) %in% x))
		dl <- length(d.list)
		for(z in 1:dl) {
			x <- as.matrix(d.list[[z]][, 3:5, with=FALSE])
			x <- cbind(x[,2], x[,3] + x[,1]/10)
			if(!a[z]){
				x.ids <- unique(x[,2])
				veg.missing <- paste0(rep(all.years, each=8), ".", 1:8)
				veg.missing <- as.numeric(veg.missing[!(veg.missing %in% x.ids)])
				x <- rbind(x, cbind(NA, rep(veg.missing, each=2*n.samples)))
				x <- x[order(x[,2]),]
				rownames(x) <- NULL
			}
			d.list[[z]] <- split(x[,1], x[,2])
			print(paste("Obtain annual vegetation samples: replicate", z, "of", dl))
		}
		d.list <- rapply(d.list, btfun, classes="numeric", how="replace", n.samples=n.samples.in, n.boot=1000, interp=TRUE, n.interp=1000)
		d.list <- unlist(d.list, recursive=FALSE)
		p <- length(d.list)
		m <- p/length(d.names)
		system.time( for(k in 1:m) { d.list[[k]] <- denFun(unlist(d.list[seq(1, p, by=m) + k - 1]), n=n.samples); print(paste("Empirical density estimation: sample", k, "of", m)) } )
		d.list <- d.list[1:m]
		gc()
		d.tmp <- data.frame(Var="Age", Location=loc, Vegetation=rep(1:8, each=2*n.samples), Vals_Probs=unlist(d.list), Year=rep(all.years, each=8*2*n.samples), stringsAsFactors=FALSE)
		rm(d.list)
		gc()
		ind.val <- as.numeric(mapply("+", seq(1, nrow(d.tmp), by=2*n.samples), MoreArgs=list(0:(n.samples-1))))
		d.tmp$Vals_Probs[ind.val] <- round(d.tmp$Vals_Probs[ind.val])
		mod.scen <- unlist(strsplit(basename(dirs[i]), "\\."))
        d.tmp <- data.table(d.tmp)
		d.tmp[, Model:=swapModelName(mod.scen[1])]
		d.tmp[, Scenario:=swapScenarioName(mod.scen[2])]
		d.tmp[, Phase:=getPhase(mod.scen[2])]
		setcolorder(d.tmp, names(d.tmp)[c(8:6,1:5)])
		if(i==1) d.alf.age <- d.tmp else d.alf.age <- rbindlist(list(d.alf.age, d.tmp))
		print(i)
	}
	
	rm(d.tmp)
	d.alf.age[, vals.ind:=rep(rep(c(T,F), each=n.samples), length=nrow(d.alf.age))]
	d.alf.age[, dcastIsDumb:=rep(1:n.samples, length=nrow(d.alf.age))]
	d.alf.age <- data.table(dcast(d.alf.age, Phase + Scenario + Model + Var + Location + Vegetation + Year + dcastIsDumb ~ vals.ind, value.var="Vals_Probs"))
    d.alf.names <- names(d.alf.age)
    d.alf.names[10:9] <- c("Val", "Prob")
	setnames(d.alf.age, d.alf.names)
	d.alf.age[, Scenario:=factor(Scenario, levels=scen.levels)]
	d.alf.age[, Decade:=paste0(substr(Year,1,3),0)]
	d.alf.age <- d.alf.age[order(paste(Year, Vegetation)), c(1:5,10:9, 6:7, 11), with=FALSE]
    d.alf.age[, Vegetation:=veg.labels[d.alf.age[, Vegetation]]]
	#dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/samples", loc.grp, loc), recursive=T, showWarnings=F)
	d.alf.age[, Val:=as.integer(round(Val))]
	d.alf.age[, Prob:=round(Prob,16)*multiplier[2]] # round to 16 seems decent for now
	r.alf.age <- unlist(d.alf.age[, samples.index, with=FALSE])
	names(r.alf.age) <- NULL
	if(j==1) x <- d.alf.age else x <- NULL
	#if(j==1) y <- region.dat else y <- NULL
	save(r.alf.age, file=paste0(samplesDir, "/", "vegetationAge.RData"))
	gc()
	#region.dat <- unlist(region.dat[,stats.index])
	#names(region.dat) <- NULL
	#save(region.dat, file=file.path(statsDir, "stats_age.RData"))
	return(list(alf.ageSamples.df=x))#, alf.fireStats.df=y))
}

# @knitr get_FireStats_FireDensities
get_FireStats_FireDensities <- function(j, inDir, n.samples=1000, stats.index, samples.index, multiplier, scen.levels){
	files <- list.files(inDir, full=T, pattern="^abfc__.*.RData$")
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
	files <- files[which(files.locs %in% loc)]
	for(i in 1:length(files)){
		load(files[i], envir=environment())
		d.tmp <- data.table(d.abfc)
		rm(d.abfc)
		gc()
		loc.grp <- d.tmp$LocGroup[1]
		loc <- d.tmp$Location[1]
		ind <- with(d.tmp, paste(Var, Year))
		v <- unlist(tapply(d.tmp$Val, ind, denFun, n=n.samples, fire=TRUE))
		d.tmp <- d.tmp[Replicate==d.tmp$Replicate[1], c(1:3,6,5,7,8), with=FALSE]
		stats.tmp <- d.tmp
		d.tmp <- data.table(data.frame(lapply(d.tmp, rep, n.samples), stringsAsFactors=FALSE))
		d.tmp[, Val:=as.integer(round(v[rep(c(T,F), each=n.samples)]))]
		d.tmp[, Prob:=v[rep(c(F,T), each=n.samples)]]
		setcolorder(d.tmp, names(d.tmp)[c(1:6,8,7)])
		v <- tapply(v, rep(1:(length(v)/(2*n.samples)), each=2*n.samples), btfun, n.samples=n.samples, n.boot=10000, interp=TRUE)
		qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
		n.stats <- 2 + length(qtiles)
		v <- do.call(rbind, lapply(v, getStats, seq.q=qtiles))
		stats.tmp <- data.table(stats.tmp, v)
        stats.tmp[, Val:=NULL]
		setcolorder(stats.tmp, names(stats.tmp)[c(1:5,6+1:ncol(v),6)])
		if(i==1) d.alf.fire <- d.tmp else d.alf.fire <- rbind(d.alf.fire, d.tmp)
		if(i==1) region.dat <- stats.tmp else region.dat <- rbind(region.dat, stats.tmp)
	}
	rm(d.tmp, stats.tmp)
	gc()
	d.alf.fire[, Scenario:=factor(Scenario, levels=scen.levels)]
	d.alf.fire[, Decade:=paste0(substr(Year,1,3),0)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	region.dat[, Decade:=paste0(substr(Year,1,3),0)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/samples", loc.grp, loc), recursive=T, showWarnings=F)
	d.alf.fire[, Prob:=round(Prob,16)*multiplier[2]] # round to 16 seems decent for now
	r.alf.fire <- unlist(d.alf.fire[Var=="Burn Area", samples.index, with=FALSE])
	names(r.alf.fire) <- NULL
	if(j==1) x <- d.alf.fire[Var=="Burn Area",] else x <- NULL
	#if(j==1) y <- subset(region.dat, Var=="Burn Area") else y <- NULL
	#if(j==1) x <- d.alf.fire else x <- NULL
	if(j==1) y <- region.dat else y <- NULL
	save(r.alf.fire, file=paste0(samplesDir, "/", "burnArea.RData"))
	gc()
	r.alf.fire <- unlist(d.alf.fire[Var=="Fire Count", samples.index, with=FALSE])
	names(r.alf.fire) <- NULL
	save(r.alf.fire, file=file.path(samplesDir, "fireCount.RData"))
	gc()
	region.dat <- unlist(region.dat[, stats.index, with=FALSE])
	names(region.dat) <- NULL
	save(region.dat, file=file.path(statsDir, "stats_fire.RData"))
	return(list(alf.fireSamples.df=x, alf.fireStats.df=y))
}

# @knitr get_VegStats_VegDensities
get_VegStats_VegDensities <- function(j, inDir, n.samples=1000, stats.index, samples.index, multiplier, veg.labels, scen.levels){
	files <- list.files(inDir, full=T, pattern="^veg__.*.RData$")
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
	files <- files[which(files.locs %in% loc)]
	for(i in 1:length(files)){
		load(files[i], envir=environment())
		#locs <- unique(v.dat$Location)
		#d.tmp <- v.dat[v.dat$Location==locs[j],]
		d.tmp <- data.table(d.veg[1, 1:6, with=FALSE], Val=0L, Vegetation=1:8, Year=rep(unique(d.veg$Year), each=8), Replicate=rep(unique(d.veg$Replicate), each=8*length(unique(d.veg$Year))))
        d.tmp[, VegID:=NULL]
		obs.ind <- paste(d.tmp$Year, d.tmp$Vegetation, d.tmp$Replicate) %in% paste(d.veg$Year, d.veg$VegID, d.veg$Replicate)
        d.tmp <- d.tmp[obs.ind,]
		d.tmp[, Val:=d.veg$Val[obs.ind]]
		rm(d.veg)
		gc()
		loc.grp <- d.tmp$LocGroup[1]
		loc <- d.tmp$Location[1]
		ind <- with(d.tmp, paste(Year, Vegetation))
		v <- unlist(tapply(d.tmp$Val, ind, denFun, n=n.samples))
		d.tmp <- d.tmp[Replicate==d.tmp$Replicate[1], c(1:3,5:8), with=FALSE]
		stats.tmp <- d.tmp
        d.tmp <- data.table(data.frame(lapply(d.tmp, rep, n.samples), stringsAsFactors=FALSE))
		d.tmp <- d.tmp[order(paste(Year, Vegetation)),]
        stats.tmp[, Vegetation:=veg.labels[stats.tmp[, Vegetation]]]
        d.tmp[, Vegetation:=veg.labels[d.tmp[, Vegetation]]]
		d.tmp[, Val:=as.integer(round(v[rep(c(T,F), each=n.samples)]))]
		d.tmp[, Prob:=v[rep(c(F,T), each=n.samples)]]
		setcolorder(d.tmp, names(d.tmp)[c(1:5,8,6,7)])
		v <- tapply(v, rep(1:(length(v)/(2*n.samples)), each=2*n.samples), btfun, n.samples=n.samples, n.boot=10000, interp=TRUE)
		qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
		n.stats <- 2 + length(qtiles)
		v <- do.call(rbind, lapply(v, getStats, seq.q=qtiles))
		stats.tmp <- data.table(stats.tmp, v)
		setcolorder(stats.tmp, names(stats.tmp)[c(1:5,7+1:ncol(v),6,7)])
		if(i==1) d.alf.veg <- d.tmp else d.alf.veg <- rbind(d.alf.veg, d.tmp)
		if(i==1) region.dat <- stats.tmp else region.dat <- rbind(region.dat, stats.tmp)
	}
	rm(d.tmp, stats.tmp)
	gc()
	d.alf.veg[, Scenario:=factor(Scenario, levels=scen.levels)]
	d.alf.veg[, Decade:=paste0(substr(Year,1,3),0)]
	region.dat[, Scenario:=factor(Scenario, levels=scen.levels)]
	region.dat[, Decade:=paste0(substr(Year,1,3),0)]
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/samples", loc.grp, loc), recursive=T, showWarnings=F)
	d.alf.veg[, Prob:=round(Prob,16)*multiplier[2]] # round to 16 seems decent for now
	r.alf.veg <- unlist(d.alf.veg[, samples.index, with=FALSE])
	names(r.alf.veg) <- NULL
	if(j==1) x <- d.alf.veg else x <- NULL
	if(j==1) y <- region.dat else y <- NULL
	save(r.alf.veg, file=file.path(samplesDir, "vegetatedArea.RData"))
	gc()
	region.dat <- unlist(region.dat[, stats.index, with=FALSE])
	names(region.dat) <- NULL
	save(region.dat, file=file.path(statsDir, "stats_veg.RData"))
	return(list(alf.vegSamples.df=x, alf.vegStats.df=y))
}

# @knitr get_fsvStats_fsvDensities
get_fsvStats_fsvDensities <- function(j, inDir, n.samples=1000, stats.index, samples.index, multiplier, veg.labels, scen.levels){
	files <- list.files(inDir, full=T, pattern="^fsv__.*.RData$")
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
	files <- files[which(files.locs %in% loc)]
	for(i in 1:length(files)){
		load(files[i], envir=environment())
		#locs <- unique(v.dat$Location)
		#d.tmp <- v.dat[v.dat$Location==locs[j],]
		d.fsv %>% filter(VegID %in% 1:6) %>% mutate(Vegetation=VegID) %>% select(Phase, Scenario, Model, Location, Vegetation, Val, Year) -> d.tmp
		rm(d.fsv)
		gc()
		loc.grp <- d.tmp$LocGroup[1]
		loc <- d.tmp$Location[1]
		ind <- with(d.tmp, paste(Year, Vegetation))
		v <- unlist(tapply(d.tmp$Val, ind, denFun, n=n.samples))
		d.tmp %>% group_by(Phase, Scenario, Model, Location, Vegetation, Year) %>% summarise(Val=0L) -> d.tmp
		stats.tmp <- d.tmp
		d.tmp <- data.table(data.frame(lapply(d.tmp, rep, n.samples), stringsAsFactors=FALSE))
		d.tmp <- d.tmp[order(paste(Year, Vegetation)),]
        stats.tmp[, Vegetation:=veg.labels[stats.tmp[, Vegetation]]]
        d.tmp[, Vegetation:=veg.labels[d.tmp[, Vegetation]]]
        d.tmp[, Val:=as.integer(round(v[rep(c(T,F), each=n.samples)]))]
		d.tmp[, Prob:=v[rep(c(F,T), each=n.samples)]]
		setcolorder(d.tmp, names(d.tmp)[c(1:5,7,8,6)])
		v <- tapply(v, rep(1:(length(v)/(2*n.samples)), each=2*n.samples), btfun, n.samples=n.samples, n.boot=10000, interp=TRUE)
		qtiles <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)
		n.stats <- 2 + length(qtiles)
		v <- do.call(rbind, lapply(v, getStats, seq.q=qtiles))
		stats.tmp <- data.table(stats.tmp, v)[, Val:=NULL]
		setcolorder(stats.tmp, names(stats.tmp)[c(1:5,6+1:ncol(v),6)])
        #d.fsv %>% filter(Year==2004) %>% group_by(VegID) %>% summarise(Mean=mean(Val), SD=sd(Val), Pct05=quantile(Val, 0.05), Pct25=quantile(Val, 0.25),Pct50=quantile(Val, 0.5),Pct75=quantile(Val, 0.75),Pct95=quantile(Val, 0.95), Max=max(Val))
        #stats.tmp[Year==2004,]
		if(i==1) d.alf.veg <- d.tmp else d.alf.veg <- rbind(d.alf.veg, d.tmp)
		if(i==1) region.dat <- stats.tmp else region.dat <- rbind(region.dat, stats.tmp)
	}
	rm(d.tmp, stats.tmp)
	gc()
	d.alf.veg$Scenario <- factor(d.alf.veg$Scenario, levels=scen.levels)
	d.alf.veg$Decade <- paste0(substr(d.alf.veg$Year,1,3),0)
	region.dat$Scenario <- factor(region.dat$Scenario, levels=scen.levels)
	region.dat$Decade <- paste0(substr(region.dat$Year,1,3),0)
	dir.create(statsDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/region_files_GCM/samples", loc.grp, loc), recursive=T, showWarnings=F)
	d.alf.veg$Prob <- round(d.alf.veg$Prob,16)*multiplier[2] # round to 16 seems decent for now
	r.alf.veg <- unlist(d.alf.veg[,samples.index])
	names(r.alf.veg) <- NULL
	if(j==1) x <- d.alf.veg else x <- NULL
	if(j==1) y <- region.dat else y <- NULL
	save(r.alf.veg, file=file.path(samplesDir, "fireSize.RData"))
	gc()
	region.dat <- unlist(region.dat[,stats.index])
	names(region.dat) <- NULL
	save(region.dat, file=file.path(statsDir, "stats_fsv.RData"))
	return(list(alf.fsvSamples.df=x, alf.fsvStats.df=y))
}

# @knitr proc_setup
# Processing
agg.stat.colnames <- c("Mean", "SD", paste0("Pct_", c("05", 10, 25, 50, 75, 90, 95)), "Max")
stats.columns <- 5+1:length(agg.stat.colnames)
samples.columns <- 6:7
samples.multipliers.alf <- c(1e1, 1e16) # First not used

# @knitr proc_fire
system.time( out.fire <- mclapply(1:n.regions, get_FireStats_FireDensities, inDir=abfcDir, n.samples=n.samples, stats.index=stats.columns, samples.index=samples.columns, multiplier=samples.multipliers.alf, mc.cores=n.cores) )
alf.fireStats.df <- out.fire[[1]]$alf.fireStats.df
alf.fireSamples.df <- out.fire[[1]]$alf.fireSamples.df
rm(out.fire)
gc()

nam <- names(alf.fireStats.df)
alf.fireStats.df[, nam[stats.columns]:=NA]
alf.fireStats.df[, Location:=NA]
nam[stats.columns] <- agg.stat.colnames
setnames(alf.fireStats.df, nam)

nam <- names(alf.fireSamples.df)
alf.fireSamples.df[, nam[samples.columns]:=NA]
alf.fireSamples.df[, Var:=NA]
alf.fireSamples.df[, Location:=NA]
rm(nam)

# @knitr proc_veg
system.time( out.veg <- mclapply(1:n.regions, get_VegStats_VegDensities, inDir=vegDir, n.samples=n.samples, stats.index=stats.columns, samples.index=samples.columns, multiplier=samples.multipliers.alf, veg.labels=veg.labels, mc.cores=n.cores) )
alf.vegStats.df <- out.veg[[1]]$alf.vegStats.df
alf.vegSamples.df <- out.veg[[1]]$alf.vegSamples.df
rm(out.veg)
gc()

nam <- names(alf.vegStats.df)
alf.vegStats.df[, nam[stats.columns]:=NA]
alf.vegStats.df[, Location:=NA]
nam[stats.columns] <- agg.stat.colnames
setnames(alf.vegStats.df, nam)

nam <- names(alf.vegSamples.df)
alf.vegSamples.df[, nam[samples.columns]:=NA]
alf.vegSamples.df[, Var:=NA]
alf.vegSamples.df[, Location:=NA]
rm(nam)

# @knitr proc_fsv
system.time( out.fsv <- mclapply(1:n.regions, get_fsvStats_fsvDensities, inDir=fsvDir, n.samples=n.samples, stats.index=stats.columns, samples.index=samples.columns, multiplier=samples.multipliers.alf, veg.labels=veg.labels, mc.cores=n.cores) )
alf.fsvStats.df <- out.veg[[1]]$alf.fsvStats.df
alf.fsvSamples.df <- out.veg[[1]]$alf.fsvSamples.df
rm(out.veg)
gc()

alf.fsvStats.df[,stats.columns] <- NA
alf.fsvStats.df$Location <- NA
names(alf.fsvStats.df)[stats.columns] <- agg.stat.colnames

alf.fsvSamples.df[,samples.columns] <- NA
alf.fsvSamples.df$Var <- NA
alf.fsvSamples.df$Location <- NA

# @knitr proc_age
system.time( out.age <- mclapply(1:n.regions, get_AgeDensities, dirs=ageDirs, n.samples=n.samples, n.samples.in=n.samples.in, samples.index=samples.columns, multiplier=samples.multipliers.alf, veg.labels=veg.labels, mc.cores=n.cores) )
#alf.ageStats.df <- out.age[[1]]$alf.ageStats.df
alf.ageSamples.df <- out.age[[1]]$alf.ageSamples.df
rm(out.age)
gc()

#alf.ageStats.df[,stats.columns] <- NA
#alf.ageStats.df$Location <- NA
#names(alf.ageStats.df)[stats.columns] <- agg.stat.colnames

alf.ageSamples.df[, Val:=NA]
alf.ageSamples.df[, Prob:=NA]
alf.ageSamples.df[, Var:=NA]
alf.ageSamples.df[, Location:=NA]

# @knitr save_metadata
# Remove unwanted objects, load metadata workspace, save along with age metadata
rm(ageDirs, agg.stat.colnames, btfun, denFun, abfcDir, fsvDir, get_AgeDensities, get_FireStats_FireDensities, get_fsvStats_fsvDensities, get_VegStats_VegDensities, getPhase,
getSamples, getStats, lapply_C, n.cores, n.regions, n.samples, samples.columns, scen.levels, stats.columns, swapModelName, swapScenarioName, vegDir)
# Create a backup of the meta.RData file, load and save over that.
# Inclusion of ALFRESCO output in QAQC R Shiny app is under early development.
# All but one shiny-apps development repo branch should not include this version of meta.RData in the cmip3_cmip5 app
load("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/meta_backup.RData")
save.image("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final/meta_backup.RData")
