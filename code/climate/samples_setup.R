# @knitr setup
setwd("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/regional/samples")

library(parallel)
library(reshape2)
library(data.table)

files <- list.files(pattern="regions_samples.RData$")
files <- files[substr(files, 1, 3) != "CRU"]
lab <- do.call(rbind, strsplit(files, "_"))[,1:4]
scen.levels <- c("SRES B1","SRES A1B","SRES A2","RCP 4.5","RCP 6.0","RCP 8.5")

swapVar <- function(x) switch(x, pr="Precipitation", tas="Temperature")
swapScen <- function(x) switch(x, sresb1="SRES B1", sresa1b="SRES A1B", sresa2="SRES A2", rcp45="RCP 4.5", rcp60="RCP 6.0", rcp85="RCP 8.5")
swapMod <- function(x){
	switch(x,
		"cccma-cgcm3-1-t47"="CCCMAcgcm31", "gfdl-cm2-1"="GFDLcm21", "miroc3-2-medres"="MIROC32m", "mpi-echam5"="MPIecham5", "ukmo-hadcm3"="ukmoHADcm3",
		"CCSM4"="CCSM4", "GFDL-CM3"="GFDLcm3", "GISS-E2-R"="GISSe2-r", "IPSL-CM5A-LR"="IPSLcm5a-lr", "MRI-CGCM3"="MRIcgcm3"
	)
}

lab[,2:4] <- cbind(sapply(lab[,2], swapScen), sapply(lab[,3], swapMod), sapply(lab[,4], swapVar))
rm(swapVar, swapScen, swapMod)

# @knitr func_organize_save
f <- function(i, d, n, grp.names, s.names, scen.levels, outfile, multiplier=c(1,1), decimals=NULL){
	AR <- strsplit(outfile, "_")[[1]][1]
	scen <- strsplit(outfile, "_")[[1]][2]
	outfile.h <- gsub(scen, "Hist", outfile)
	yr1.projected <- if(AR=="AR4") 2001 else 2006
	
	rsd1 <- subset(d, Location==s.names[i])
	rsd1 <- melt(rsd1, id.vars=names(rsd1)[1:5], measure.vars=names(rsd1)[-c(1:5)], variable.name="Time", value.name="Vals_Probs")
	rsd1$vals.ind <- rep(rep(c(T,F), each=n), length=nrow(rsd1))
	rsd1$dcastIsDumb <- rep(1:n, length=nrow(rsd1))
	rsd1 <- dcast(rsd1, Phase + Scenario + Model + Var + Location + Time + dcastIsDumb ~ vals.ind, value.var="Vals_Probs")
	names(rsd1)[9:8] <- c("Val", "Prob")
	rsd1$Year <- substr(as.character(rsd1$Time), 1, 4)
	rsd1$Month <- substr(as.character(rsd1$Time), 6, 8)
	rsd1$Month <- factor(rsd1$Month, levels=month.abb)
	rsd1$Scenario <- factor(rsd1$Scenario, levels=scen.levels)
	rsd1$Decade <- paste0(substr(rsd1$Year,1,3),0)
	rsd1 <- rsd1[c(1:5,9,8,10:12)]
	grp <- rep(names(grp.names), times=sapply(grp.names, length))[i]
	dir.create(outDir <- file.path("../../final/region_files_GCM/samples", grp, s.names[i], "climate"), recursive=T, showWarnings=F)
	if(!is.null(decimals)) rsd1$Val <- round(rsd1$Val, decimals[1])*multiplier[1]
	if(!is.null(decimals)) rsd1$Prob <- round(rsd1$Prob, decimals[2])*multiplier[2]
	
	rsd.h <- subset(rsd1, Year < yr1.projected)
	rsd <- subset(rsd1, Year >= yr1.projected)
	rownames(rsd) <- rownames(rsd.h) <- NULL
	rsd.h <- data.table(rsd.h)
	rsd <- data.table(rsd)
	if(!any(is.na(rsd.h$Val))) save(rsd.h, file=file.path(outDir, outfile.h))
	save(rsd, file=file.path(outDir, outfile))
	print(i)
}

# @knitr proc
if(file.exists("../../final/meta.RData")) load("../../final/meta.RData")
samples.multipliers <- c(1e1, 1e16) # see function arguments 'multiplier' and 'decimals'

for(i in 1:length(files)){
	outfile <- paste0(gsub(" ", "", paste(lab[i,], collapse="_")), ".RData")
	load(files[i])
	d <- rbindlist(samples.out)
	d <- data.frame(Phase=lab[i,1], Scenario=lab[i,2], Model=lab[i,3], Var=lab[i,4], d, stringsAsFactors=F)
	names(d)[6:ncol(d)] <- gsub("X", "", names(d)[6:ncol(d)])
	mclapply(1:length(samples.names), f, d=d, n=n.samples, grp.names=region.names.out, s.names=samples.names, scen.levels=scen.levels, outfile=outfile, multiplier=samples.multipliers, decimals=c(1, 16), mc.cores=32)
	print(i)
}

rm(lab, samples.out, i, files, d, outfile, f)
save.image("../../final/meta.RData")
