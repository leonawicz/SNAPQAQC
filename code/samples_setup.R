# @knitr setup
setwd("/workspace/UA/mfleonawicz/leonawicz/Projects/active/AR4_AR5_comparisons/data/regional/samples")

library(data.table)

files <- list.files(pattern="regions_samples.RData$")
files <- files[substr(files, 1, 3) != "CRU"]

models <- list(
	c("CCCMAcgcm31","CCSM4"),
	c("GFDLcm21","GFDLcm3"),
	c("MIROC32m","GISSe2-r"),
	c("MPIecham5","IPSLcm5a-lr"),
	c("ukmoHADcm3","MRIcgcm3")
)

# @knitr load
dlist <- vector("list", length(files))
for(i in 1:length(files)){
	load(files[i])
	m <- do.call(rbind, samples.out)
	n <- nrow(m)
	dlist[[i]] <- data.frame(
		Phase=rep(c("CMIP3","CMIP5"),each=n/(4*length(samples.out))),
		Scenario=rep(c("SRES B1","SRES A1B","SRES A2","RCP 4.5","RCP 6.0","RCP 8.5"),each=n/(12*length(samples.out))),
		Model=rep(models[[i]],each=n/(4*length(samples.out))),
		Var=rep(c("Temperature","Precipitation"),each=n/(2*length(samples.out))),
		m,
		stringsAsFactors=F
	)
	print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(samples.out,n,m,i,files,models,dlist)
names(d)[6:ncol(d)] <- gsub("X", "", names(d)[6:ncol(d)])
gc()
#save.image("../../final/region_samples_data.RData")

# @knitr organize
f <- function(i, n, index, multiplier){
	name.tmp <- samples.names[i]
	rsd1 <- subset(d, Location==name.tmp)
	scen.levels <- unique(rsd1$Scenario)
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
	rownames(rsd1) <- NULL
	grp <- rep(names(region.names.out), times=sapply(region.names.out, length))[i]
	dir.create(outDir <- file.path("../../final/region_files_GCM/samples", grp, name.tmp), recursive=T, showWarnings=F)
	rsd1$Val <- round(rsd1$Val,1)*multiplier[1]
	rsd1$Prob <- round(rsd1$Prob,8)*multiplier[2]
	if(i==1) x <- subset(rsd1, Var=="Precipitation") else x <- NULL
	gc()
	rsd <- unlist(subset(rsd1, Var=="Precipitation")[,index])
	names(rsd) <- NULL
	save(rsd, file=paste0(outDir, "/", "precipitation.RData"))
	gc()
	rsd <- unlist(subset(rsd1, Var=="Temperature")[,index])
	names(rsd) <- NULL
	save(rsd, file=paste0(outDir, "/", "temperature.RData"))
	gc()
	print(i)
	x
}

# @knitr save
library(reshape2)
library(parallel)
samples.columns <- 6:7
samples.multipliers <- c(1e1, 1e8)

out <- mclapply(1:length(samples.names), f, n=n.samples, index=samples.columns, multiplier=samples.multipliers, mc.cores=27)
gcm.samples.df <- out[[1]]
gcm.samples.df[,samples.columns] <- NA
gcm.samples.df$Var <- NA
gcm.samples.df$Location <- NA

rm(d, f, out)
load("../../final/meta.RData")
save.image("../../final/meta.RData")
