# @knitr setup
setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/regional/samples")

library(data.table)
library(reshape2)

cru <- "32"
files <- list.files(pattern=paste0("^CRU", cru, ".*.regions_samples.RData$"))

models <- paste0("CRU", cru)

# @knitr load
dlist <- vector("list", length(files))
for(i in 1:length(files)){
	load(files[i])
	m <- do.call(rbind, samples.out)
	n <- nrow(m)
	dlist[[i]] <- data.frame(
		Var=rep(c("Temperature","Precipitation"),each=n/(2*length(samples.out))),
		m,
		stringsAsFactors=F
	)
	print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(samples.out,n,m,i,files,models,dlist)
names(d)[3:ncol(d)] <- gsub("X", "", names(d)[3:ncol(d)])
#save.image("../../final/CRU31_region_samples_data.RData")

# @knitr organize
f <- function(i, n, index, multiplier, cru){
	name.tmp <- samples.names[i]
	rsd <- subset(d, Location==name.tmp)
	rsd <- melt(rsd, id.vars=names(rsd)[1:2], measure.vars=names(rsd)[-c(1:2)], variable.name="Time", value.name="Vals_Probs")
	rsd$vals.ind <- rep(rep(c(T,F), each=n), length=nrow(rsd))
	rsd$dcastIsDumb <- rep(1:n, length=nrow(rsd))
	rsd <- dcast(rsd, Var + Location + Time + dcastIsDumb ~ vals.ind, value.var="Vals_Probs")
	names(rsd)[6:5] <- c("Val", "Prob")
	rsd$Year <- substr(as.character(rsd$Time), 1, 4)
	rsd$Month <- substr(as.character(rsd$Time), 6, 8)
	rsd$Month <- factor(rsd$Month, levels=month.abb)
	rsd$Decade <- paste0(substr(rsd$Year,1,3),0)
	rownames(rsd) <- NULL
	rsd <- rsd[c(1:2,6,5,7:9)]
	grp <- rep(names(region.names.out), times=sapply(region.names.out, length))[i]
	dir.create(outDir <- paste0("../../final/region_files_CRU", cru, "/samples/", grp, "/", name.tmp), recursive=T, showWarnings=F)
	rsd$Val <- round(rsd$Val,1)*multiplier[1]
	rsd$Prob <- round(rsd$Prob,8)*multiplier[2]
	if(i==1) x <- subset(rsd, Var=="Precipitation") else x <- NULL
	gc()
	rsd.cru <- unlist(subset(rsd, Var=="Precipitation")[,index])
	names(rsd.cru) <- NULL
	save(rsd.cru, file=paste0(outDir, "/", "precipitation.RData"))
	gc()
	rsd.cru <- unlist(subset(rsd, Var=="Temperature")[,index])
	names(rsd.cru) <- NULL
	save(rsd.cru, file=paste0(outDir, "/", "temperature.RData"))
	gc()
	print(i)
	x
}

# @knitr save
library(parallel)

samples.columns.cru <- 3:4
samples.multipliers.cru <- c(1e1, 1e8)

out <- mclapply(1:length(samples.names), f, n=n.samples, index=samples.columns.cru, multiplier=samples.multipliers.cru, cru=cru, mc.cores=32)
cru.samples.df <- out[[1]]
cru.samples.df[,samples.columns.cru] <- NA
cru.samples.df$Var <- NA
cru.samples.df$Location <- NA
if(as.numeric(cru)==31) cru31.samples.df <- cru.samples.df else if(as.numeric(cru)==32) cru32.samples.df <- cru.samples.df
rm(cru, cru.samples.df, d, f, out)
load("../../final/meta.RData")
save.image("../../final/meta.RData")
