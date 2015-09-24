# @knitr setup
setwd("/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/regional/stats")

library(data.table)

cru <- "32"
files <- list.files(pattern=paste0("^CRU", cru, ".*.regions_stats.RData$"))

models <- paste0("CRU", cru)

# @knitr load
dlist <- vector("list", length(files))
for(i in 1:length(files)){
	load(files[i])
	m <- do.call(rbind,stats.out)
	n <- nrow(m)
	p <- length(stats.out)
	dlist[[i]] <- data.frame(
		Var=rep(c("Temperature","Precipitation"),each=n/(2*p)),
		Location=rep(names(stats.out),each=n/p),
		m,
		stringsAsFactors=F
	)
	print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(stats.out,n,m,i,p,files,models,dlist)

# @knitr organize
d[d$Var=="Temperature", 3:ncol(d)] <- round(d[d$Var=="Temperature",  3:ncol(d)],1)
d[d$Var=="Precipitation", 3:ncol(d)] <- round(d[d$Var=="Precipitation",  3:ncol(d)])
d$Year <- results.years
rm(results.years)
d$Month <- month.abb
d$Month <- factor(d$Month, levels=month.abb)
d$Decade <- paste0(substr(d$Year,1,3),0)

#agg.stat.names <- c("Mean", "Std. Dev.", "5th percentile","10th percentile", "25th percentile", "50th (Median)", "75th percentile", "90th percentile", "95th percentile")
stats.columns.cru <- seq(which(names(d)=="Mean"), length.out=length(agg.stat.names))

# @knitr save
library(parallel)

f <- function(i, cru){
	name.tmp <- as.character(unlist(region.names.out))[i]
	region.cru.dat <- subset(d, Location==name.tmp)
	names(region.cru.dat)[stats.columns.cru] <- agg.stat.colnames
	grp <- rep(names(region.names.out), times=sapply(region.names.out, length))[i]
	dir.create(outDir <- paste0("../../final/region_files_CRU", cru, "/stats/", grp, "/", name.tmp), recursive=T, showWarnings=F)
	save(region.cru.dat, file=file.path(outDir, "stats_climate.RData"))
	print(i)
}

mclapply(1:length(unique(d$Location)), f, cru=cru, mc.cores=32)
