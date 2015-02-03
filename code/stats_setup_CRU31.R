# @knitr setup
setwd("/workspace/UA/mfleonawicz/leonawicz/Projects/active/AR4_AR5_comparisons/data/regional/stats")
files <- list.files(pattern="^CRU31.*.regions_stats.RData$")

models <- "CRU31"

# @knitr load
for(i in 1:length(files)){
	load(files[i])
	m <- do.call(rbind,stats.out)
	n <- nrow(m)
	p <- length(stats.out)
	d.tmp <- data.frame(
		Var=rep(c("Temperature","Precipitation"),each=n/(2*p)),
		Location=rep(names(stats.out),each=n/p),
		m,
		stringsAsFactors=F
	)
	if(i==1) d <- d.tmp else d <- rbind(d,d.tmp)
	print(i)
}

rm(stats.out,n,m,i,p,files,models,d.tmp)

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
#save.image("../../final/CRU31_region_stats_data.RData")

# @knitr save
library(parallel)

f <- function(i){
	name.tmp <- as.character(unlist(region.names.out))[i]
	#assign(name.tmp, subset(d.cities, Location==cities.meta$Location[i]))
	#save(list=c(name.tmp), file=paste0("../final/city_files_GCM/", gsub(", ", "--", name.tmp), ".RData"))
	region.cru.dat <- subset(d, Location==name.tmp)
	names(region.cru.dat)[stats.columns.cru] <- agg.stat.colnames
	grp <- rep(names(region.names.out), times=sapply(region.names.out, length))[i]
	dir.create(outDir <- file.path("../../final/region_files_CRU/stats", grp, name.tmp), recursive=T, showWarnings=F)
	save(region.cru.dat, file=file.path(outDir, "stats_climate.RData"))
	print(i)
}

mclapply(1:length(unique(d$Location)), f, mc.cores=32)
