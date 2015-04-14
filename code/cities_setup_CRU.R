# @knitr setup
comArgs <- commandArgs(TRUE)
if(length(comArgs)) for(i in 1:length(comArgs)) eval(parse(text=comArgs[[i]]))
if(!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if(!exists("cities.batch")) cities.batch <- 1
if(!exists("cru")) cru <- "32"

setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/cities")

library(data.table)

files <- list.files(pattern=paste0("^CRU", cru, ".*.batch", cities.batch, "_", domain, ".RData$"))

models <- paste0("CRU", cru)

# @knitr load
dlist <- vector("list", length(files))
for(i in 1:length(files)){
	load(files[i])
	cities.meta <- d.cities
	m <- as.numeric(d)
	n <- length(m)
	dlist[[i]] <- data.frame(
		Var=rep(c("Temperature","Precipitation"), each=n/(2*nrow(cities.meta))),
		Location=rep(paste0(cities.meta$loc, ", ", cities.meta$region), each=n/nrow(cities.meta)),
		Val=m,
		stringsAsFactors=F
	)
	gc()
	print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(n, m, i, files, models, dlist, cities.batch)

# @knitr organize
d$Val[d$Var=="Temperature"] <- round(d$Val[d$Var=="Temperature"],1)
d$Val[d$Var=="Precipitation"] <- round(d$Val[d$Var=="Precipitation"])
d$Month <- month.abb
d$Year <- results.years #rep(1901:2009,each=12)
d$Month <- factor(d$Month, levels=d$Month[1:12])
cities.meta$loc <- paste0(cities.meta$loc, ", ", cities.meta$region)
cities.meta <- cities.meta[c(1,3:6)]
names(cities.meta) <- c("Country", "Location", "Population", "Lat", "Lon")
d$Decade <- as.character(d$Year - d$Year %% 10)
gc()
d.cities.cru <- d
rm(d, results.years)
gc()
#save.image("../final/data_cities_CRU31.RData")

# @knitr save
library(parallel)
f <- function(i, overwrite=FALSE){
	name.tmp <- gsub("\\.", "PER", gsub("/", "FSLASH", gsub("`", "", gsub("~", "", gsub("?", "", gsub("\\'", "APOS", cities.meta$Location[i]))))))
	filename <- paste0("../final/city_files_CRU", cru, "/", domain, "/", gsub(", ", "--", name.tmp), "__", domain, ".RData")
	if(overwrite | !file.exists(filename)){
		city.cru.dat <- subset(d.cities.cru, Location==cities.meta$Location[i])
		save(city.cru.dat, file=filename)
	}
	print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores=32)

rm(cru, d.cities.cru, f, domain)
save(cities.meta, file="../final/cities_meta.RData") # only necessary one time out of all versions of CRU and GCMs, 10-minute resolution inputs provide larger city set (NWT)
#load("../final/meta.RData")
#save.image("../final/meta.RData")
