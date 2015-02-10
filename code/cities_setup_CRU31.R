# @knitr setup
comArgs <- commandArgs(TRUE)
if(length(comArgs)) for(i in 1:length(comArgs)) eval(parse(text=comArgs[[i]]))
if(!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if(!exists("cities.batch")) cities.batch <- ""

setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/cities")

library(data.table)

files <- list.files(pattern=paste0("^CRU31.*.batch", cities.batch, ".*.", domain, ".RData$"))

models <- "CRU31"

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
d$Decade <- paste0(substr(d$Year,1,3),0)
gc()
d.cities.cru31 <- d
rm(d, results.years)
gc()
#save.image("../final/data_cities_CRU31.RData")

# @knitr save
library(parallel)
f <- function(i){
	name.tmp <- gsub("`", "", gsub("~", "", gsub("?", "", gsub("\\'", "APOS", cities.meta$Location[i]))))
	city.cru.dat <- subset(d.cities.cru31, Location==cities.meta$Location[i])
	save(city.cru.dat, file=paste0("../final/city_files_CRU/", gsub(", ", "--", name.tmp), "__", domain, ".RData"))
	print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores=32)

assign(paste0("cities.meta.", domain), cities.meta)

rm(d.cities.cru31, f, cities.meta, domain)
save(cities.meta.akcan2km, file="../final/cities_meta_akcan2km.RData")
#load("../final/meta.RData")
#save.image("../final/meta.RData")
