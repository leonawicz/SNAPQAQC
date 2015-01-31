# @knitr setup
comArgs <- commandArgs(TRUE)
if(any(comArgs=="akcan2km")) domain <- "akcan2km" else if(any(comArgs=="world10min")) domain <- "world10min" else stop("Unquoted 'akcan2km' or 'world10min' argument not supplied.")

setwd("/workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/data/cities")

files <- list.files(pattern=paste0("^CRU31.*.", domain,".RData$")

models <- "CRU31"

# @knitr load
for(i in 1:length(files)){
	load(files[i])
	cities.meta <- d.cities
	m <- as.numeric(d)
	n <- length(m)
	d.tmp <- data.frame(
		Var=rep(c("Temperature","Precipitation"), each=n/(2*nrow(cities.meta))),
		Location=rep(paste0(cities.meta$loc, ", ", cities.meta$region), each=n/nrow(cities.meta)),
		Val=m,
		stringsAsFactors=F
	)
	if(i==1) d.hold <- d.tmp else d.hold <- rbind(d.hold, d.tmp)
	gc()
	print(i)
}
d <- d.hold

rm(d.hold, n, m, i, files, models, d.tmp)

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
rm(d)
gc()
#save.image("../final/data_cities_CRU31.RData")

# @knitr save
library(parallel)
f <- function(i){
	name.tmp <- gsub("`", "", gsub("~", "", gsub("?", "", gsub("\\'", "", cities.meta$Location[i]))))
	city.cru.dat <- subset(d.cities.cru31, Location==cities.meta$Location[i])
	save(city.cru.dat, file=paste0("../final/city_files_CRU/", gsub(", ", "--", name.tmp), "_", domain, ".RData"))
	print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores=32)

assign(pste0("cities.meta.", domain), cities.meta)

rm(d.cities, f, cities.meta)
load("../final/meta.RData")
save.image("../final/meta.RData")
