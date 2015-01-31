# @knitr setup
comArgs <- commandArgs(TRUE)
if(any(comArgs=="akcan2km")) domain <- "akcan2km" else if(any(comArgs=="world10min")) domain <- "world10min" else stop("Unquoted 'akcan2km' or 'world10min' argument not supplied.")

setwd("/workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/data/cities")

files <- list.files(pattern=paste0(domain, ".RData$"))
files <- files[substr(files, 1, 3) != "CRU"]

models <- list(
	c("CCCMAcgcm31","CCSM4"),
	c("GFDLcm21","GFDLcm3"),
	c("MIROC32m","GISSe2-r"),
	c("MPIecham5","IPSLcm5a-lr"),
	c("ukmoHADcm3","MRIcgcm3")
)

# @knitr load
for(i in 1:length(files)){
	load(files[i])
	cities.meta <- d.cities
	m <- as.numeric(d)
	n <- length(m)
	d.tmp <- data.frame(
		Phase=rep(c("CMIP3","CMIP5"), each=n/4),
		Scenario=rep(c("SRES B1","SRES A1B","SRES A2","RCP 4.5","RCP 6.0","RCP 8.5"),each=n/12),
		Model=rep(models[[i]], each=n/4),
		Var=rep(c("Temperature","Precipitation"), each=n/2),
		Location=rep(paste0(cities.meta$loc, ", ", cities.meta$region), each=n/(12*nrow(cities.meta))),
		Val=m,
		stringsAsFactors=F
	)
	if(i==1) d.hold <- d.tmp else d.hold <- rbind(d.hold, d.tmp)
	gc()
	print(i)
}
d <- d.hold
rm(d.hold, n, m, i, files, models, d.tmp, domain)

# @knitr organize
d$Val[d$Var=="Temperature"] <- round(d$Val[d$Var=="Temperature"], 1)
d$Val[d$Var=="Precipitation"] <- round(d$Val[d$Var=="Precipitation"])
d$Month <- month.abb
d$Year <- results.years #rep(1870:2099,each=12)
d$Month <- factor(d$Month, levels=d$Month[1:12])
d$Scenario <- factor(d$Scenario, levels=unique(d$Scenario))
cities.meta$loc <- paste0(cities.meta$loc, ", ", cities.meta$region)
cities.meta <- cities.meta[c(1,3:6)]
names(cities.meta) <- c("Country", "Location", "Population", "Lat", "Lon")
d$Decade <- paste0(substr(d$Year,1,3),0)
gc()
d.cities <- d
rm(d)
gc()
#save.image("../final/data_cities.RData")

# @knitr save
library(parallel)

f <- function(i){
	name.tmp <- gsub("`", "", gsub("~", "", gsub("?", "", gsub("\\'", "", cities.meta$Location[i]))))
	city.dat <- subset(d.cities, Location==cities.meta$Location[i])
	save(city.dat, file=paste0("../final/city_files_GCM/", gsub(", ", "--", name.tmp), "_", domain, ".RData"))
	print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores=32)
