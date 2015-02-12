# @knitr setup
comArgs <- commandArgs(TRUE)
if(length(comArgs)) for(i in 1:length(comArgs)) eval(parse(text=comArgs[[i]]))
if(!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if(!exists("cities.batch")) cities.batch <- 1

setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/cities")

library(data.table)

files <- list.files(pattern=paste0("batch", cities.batch, "_", domain, ".RData$"))
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
	cities.meta <- d.cities
	m <- as.numeric(d)
	n <- length(m)
	dlist[[i]] <- data.frame(
		Phase=rep(c("CMIP3","CMIP5"), each=n/4),
		Scenario=rep(c("SRES B1","SRES A1B","SRES A2","RCP 4.5","RCP 6.0","RCP 8.5"),each=n/12),
		Model=rep(models[[i]], each=n/4),
		Var=rep(c("Temperature","Precipitation"), each=n/2),
		Location=rep(paste0(cities.meta$loc, ", ", cities.meta$region), each=n/(12*nrow(cities.meta))),
		Val=m,
		stringsAsFactors=F
	)
	gc()
	print(i)
}
d <- as.data.frame(rbindlist(dlist))
rm(n, m, i, files, models, dlist)
gc()

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
d$Decade <- as.character(10*(d$Year %/% 10))
gc()
d.cities <- d
rm(d, results.years)
gc()
#save.image("../final/data_cities.RData")

# @knitr save
# Save individual R workspace files for each city for use in master QAQC Shiny app
library(parallel)

f <- function(i){
	name.tmp <- gsub("\\.", "PER", gsub("/", "FSLASH", gsub("`", "", gsub("~", "", gsub("?", "", gsub("\\'", "APOS", cities.meta$Location[i]))))))
	city.dat <- subset(d.cities, Location==cities.meta$Location[i])
	save(city.dat, file=paste0("../final/city_files_GCM/", gsub(", ", "--", name.tmp), "__", domain, ".RData"))
	print(i)
}

mclapply(1:length(cities.meta$Location), f, mc.cores=32)
