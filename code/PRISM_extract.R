# @knitr setup
setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final")
library(raster)
domain <- "akcan2km"
#domain <- "world10min"

mainDir <- "/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs"
pDir <- if(domain=="akcan2km") "pr/ak83albers" else "pr/wgs84"
tDir <- if(domain=="akcan2km") "tas/ak83albers" else "tas/wgs84"

locs <- read.csv("/workspace/Shared/Users/mfleonawicz/github/statistics/AR5_scripts/AR5_QAQC/locs.csv")
s.p <- stack(list.files(file.path(mainDir, pDir), pattern=".tif$", full=TRUE))
s.t <- stack(list.files(file.path(mainDir, tDir), pattern=".tif$", full=TRUE))

# @knitr cities
# Select cities and extract data
locs <- locs[locs$pop > 10,]
l <- paste(locs$region, locs$loc)
lu <- unique(l)
dup <- which(duplicated(l))
l.dup.u <- unique(l[dup])
drp <- c()
for(i in 1:length(l.dup.u)){
	ind <- which(l==l.dup.u[i])
	ind <- ind[-which.max(locs$pop[ind])[1]]
	drp <- c(drp, ind)
}
cities <- locs[-drp,]
prism.cities <- paste(cities$loc, cities$region, sep=", ")
cities <- if(domain=="akcan2km") cbind(cities$lon_albers, cities$lat_albers) else if(domain=="world10min") cbind(cities$lon, cities$lat)

prism.p <- extract(s.p, cities)
prism.t <- round(extract(s.t, cities), 1)

#rm(list=ls()[!(ls() %in% c("prism.p", "prism.t", "prism.cities"))])

# @knitr save
# Save results
load("cc4lite/cc4lite.RData")
save(d, locs, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_prism.RData")
