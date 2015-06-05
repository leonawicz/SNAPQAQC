# @knitr setup
setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/final")
library(raster)
# At this time, no buffered 10-minute resolution WGS84 mean value version extracted for domain='world10min'
domain <- "akcan2km"
#domain <- "world10min"

mainDir <- "/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs"
pDir <- if(domain=="akcan2km") "pr/ak83albers" else "pr/wgs84"
tDir <- if(domain=="akcan2km") "tas/ak83albers" else "tas/wgs84"

locs <- read.csv("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/locs.csv")
s.p <- stack(list.files(file.path(mainDir, pDir), pattern=".tif$", full=TRUE))
s.t <- stack(list.files(file.path(mainDir, tDir), pattern=".tif$", full=TRUE))

# @knitr cities
# Select cities and extract data
locs <- subset(locs, region!="Northwest Territories")
#locs <- locs[locs$pop > 10,]
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
file1 <- "cc4lite/cc4lite_cru31_prism_akcan2km.RData"
file2 <- "cc4lite/cc4lite_cru32_prism_akcan2km.RData"
file3 <- "cc4lite/cc4lite_cru3132_prism_akcan2km.RData"
file4 <- "cc4lite/cc4lite_cru31_prism_world10min.RData"
file5 <- "cc4lite/cc4lite_cru32_prism_world10min.RData"
file6 <- "cc4lite/cc4lite_cru3132_prism_world10min.RData"
if(domain=="akcan2km"){
	load(paste0("cc4lite/cc4lite_cru31_akcan2km.RData"))
	load(paste0("cc4lite/cc4lite_cru32_akcan2km.RData"))
	save(d.2km, locs, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_prism_akcan2km.RData")
	save(d.2km, d.cru31.2km, locs, prism.p, prism.t, prism.cities, file=file1)
	save(d.2km, d.cru32.2km, locs, prism.p, prism.t, prism.cities, file=file2)
	save(d.2km, d.cru31.2km, d.cru32.2km, locs, prism.p, prism.t, prism.cities, file=file3)
} else if(domain=="world10min") {
	load(paste0("cc4lite/cc4lite_cru31_world10min.RData"))
	load(paste0("cc4lite/cc4lite_cru32_world10min.RData"))
	save(d.10min, locs, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_prism_world10min.RData")
	save(d.10min, d.cru31.10min, locs, prism.p, prism.t, prism.cities, file=file4)
	save(d.10min, d.cru32.10min, locs, prism.p, prism.t, prism.cities, file=file5)
	save(d.10min, d.cru31.10min, d.cru32.10min, locs, prism.p, prism.t, prism.cities, file=file6)
}

if(all(file.exists(file3, file6))){
	load(file3)
	locs.2km <- locs
	load(file6)
	locs.10min <- locs
	save(d.2km, d.10min, locs.2km, locs.10min, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_prism_2km10min.RData")
	save(d.2km, d.10min, d.cru31.2km, d.cru31.10min, locs.2km, locs.10min, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_cru31_prism_2km10min.RData")
	save(d.2km, d.10min, d.cru32.2km, d.cru32.10min, locs.2km, locs.10min, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_cru32_prism_2km10min.RData")
	save(d.2km, d.10min, d.cru31.2km, d.cru31.10min, d.cru32.2km, d.cru32.10min, locs.2km, locs.10min, prism.p, prism.t, prism.cities, file="cc4lite/cc4lite_cru3132_prism_2km10min.RData")
}
