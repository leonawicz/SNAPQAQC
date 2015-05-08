##############################################################
#### Spatial autocorrelation and clump size distributions ####
##############################################################

#### Script author:  Matthew Leonawicz ####
#### Maintainted by: Matthew Leonawicz ####
#### Last updated:   05/08/2015        ####

# @knitr setup1
# server
library(raster)
library(parallel)
library(data.table)

setwd("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/workspaces")
r <- readAll(raster("../../DataExtraction/data/tas_mean_C_AR5_GFDL-CM3_rcp60_01_2062.tif"))

f <- function(i, r, v){
	rc <- clump(Which(r==v[i]))
	frq <- freq(rc)
	tb <- table(frq[!is.na(frq[, "value"]), "count"])
	data.table(Area=4*as.integer(names(tb)), Freq=as.integer(tb), Val=v[i]) # 2 x 2 km grid cells
}

system.time(d <- rbindlist(mclapply(1:length(v), f, r=r, v=v, mc.cores=32)))
setcolorder(d, c("Val", "Area", "Freq"))

d <- rbind(data.table(Val=unique(d[, Val]), Area=0, Freq=0), d)
d[, RelSum := cumsum(Freq)/sum(Freq), by=Val]
d[, Binned := cut(Val, seq(min(Val), max(Val), length.out=10), include.lowest=T)]
d[, TotalArea := Area*Freq]
d[, RelTotalArea := cumsum(TotalArea)/sum(TotalArea), by=Val]
save(d, file="dt_clump_t2km.RData")

# @knitr setup2
# local
library(raster)
library(data.table)
library(ggplot2)

setwd("C:/github/SNAPQAQC/workspaces")
load("dt_clump_t2km.RData")

r <- raster("C:/github/DataExtraction/data/tas_mean_C_AR5_GFDL-CM3_rcp60_01_2062.tif")

# @knitr autocorr
m_i <- Moran(r)
g_c <- Geary(r)
r.ml <- MoranLocal(r)
plot(r.ml)

# @knitr values
v <- sort(unique(r[]))
v <- v[!is.na(v)]
n <- length(v)
n.dat <- length(which(!is.na(r[])))

# @knitr clump_size_plots
# Scaled cumulative clump size frequency by clump size per temperature value
xlb <- expression("Area ("~km^2~")")
ylb <- "Relative cumulative sum"
g1 <- ggplot(data=d, aes(x=Area, y=RelSum, group=Val)) + geom_step() + ggtitle("Relative cumulative clump size frequency by clump size\nper temperature value") + xlab(xlb) + ylab(ylb)

# Scaled cumulative clump size frequency by clump size given temperature 
g2 <- g1 + facet_wrap(~ Binned, scales="free")

# Scaled cumulative coverage area by clump size per temperature value
g3 <- ggplot(data=d, aes(x=Area, y=RelTotalArea, group=Val)) + geom_step() + ggtitle("Relative cumulative total area\nper temperature value (C)") + xlab(xlb) + ylab(ylb)

# Scaled cumulative coverage area by clump size given temperature
g4 <- g3 + facet_wrap(~ Binned, scales="free")

# @knitr clump_size_plot1
g1
# @knitr clump_size_plot2
g2
# @knitr clump_size_plot3
g3
# @knitr clump_size_plot4
g4

# @knitr summarize
d2 <- d[, sum(TotalArea), by=Area]
setkey(d2, Area)
d2[, CS := cumsum(V1)]

# Cumulative coverage area by clump size < 100 km^2, integrated across all temperatures
g5 <- ggplot(data=d2[Area <= 100,], aes(x=Area, y=CS)) + geom_step() + ggtitle("Relative cumulative total area by temperature (C)") + xlab(xlab) + ylab(ylab)
n.cells <- sum(d2[, V1])/4
p <- sum(d2[Area > 4, V1])/4

# @knitr clump_size_plot5
g5
