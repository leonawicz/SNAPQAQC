# @knitr setup
library(data.table)
library(raster)
library(parallel)

# Load nested list of cell indices defining groups of shapefile polygons
load("/workspace/UA/mfleonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN1km.RData")
load("/workspace/UA/mfleonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN1km_rmNA.RData")

model.index <- 1
dirs <- list.files("/big_scratch/apbennett/Calibration/FinalCalib", pattern=".*.sres.*.", full=T)
repid <- 1:200
years <- 1950
mainDirs <- rep(paste0(dirs,"/Maps")[model.index], each=length(repid))
modnames <- basename(dirname(mainDirs)) # Should only be one name at a time due to $ScenMod column naming below
mainDir <- mainDirs[model.index]
outDir <- "/workspace/UA/mfleonawicz/projects/SNAPQAQC/plots/ageDensityCalibration"

getAgeVegStats <- function(i, mainDir, years=NULL, cells.list, loc, n.samples, v.id=1:6, ...){
	rep.lab <- paste0("_",c(0:199),"_")[i]
	patA <- gsub("expression","",paste(bquote(expression("^A.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesA <- list.files(mainDir, pattern=patA, full=T)
	patV <- gsub("expression","",paste(bquote(expression("^V.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesV <- list.files(mainDir, pattern=patV, full=T)
	files.years <- as.numeric(substr(filesA,nchar(filesA)-7,nchar(filesA)-4))
	if(is.numeric(years)) { p <- files.years %in% years; if(any(p)) { filesA <- filesA[which(p)]; filesV <- filesV[which(p)]; files.years <- files.years[which(p)] } }
	years <- unique(files.years)
	r <- getValues(raster(filesA[1])) # use as a template
	data.ind <- which(!is.na(r))
	
	# Density estimation
	denFun <- function(x, n, adj=0.25, min.zero=TRUE, diversify=FALSE){
		x <- x[!is.na(x)]
		lx <- length(x)
		if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
		if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
		dif <- diff(range(x))
		z <- density(x, adjust=adj, n=n, from=min(x)-max(1, 0.05*dif), to=max(x)+max(1, 0.05*dif))
		if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+max(1, 0.05*dif))
		as.numeric(c(z$x, z$y))
	}
	
	xl <- xl.all <- vector("list", length(years))
	for(j in 1:length(years)){
		a <- getValues(raster(filesA[j]))[data.ind]
		v <- getValues(raster(filesV[j]))[data.ind]
		a.tmp <- a[cells.list]
		v.tmp <- v[cells.list]
		v.tmp.n <- tapply(cells.list, v.tmp, length)
		den.tmp <- tapply(a.tmp, v.tmp, denFun, n=n.samples)
		xl[[j]] <- data.frame(Location=loc, VegID=rep(as.numeric(names(den.tmp)), each=2*n.samples), VegArea=rep(v.tmp.n, each=2*n.samples), Val=unlist(den.tmp), Year=years[j], stringsAsFactors=FALSE)
		xl.all[[j]] <- tapply(a.tmp, v.tmp, c)
		print(paste("Age and veg maps:", years[j]))
	}
	
	x <- rbindlist(xl)
	x$Replicate <- i
	x <- x[x$VegID %in% v.id,]
	list(x, xl.all)
}

set.seed(47)
system.time( x <- mclapply(repid, getAgeVegStats, mainDir=mainDir, years=years, cells.list=cells_shp_list_rmNA[[5]][[1]], loc=region.names.out[[5]][[1]], n.samples=100, v.id=1, mc.cores=32) )
ind <- sample(1:200, 30)

# Density
x.den <- lapply(x, function(x) x[[1]]$Val)[ind]
n <- length(x.den[[1]])/2

# All values
x.all <- unlist(lapply(x, function(x) x[[2]][[1]][[1]])[ind])

# Second density estimation
denFun <- function(x, n=20, adj=0.25, min.zero=TRUE, diversify=FALSE, missing.veg.NA=TRUE, fire=FALSE){
	if(all(is.na(x))) return(rep(NA, 2*n))
	x <- x[!is.na(x)]
	lx <- length(x)
	if(sum(x)==0 & missing.veg.NA & !fire) return(rep(NA, 2*n))
	if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
	if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
	dif <- diff(range(x))
	z <- density(x, adjust=adj, n=n, from=min(x)-max(1, 0.05*dif), to=max(x)+max(1, 0.05*dif))
	if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+max(1, 0.05*dif))
	as.numeric(c(z$x, z$y))
}

# Bootstrapping
btfun <- function(p, n.samples=length(p)/2, n.boot=10000, interp=FALSE, n.interp=1000, ...){
	if(!length(p)) return(p)
	if(all(is.na(p))) return(rep(NA, n.boot))
	p <- list(x=p[1:n.samples], y=p[(n.samples+1):(2*n.samples)])
	if(interp && length(unique(p[1:n.samples])) > 1) p <- approx(p$x, p$y, n=n.interp)
	p <- round(sample(p$x, n.boot, prob=p$y, rep=T), ...)
	p
}

x.boot1 <- unlist(lapply(x.den, btfun, n.boot=100, interp=TRUE, n.interp=1000))

x.den2 <- denFun(x.boot1, n=n)

x.boot2 <- btfun(x.den2, n.boot=100000, interp=TRUE, n.interp=1000)

png(paste0(outDir, "/test2_", length(ind), "reps.png"), res=200, width=3000, height=1000)
layout(matrix(1:3, 1))
hist(x.all, breaks=20, ylim=c(0, 0.004), xlab="", main=paste("Histogram of all pixels,", length(ind), "replicates"), col="gray40", freq=F)
for(i in 1:length(x.den)) lines(x.den[[i]][1:n], x.den[[i]][(n+1):(2*n)], col="#FF000015", lwd=1)
legend("topleft", "Probability density estimation", lwd=2, col="red", bty="n")
hist(x.boot1, breaks=20, ylim=c(0, 0.004), xlab="Vegetation age (years)", ylab="", main=paste(length(ind), "pooled boostrap samples and mean density"), col="gray40", freq=F)
lines(x.den2[1:n], x.den2[(n+1):(2*n)], col="red", lwd=2)
hist(x.boot2, breaks=20, ylim=c(0, 0.004), xlab="", ylab="", main="Final bootstrap sample", col="gray40", freq=F)
dev.off()


