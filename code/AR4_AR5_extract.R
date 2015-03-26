# @knitr clargs
comArgs <- commandArgs(TRUE)
if(length(comArgs)) for(i in 1:length(comArgs)) eval(parse(text=comArgs[[i]]))

if(!exists("batch")) stop("Model pair batch integer not supplied.")
if(!exists("domain")) stop("domain argument not provided. Must be either 'akcan2km' or 'world10min'")
if(!exists("regions")) regions <- FALSE
if(!exists("cities")) cities <- FALSE
if(!(regions | cities)) stop("regions and cities both FALSE. Nothing to process.")

# @knitr packages
library(raster)
library(maptools)
library(parallel)
library(plyr)

# @knitr source
if(domain=="akcan2km"){ # For regions and/or cities
	topDir <- file.path("/Data/Base_Data/Climate/AK_CAN_2km",c("historical","projected")) # files are not read, but metadata parsed from filenames list
	if(regions){
		load("/workspace/UA/mfleonawicz/leonawicz/projects/DataExtraction/workspaces/shapes2cells_AKCAN2km_5pct.RData")
	} else cells_shp_list_5pct <- region.names.out <- n.shp <- NULL
} else if(domain=="world10min") { # Currently for cities only
	topDir <- file.path("/Data/Base_Data/Climate/World/World_10min",c("historical","projected")) # files are not read, but metadata parsed from filenames list
	cells_shp_list_5pct <- region.names.out <- n.shp <- NULL
	#### Need to insert a load() command analogous to that above for regions
}

locs <- read.csv("/workspace/Shared/Users/mfleonawicz/github/statistics/AR5_scripts/AR5_QAQC/locs.csv")

# @knitr setup
varid <- rep(c("tas","pr"),each=6)
scenid <- if(length(topDir)==2){
			c(rep("",12), rep(c("sresb1","sresa1b","sresa2","rcp45","rcp60","rcp85"),2))
		} else if(basename(topDir)=="projected") {
			rep(c("sresb1","sresa1b","sresa2","rcp45","rcp60","rcp85"),2)
		} else if(basename(topDir)=="historical") {
			rep("",12)
		}
arid <- rep(rep(c("AR4_CMIP3_models","AR5_CMIP5_models"),each=3),2)

model.pairs <- list(
	c("cccma-cgcm3-1-t47","CCSM4"),
	c("gfdl-cm2-1","GFDL-CM3"),
	c("miroc3-2-medres","GISS-E2-R"),
	c("mpi-echam5","IPSL-CM5A-LR"),
	c("ukmo-hadcm3","MRI-CGCM3")
)

if(exists("years")) yr1 <- years[1] else yr1 <- 1870
if(exists("years")) yr2 <- tail(years, 1) else yr2 <- 2099
years <- yr1:yr2

#if(!is.null(cities)){
#	cities <- subset(locs, pop >= 1000, c("lon_albers","lat_albers"))
#	d.cities <- subset(locs, pop >= 1000)[-c(which(names(locs) %in% c("lon_albers","lat_albers")))]
#}
if(cities){
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
	if(exists("cities.batch")){
		batch.bounds <- round(seq(1, nrow(cities) + 1, length=11))[c(cities.batch, cities.batch + 1)] - c(0,1)
		cities <- cities[batch.bounds[1]:batch.bounds[2],]
	} else cities.batch <- 1
	d.cities <- cities[-c(which(names(locs) %in% c("lon_albers","lat_albers")))]
	cities <- if(domain=="akcan2km") cbind(cities$lon_albers, cities$lat_albers) else if(domain=="world10min") cbind(cities$lon+360, cities$lat) # +360 for PC lat lon rasters
}

n.samples <- 20
n2 <- 2*n.samples
agg.stat.colnames <- c("Mean", "SD", paste0("Pct_", c("05", 10, 25, 50, 75, 90, 95)))
agg.stat.names <- c("Mean", "Std Dev", paste0(c(5,10,25,50,75,90,95), "th percentile"))
agg.stat.names[agg.stat.names=="50th percentile"] <- "Median"

# @knitr functions1
# Density estimation
denFun <- function(x, n, variable){
	x <- x[!is.na(x)]
	dif <- diff(range(x))
	z <- density(x, adjust=2, n=n, from=min(x)-0.05*dif, to=max(x)+0.05*dif)
	if(variable=="pr" && any(z$x < 0)) z <- density(x, adjust=2, n=n, from=0, to=max(x)+0.05*dif)
	as.numeric(c(z$x, z$y))
}

# @knitr functions2
# Processing function
getData <- function(i, model, cells.list=NULL, shp.names=NULL, n.shp=NULL, seed=232, regions=TRUE, n.samples=512, cities=NULL, start.year=NULL, end.year=NULL){
	print(i)
	if(length(model)==1) model <- rep(model, i)
	for(p in 1:length(topDir)){
		if(p==2) print(paste("Process",i,"beginning projected data extraction."))
		if(p==1) scenid.tmp <- scenid[1:12] else if(p==2) scenid.tmp <- scenid[13:24]
		path <- gsub("/remove","",gsub("//","/",file.path(topDir[p],arid[i],scenid.tmp[i],model[i],varid[i])))
		files <- list.files(path,pattern=".tif$",full=T)
		if(!is.null(start.year)) files <- files[substr(files, nchar(files)-7, nchar(files)-4) >= start.year]
		if(!is.null(end.year)) files <- files[substr(files, nchar(files)-7, nchar(files)-4) <= end.year]
		
		if(!is.null(cities)){
			r <- readAll(raster(files[1])) # template done
			cells_cities <- extract(r, cities, cellnumbers=T)[,1]
			print("Raster cell indices for point locations obtained.")
		}
		
		mo.tmp <- substr(files,nchar(files)-10,nchar(files)-9)
		yr.tmp <- substr(files,nchar(files)-7,nchar(files)-4)
		yr.mo.tmp <- paste0(yr.tmp,mo.tmp)
		files <- files[order(yr.mo.tmp)]
		yr.tmp <- substr(files,nchar(files)-7,nchar(files)-4)
		yr.tmp <- substr(files,nchar(files)-7,nchar(files)-4)
		if(regions){
			seq.q <- c(0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95) #seq(0, 1, by=0.05) #Hard coded
			m <- matrix(NA, nrow=length(files), ncol=length(seq.q) + 2) # +2 for mean and SD below

			mden <- matrix(NA, nrow=2*n.samples, ncol=length(files))
			samples.list <- rapply(cells.list, f=function(x, m) m, classes="integer", how="replace", m=mden)
			names(samples.list) <- names(shp.names)
			for(l1 in 1:length(shp.names)) names(samples.list[[l1]]) <- shp.names[[l1]]
		}
		
		if(regions) for(l1 in 1:length(shp.names)) for(l2 in 1:length(shp.names[[l1]])) assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], sep="__"), m)
		if(regions) rm(m)
		gc()
		if(!is.null(cities)) m.cities <- matrix(NA, nrow=length(cells_cities), ncol=length(files))
		if((i %in% c(1,4,7,10)) | p==2){ # no need to read historical files more than once for each CMIP phase, scenario, and variable.
			for(b in 1:length(unique(yr.tmp))){
				pat <- gsub("expression","",paste(bquote(expression(".",.(unique(yr.tmp)[b]),".tif$")),collapse=""))
				files.sub <- list.files(path, pattern=pat, full=T)
				mat <- getValues(stack(files.sub))
				n <- ncol(mat)
				if(regions){
					for(l1 in 1:length(shp.names)){
						for(l2 in 1:length(shp.names[[l1]])){
							cells.tmp <- cells.list[[ names(shp.names)[l1] ]][[ shp.names[[l1]][l2] ]]
							m.tmp <- get(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], sep="__"))
							m.tmp[(1:n)+(12*(b-1)), 1] <- as.numeric(colMeans(mat[cells.tmp,], na.rm=T))
							m.tmp[(1:n)+(12*(b-1)), 2] <- as.numeric(apply(mat[cells.tmp,], 2, sd, na.rm=T))
							m.tmp[(1:n)+(12*(b-1)), 2+(1:length(seq.q))] <- t(apply(mat[cells.tmp,], 2, quantile, probs=seq.q, na.rm=T))
							assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], sep="__"), m.tmp)
							samples.tmp <- apply(mat[cells.tmp,], 2, denFun, n=n.samples, variable=varid[i]) # Time across columns
							samples.list[[ names(shp.names)[l1] ]][[ shp.names[[l1]][l2] ]][, (1:n)+(12*(b-1))] <- samples.tmp
						}
					}
					gc()
				}
				if(!is.null(cities)) m.cities[,(1:n)+(12*(b-1))] <- mat[cells_cities,]
				print(paste0("Process",i,": ",unique(yr.tmp)[b]))
			}
		} else print(paste("Process",i,"skipping historical files to avoid redundancy."))
		if(p==1){
			if(regions){
				for(l1 in 1:length(shp.names)) for(l2 in 1:length(shp.names[[l1]])) assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], "hold", sep="__"), get(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], sep="__")))
				samples.list.hold <- samples.list
			}
			if(!is.null(cities)) m.cities.hold <- m.cities
			yr.hold <- as.numeric(yr.tmp)
		} else {
			if(regions){
				for(l1 in 1:length(shp.names)){
					for(l2 in 1:length(shp.names[[l1]])){
						assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], "hold", sep="__"),
							rbind(get(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], "hold", sep="__")), get(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], sep="__"))))
					}
				}
				for(l1 in 1:length(shp.names)){
					for(l2 in 1:length(shp.names[[l1]])){
						samples.list.hold[[ names(shp.names)[l1] ]][[ shp.names[[l1]][l2] ]] <-
							cbind(samples.list.hold[[ names(shp.names)[l1] ]][[ shp.names[[l1]][l2] ]], samples.list[[ names(shp.names)[l1] ]][[ shp.names[[l1]][l2] ]]) # Time across columns
					}
				}
			}
			if(!is.null(cities)) m.cities.hold <- cbind(m.cities.hold, m.cities)
			yr.hold <- c(yr.hold,as.numeric(yr.tmp))
		}
	}
	if(regions){
		output.names <- ls(pattern="^m__.*.__hold$")
		m.list <- mget(output.names)
		regions.list <- c(m.list, list(yr=yr.hold), list(samples.list.hold))
		names(regions.list) <- c(sapply(strsplit(output.names, "__"), "[[", 3), "yr", "samples")
	}
	if(regions & is.null(cities))	l <- list(regions.list, list("No cities"))
	if(regions & !is.null(cities))	l <- list(regions.list, list(m.cities=m.cities.hold))
	if(!regions & !is.null(cities))	l <- list(c(list(yr=yr.hold), list("empty")), list(m.cities=m.cities.hold))
	l
}

# @knitr run
models <- rep(rep(model.pairs[[batch]],each=3),2)
results <- mclapply(X=1:12, FUN=getData, model=models, cells.list=cells_shp_list_5pct, shp.names=region.names.out, n.shp=n.shp,
	regions=regions, n.samples=n.samples, cities=cities, start.year=yr1, end.year=yr2, mc.cores=12)

# @knitr save
# Organize and save results
for(k in 1:2){
	if(!regions & k==1) next
	stats <- lapply(results, "[[", k)
	if(is.character(stats[[1]])) next
	len <- ifelse(k==1, 2, 0)
	if(k==1){
		#regional stats
		for(i in 1:length(stats)){
			results.years <- results[[i]][[1]][[length(results[[i]][[1]])-1]]
			for(j in 1:(length(stats[[i]])-len)) stats[[i]][[j]] <- stats[[i]][[j]][ which(results.years %in% years), ] # matrix (regional stats)
		}
		stats.out.names <- names(results[[1]][[1]])
		stats.out.names <- stats.out.names[-which(stats.out.names %in% c("yr","samples"))]
		stats.out <- list()
		for(i in 1:length(stats.out.names)) stats.out[[i]] <- do.call(rbind, lapply(stats, "[[", i))
		names(stats.out) <- stats.out.names
		time.seq <- rep(seq(years[1],length.out=nrow(stats.out[[1]])/(12*12)), each=12)
		for(i in 1:length(stats.out)){
			colnames(stats.out[[i]]) <- agg.stat.colnames
			rownames(stats.out[[i]]) <- rep(time.seq, 12)
			for(block in c(2,5,8,11)){
				block.rows <- 1:length(time.seq) + length(time.seq)*(block-1)
				na.rows <- which(is.na(stats.out[[i]][block.rows, 1]))
				stats.out[[i]][block.rows[na.rows],] <- stats.out[[i]][block.rows[na.rows]+length(time.seq),] <- stats.out[[i]][block.rows[na.rows]-length(time.seq),]
			}
		print(length(stats.out)-i)
		}
		save(stats.out, results.years, region.names.out, agg.stat.names, agg.stat.colnames,
			file=paste0("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/regional/stats/",models[1],"_",models[4],"_annual_regions_stats.RData"))
		
		# regional samples
		samples <- lapply(results, function(x) x[[1]][["samples"]])
		list2df <- function(x){
			x <- ldply(x, data.frame)
			x[,".id"] <- paste(x[,".id"], names(x)[2], sep="_")
			names(x) <- c("Location", paste("T", 1:(ncol(x)-1), sep="_"))
			x
		}
		for(l1 in 1:length(samples)) for(l2 in 1:length(samples[[l1]])) samples[[l1]][[l2]] <- list2df(samples[[l1]][[l2]])
		for(l1 in 1:length(samples)) samples[[l1]] <- do.call(rbind, samples[[l1]])
		samples.names <- sapply(strsplit(unique(samples[[1]]$Location), "_"), "[[", 1)
		samples <- do.call(rbind, samples)
		samples$Location <- rep(rep(samples.names, each=n2), 12)
		names(samples)[-1] <- paste(results.years, month.abb, sep="_")
		rownames(samples) <- NULL
		samples.out <- list()
		for(i in 1:length(samples.names)){
			samples.out[[i]] <- subset(samples, Location==samples.names[i])
			rownames(samples.out[[i]]) <- NULL
		}
		names(samples.out) <- samples.names
		nc <- ncol(samples.out[[1]])
		for(i in 1:length(samples.out)){
			for(cols in 2:nc){
				if(all(is.na(samples.out[[i]][(n2+1):(2*n2), cols]))) samples.out[[i]][(n2+1):(2*n2), cols] <- samples.out[[i]][(2*n2+1):(3*n2), cols] <- samples.out[[i]][1:n2, cols]
				if(all(is.na(samples.out[[i]][(4*n2+1):(5*n2), cols]))) samples.out[[i]][(4*n2+1):(5*n2), cols] <- samples.out[[i]][(5*n2+1):(6*n2), cols] <- samples.out[[i]][(3*n2+1):(4*n2), cols]
				if(all(is.na(samples.out[[i]][(7*n2+1):(8*n2), cols]))) samples.out[[i]][(7*n2+1):(8*n2), cols] <- samples.out[[i]][(8*n2+1):(9*n2), cols] <- samples.out[[i]][(6*n2+1):(7*n2), cols]
				if(all(is.na(samples.out[[i]][(10*n2+1):(11*n2), cols]))) samples.out[[i]][(10*n2+1):(11*n2), cols] <- samples.out[[i]][(11*n2+1):(12*n2), cols] <- samples.out[[i]][(9*n2+1):(10*n2), cols]
			}
		print(length(samples.out)-i)
		}
		save(samples.out, samples.names, region.names.out, n.samples, file=paste0("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/regional/samples/",models[1],"_",models[4],"_annual_regions_samples.RData"))
	}
	if(k==2 & !is.null(cities)){
		for(i in 1:length(stats)){
			results.years <- results[[i]][[1]][[length(results[[i]][[1]])-1]]
			for(j in 1:(length(stats[[i]])-len)) stats[[i]][[j]] <- stats[[i]][[j]][, which(results.years %in% years) ] # matrix (multiple cities)
		}
		stats.out <- list(cities=sapply(stats, function(x) t(x[[1]])))
		for(i in 1:length(stats.out)){
			rownames(stats.out[[i]]) <- rep( rep(seq(yr1, yr2), each=12), nrow(cities) )
			na.rows1 <- which(is.na(stats.out[[i]][,2]))
			na.rows2 <- which(is.na(stats.out[[i]][,5]))
			stats.out[[i]][na.rows1,2] <- stats.out[[i]][na.rows1,3] <- stats.out[[i]][na.rows1,1]
			stats.out[[i]][na.rows2,5] <- stats.out[[i]][na.rows2,6] <- stats.out[[i]][na.rows2,4]
			stats.out[[i]][na.rows1,8] <- stats.out[[i]][na.rows1,9] <- stats.out[[i]][na.rows1,7]
			stats.out[[i]][na.rows2,11] <- stats.out[[i]][na.rows2,12] <- stats.out[[i]][na.rows2,10]
		}
		d <- stats.out[[1]]
		save(d, d.cities, results.years, file=paste0("/workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/data/cities/",models[1],"_",models[4],"_annual_cities_batch", cities.batch, "_", domain, ".RData"))
	}
}
