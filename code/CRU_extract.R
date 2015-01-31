# @knitr clargs
comArgs <- commandArgs(TRUE)
if(length(comArgs)) for(i in 1:length(comArgs)) eval(parse(text=comArgs[[i]]))

if(any(comArgs=="akcan2km")) domain <- "akcan2km" else if(any(comArgs=="world10min")) domain <- "world10min" else stop("Unquoted 'akcan2km' or 'world10min' argument not supplied.")
if(any(comArgs=="regions")) regions <- TRUE else regions <- FALSE
if(any(comArgs=="cities")) cities <- TRUE else cities <- FALSE
if(!(regions | cities)) stop("Unquoted 'regions' or 'cities' argument supplied. Must provide at least one.")

# @knitr packages
library(raster)
library(maptools)
library(parallel)
library(plyr)

# @knitr source
if(domain=="akcan2km"){ # For regions and/or cities
	topDir <- file.path("/Data/Base_Data/Climate/AK_CAN_2km",c("historical/singleBand/CRU/cru_TS31/historical"))
	load("/workspace/UA/mfleonawicz/Leonawicz/Projects/2014/DataExtraction/workspaces/shapes2cells_AKCAN2km_5pct.RData")
	locs <- read.csv("/workspace/Shared/Users/mfleonawicz/github/statistics/AR5_scripts/AR5_QAQC/locs.csv")
} else if(domain=="world10min") { # Currently for cities only
	topDir <- file.path("/Data/Base_Data/Climate/World/World_10min",c("historical/CRU/CRU_TS31")) # files are not read, but metadata parsed from filenames list
	#### Need to insert a load() command for a locs object analogous to that above
}

# @knitr setup
if(exists("years")) yr1 <- years[1] else yr1 <- 1901
if(exists("years")) yr2 <- tail(years, 1) else yr2 <- 2009
years <- yr1:yr2
varid <- c("tas","pr")

#if(!is.null(cities)){
#	cities <- subset(locs, pop >= 1000, c("lon_albers","lat_albers"))
#	d.cities <- subset(locs, pop >= 1000)[-c(which(names(locs) %in% c("lon_albers","lat_albers")))]
#}
if(cities){
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
	d.cities <- cities[-c(which(names(locs) %in% c("lon_albers","lat_albers")))]
	cities <- if(domain=="akcan2km") cbind(cities$lon_albers, cities$lat_albers) else if(domain=="world10min") cbind(cities$lon, cities$lat)
}

seasons <- "annual"
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
calcMeans <- function(i, varid, cells.list, shp.names, n.shp, seed=232, regions=TRUE, n.samples=512, cities=NULL, start.year=NULL, end.year=NULL){
	print(i)
	for(p in 1:length(varid)){
		path <- file.path(topDir, varid[p])
		files <- list.files(path,pattern=".tif$",full=T)
		if(!is.null(start.year)) files <- files[substr(files, nchar(files)-7, nchar(files)-4) >= start.year]
		if(!is.null(end.year)) files <- files[substr(files, nchar(files)-7, nchar(files)-4) <= end.year]

		if(!is.null(cities)){
			r <- readAll(raster(files[1])) # template done
			cells_cities <- extract(r, cities, cellnumbers=T)[,1]
			rank_cells_cities <- rank(cells_cities)
			cells_cities <- sort(cells_cities)
			print("Raster cell indices for point locations obtained.")
		}

		mo.tmp <- substr(files,nchar(files)-10,nchar(files)-9)
		yr.tmp <- substr(files,nchar(files)-7,nchar(files)-4)
		yr.mo.tmp <- paste0(yr.tmp,mo.tmp)
		files <- files[order(yr.mo.tmp)]
		yr.tmp <- substr(files,nchar(files)-7,nchar(files)-4)
		seq.q <- c(0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95) #seq(0, 1, by=0.05) # Hard coded
		m <- matrix(NA, nrow=length(files), ncol=length(seq.q) + 2) # +2 for mean and SD below
		
		mden <- matrix(NA, nrow=2*n.samples, ncol=length(files))
		samples.list <- rapply(cells.list, f=function(x, m) m, classes="integer", how="replace", m=mden)
		names(samples.list) <- names(shp.names)
		for(l1 in 1:length(shp.names)) names(samples.list[[l1]]) <- shp.names[[l1]]
		
		if(regions) for(l1 in 1:length(shp.names)) for(l2 in 1:length(shp.names[[l1]])) assign(paste("m", names(shp.names)[l1], shp.names[[l1]][l2], sep="__"), m)
		rm(m)
		gc()
		if(!is.null(cities)) m.cities <- matrix(NA, nrow=length(cells_cities), ncol=length(files))
		for(b in 1:length(unique(yr.tmp))){
			pat <- gsub("expression","",paste(bquote(expression(".",.(unique(yr.tmp)[b]),".tif$")),collapse=""))
			files.sub <- list.files(path, pattern=pat, full=T)
			mat <- getValues(stack(files.sub, quick=T))
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
						samples.tmp <- apply(mat[cells.tmp,], 2, denFun, n=n.samples, variable=varid[p]) # Time across columns
						samples.list[[ names(shp.names)[l1] ]][[ shp.names[[l1]][l2] ]][, (1:n)+(12*(b-1))] <- samples.tmp
					}
				}
				gc()
			}
			if(!is.null(cities)){
				m.cities[,(1:n)+(12*(b-1))] <- mat[cells_cities[rank_cells_cities],]
			}
			print(paste0("Process",i,": ",unique(yr.tmp)[b]))
		}
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
	if(!regions & !is.null(cities))	l <- list(list(yr=yr.hold), list(m.cities=m.cities.hold))
	l
}

# @knitr run
results <- mclapply(X=1, FUN=calcMeans, j=1, varid=varid, cells.list=cells_shp_list_5pct, shp.names=region.names.out, n.shp=n.shp,
	regions=regions, n.samples=n.samples, cities=cities, start.year=yr1, end.year=yr2, mc.cores=1)

# @knitr save
# Organize and save results
for(k in 1:2){
	if(!regions & k==1) next
	stats <- lapply(results, "[[", k)
	if(is.character(stats[[1]][[1]])) next
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
		time.seq <- rep(seq(years[1],length.out=nrow(stats.out[[1]])/(2*12)), each=12)
		for(i in 1:length(stats.out)){
			colnames(stats.out[[i]]) <- agg.stat.names
			rownames(stats.out[[i]]) <- rep(time.seq, 2)
		print(length(stats.out)-i)
		}
		save(stats.out, results.years, region.names.out, agg.stat.names, agg.stat.colnames,
			file=paste0("/workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/data/regional/stats/","CRU31","_",seasons[1],"_regions_stats.RData"))
		
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
		samples <- samples[[1]]
		mid <- (ncol(samples)-1)/2 + 1
		names.hold <- names(samples)[1:mid]
		samples <- data.frame(Location=samples$Location, rbind(as.matrix(samples[, 2:mid]), as.matrix(samples[, (mid+1):ncol(samples)])))
		samples$Location <- rep(rep(samples.names, each=n2))

		names(samples)[-1] <- paste(results.years[1:(length(results.years)/2)], month.abb, sep="_")
		samples.out <- list()
		for(i in 1:length(samples.names)){
			samples.out[[i]] <- subset(samples, Location==samples.names[i])
			rownames(samples.out[[i]]) <- NULL
		}
		names(samples.out) <- samples.names
		save(samples.out, samples.names, region.names.out, n.samples, file=paste0("/workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/data/regional/samples/","CRU31","_",seasons[1],"_regions_samples.RData"))
	
	}
	if(k==2 & !is.null(cities)){
		for(i in 1:length(stats)){
			results.years <- results[[i]][[1]][[length(results[[i]][[1]])]]
			for(j in 1:(length(stats[[i]])-len)) stats[[i]][[j]] <- stats[[i]][[j]][, which(results.years %in% years) ] # matrix (multiple cities)
		}
		stats.out <- list(cities=sapply(stats, function(x) t(x[[1]])))
		for(i in 1:length(stats.out)) rownames(stats.out[[i]]) <- rep( rep( rep(seq(yr1, yr2), each=12), nrow(cities) ), length(varid) )
		d <- stats.out[[1]]
		save(d, d.cities, results.years, file=paste0("/workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/data/cities/","CRU31_annual_cities_", domain, ".RData"))
	}
}
