# @knitr setup
library(raster)
library(data.table)
rasterOptions(chunksize=10^12,maxmemory=10^11)

# @knitr getFireStats
getFireStats <- function(i, mainDir, years=NULL, cells.list, shp.names.list, n, ...){
	rep.lab <- paste0("_",c(0:199),"_")[i]
	pat <- gsub("expression","",paste(bquote(expression("^FireS.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesF <- list.files(mainDir, pattern=pat, full=T)
	patV <- gsub("expression","",paste(bquote(expression("^V.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesV <- list.files(mainDir, pattern=patV, full=T)
	#files.years <- as.numeric(substr(files,nchar(files)-7,nchar(files)-4))
	files.years <- as.numeric(gsub("FireScar_\\d+_", "", gsub(".tif", "", basename(filesF)))) # test line
	ord <- order(files.years) # test line
	filesF <- filesF[ord] # test line
	filesV <- filesV[ord] # test line
	if(is.numeric(years)) { p <- files.years %in% years; if(any(p)) { filesF <- filesF[which(p)]; filesV <- filesV[which(p)]; files.years <- files.years[which(p)] } } # test filesV addition
	years <- unique(files.years)
	grp.names.vec <- rep(names(cells.list), times=sapply(cells.list, length))
	loc.names.vec <- as.character(unlist(lapply(cells.list, names)))
	
	getABFC <- function(x, values, index){
		#x <- intersect(index, x)
		if(length(x)){
			AB <- length(which(values[x] > 0))
			FC <- length(unique(values[x]))
		} else AB <- FC <- 0
		c(AB, FC)
	}
	
	# fire size by vegetation class
	fsByVeg <- function(cells, i, v, v.fid, f.ind){
		if(!length(cells)) return(data.frame())
		dlist <- vector("list", length(i))
		ind <- match(cells, f.ind)
		if(length(ind) < length(f.ind)) v[-ind] <- NA
		for(k in 1:length(i)){
			x <- v
			x[x!=i[k]] <- NA
			x <- v.fid[!is.na(x) & !is.na(v.fid)]
			if(length(x)){
				x <- tapply(x, x, length)
				ord <- order(as.numeric(x))
				dlist[[k]] <- data.frame(Vegetation=i[k], FS=as.numeric(x)[ord], FID=as.numeric(names(x))[ord])
			}
		}
		as.data.frame(rbindlist(dlist))
	}
	
	print("Beginning processing loop")
	dlist <- vector("list", length(years))
	for(j in 1:length(years)){
		
		# fire size by vegetation class
		v.fid <- getValues(raster(filesF[j], band=2))
		v.veg <- getValues(raster(filesV[j]))
		f.ind <- which(!is.na(v.fid))
		v.fid.rmNA <- v.fid[f.ind]
		v.veg.rmNA <- v.veg[f.ind]
		cells.rmNA <- rapply(cells.list, f=function(x, i) intersect(i, x), classes="integer", how="replace", i=f.ind)
		#v.veg[v.veg==2 | v.veg==3] <- 1 # 2 and 3 tree classes combine into class 1 to become 'forest', tundra types 4, 5, and 6 remain as before
		vid <- 1:5 #vid <- sort(unique(v.veg[!is.na(v.veg) & v.veg > 0]))
		if(!all(is.na(v.fid))){
			system.time( dl <- rapply(cells.rmNA, f=fsByVeg, classes="integer", how="replace", i=vid, v=v.veg.rmNA, v.fid=v.fid.rmNA, f.ind=f.ind) )
			n.fires <- unlist(sapply(1:length(dl), function(x, l) sapply(l[[x]], nrow), l=dl))
			d.j <- rbindlist(lapply(dl, rbindlist))
			dlist[[j]] <- data.frame(LocGroup=rep(grp.names.vec, times=n.fires), Location=rep(loc.names.vec, times=n.fires), VegID=d.j$Vegetation, FID=d.j$FID, Val=d.j$FS, Year=years[j], Replicate=i, stringsAsFactors=FALSE) # end test lines
		}
		
		# basic aggregate burn area and fire frequency
		
		if(j==1) m <- matrix(NA,length(filesF),2*n) # create matrix for multiple subdomains
		v <- rapply(cells.rmNA, f=getABFC, classes="integer", how="unlist", values=v.fid.rmNA, index=f.ind)
		m[j,1:(2*n)] <- v
		
		print(years[j])
	}
	d.fs <- as.data.frame(rbindlist(dlist))

	#d.fs$Vegetation <- v.names[d.fs.veg$Vegetation]
	varid <- rep(c("AB", "FC"), n)
	ab <- as.integer(m[, which(varid=="AB")])
	fc <- as.integer(m[, which(varid=="FC")])
	m <- data.frame(LocGroup=rep(grp.names.vec, each=length(years)), Location=rep(loc.names.vec, each=length(years)), Var=rep(c("Burn Area", "Fire Count"), each=length(ab)), Val=c(ab, fc), Year=years, Replicate=i, stringsAsFactors=FALSE)
	rownames(m) <- NULL
	#print(paste("Returning area burned and fire frequency data frame."))
	#m
	print(paste("Returning list of (1) area burned and fire frequency data frame and (2) fire size by vegetation class data frame."))
	list(m, d.fs)
}

# @knitr getAgeVegStats
getAgeVegStats <- function(i, mainDir, denDir, years=NULL, cells.list, shp.names.list, n, n.samples, ...){
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
	x <- rapply(cells.list, f=function(x) lapply(1:length(filesA), function(x) 0L), classes="integer", how="replace")
	
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
	
	print("Beginning processing loop")
	dlist <- vector("list", length(years))
	for(j in 1:length(years)){
		a <- getValues(raster(filesA[j]))[data.ind]
		v <- getValues(raster(filesV[j]))[data.ind]
		#a[a<0] <- a[a<0] + 2147483647 # temporary hack
		dlist1 <- vector("list", length(x))
		for(l1 in 1:length(x)){
			dlist2 <- vector("list", length(x[[l1]]))
			for(l2 in 1:length(x[[l1]])){
				cells.tmp <- cells.list[[l1]][[l2]]
				a.tmp <- a[cells.tmp]
				v.tmp <- v[cells.tmp]
				v.tmp.n <- tapply(cells.tmp, v.tmp, length)
				den.tmp <- tapply(a.tmp, v.tmp, denFun, n=n.samples)
				dlist2[[l2]] <- data.frame(LocGroup=names(x)[l1], Location=names(x[[l1]])[l2], VegID=rep(as.numeric(names(den.tmp)), each=2*n.samples), VegArea=rep(v.tmp.n, each=2*n.samples), Val=unlist(den.tmp), Year=years[j], stringsAsFactors=FALSE)
			}
			dlist1[[l1]] <- rbindlist(dlist2)
		}
		dlist[[j]] <- rbindlist(dlist1)
		print(paste("Age and veg maps:", years[j]))
	}
	
	x <- rbindlist(dlist)
	x$Replicate <- i
	veg.area <- x[seq(1, nrow(x), by=2*n.samples), -5]
	names(veg.area)[4] <- "Val"
	veg.area$Var <- "Vegetated Area"
	veg.area <- veg.area[,c(1,2,7,4,3,5,6)]
	locs <- unique(x$Location)
	for(j in 1:length(locs)){
		obj.name.tmp <- paste0("age__", locs[j], "__rep", i)
		assign(obj.name.tmp, x[x$Location==locs[j], -c(4,7)])
		save(list=c("locs", obj.name.tmp), file=paste0(denDir, "/", obj.name.tmp, ".RData"))
		print(paste(obj.name.tmp, "object of", length(locs), "saved."))
		rm(list=obj.name.tmp)
		gc()
	}
	rm(x)
	gc()
	print("Returning veg class areas data frame for subregions.")
	veg.area
}

# @knitr getAlfStats
getAlfStats <- function(i, mainDir, denDir, years=NULL, cells.list, shp.names.list, n, n.samples, replace.deciduous=FALSE, ...){
	rep.lab <- paste0("_",c(0:199),"_")[i]
	patF <- gsub("expression","",paste(bquote(expression("^FireS.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesF <- list.files(mainDir, pattern=patF, full=T)
	patA <- gsub("expression","",paste(bquote(expression("^A.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesA <- list.files(mainDir, pattern=patA, full=T)
	patV <- gsub("expression","",paste(bquote(expression("^V.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesV <- list.files(mainDir, pattern=patV, full=T)
	files.years <- as.numeric(substr(filesA,nchar(filesA)-7,nchar(filesA)-4))
	if(is.numeric(years)) { p <- files.years %in% years; if(any(p)) { filesF <- filesF[p]; filesA <- filesA[p]; filesV <- filesV[p]; files.years <- files.years[p] } }
	years <- unique(files.years)
	print("Vectors of annual files obtained for fire, age, and veg.")
	grp.names.vec <- rep(names(cells.list), times=sapply(cells.list, length))
	loc.names.vec <- as.character(unlist(lapply(cells.list, names)))
	f.var <- rep(c("AB", "FC"), n)
	x <- rapply(cells.list, f=function(x) lapply(1:length(filesA), function(x) 0L), classes="integer", how="replace")
	# Area burned and fire count function
	getABFC <- function(x, values, index){
		x <- intersect(index, x)
		if(length(x)){
			AB <- length(which(values[x] > 0))
			FC <- length(unique(values[x]))
		} else AB <- FC <- 0
		c(AB, FC)
	}
	# Density estimation function
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
	print("Support functions defined and setup completed. Beginning annual processing...")
	for(j in 1:length(years)){
		ff <- getValues(raster(filesF[j], band=2))
		f.ind <- which(!is.na(ff))
		ff <- rapply(cells.list, f=getABFC, classes="integer", how="unlist", values=ff, index=f.ind)
		if(j==1){
			fire <- matrix(NA,length(filesF),2*n) # create an area burned and fire count matrix for multiple subdomains
			a <- getValues(raster(filesA[1])) # use as a template
			a.ind <- which(!is.na(a))
			a <- a[a.ind]
		} else {
			a <- a + 1
			af.ind <- intersect(a.ind, f.ind)
			if(length(af.ind)) a[af.ind] <- 0
		}
		fire[j,1:(2*n)] <- ff
		v <- getValues(raster(filesV[j]))[a.ind]
		for(l1 in 1:length(x)){
			for(l2 in 1:length(x[[l1]])){
				cells.tmp <- cells.list[[l1]][[l2]]
				a.tmp <- a[cells.tmp]
				v.tmp <- v[cells.tmp]
				v.tmp.n <- tapply(cells.tmp, v.tmp, length)
				den.tmp <- tapply(a.tmp, v.tmp, denFun, n=n.samples)
				x[[l1]][[l2]][[j]] <- data.frame(LocGroup=names(x)[l1], Location=names(x[[l1]])[l2], VegID=rep(as.numeric(names(den.tmp)), each=2*n.samples), VegArea=rep(v.tmp.n, each=2*n.samples), Val=unlist(den.tmp), Year=years[j], stringsAsFactors=FALSE)
			}
		}
		print(paste("Fire, age and vegetation processed:", years[j]))
	}
	for(l1 in 1:length(x)) for(l2 in 1:length(x[[l1]])) x[[l1]][[l2]] <- do.call(rbind, x[[l1]][[l2]])
	x <- do.call(rbind, lapply(x, do.call, what=rbind))
	x$Replicate <- i
	veg.area <- x[seq(1, nrow(x), by=2*n.samples), -5]
	names(veg.area)[4] <- "Val"
	veg.area$Var <- "Vegetated Area"
	veg.area <- veg.area[,c(1,2,7,4,3,5,6)]
	rownames(veg.area) <- NULL
	print(paste("Veg class areas data frame compiled. Writing age density data frames to R workspaces..."))
	locs <- unique(x$Location)
	for(j in 1:length(locs)){
		obj.name.tmp <- paste0("veg.age.loc", j, ".rep", i)
		assign(obj.name.tmp, x[x$Location==locs[j], -c(4,7)])
		save(list=c("locs", obj.name.tmp), file=paste0(denDir, "/", obj.name.tmp, ".RData"))
		print(paste(obj.name.tmp, "object of", length(locs), "saved."))
		rm(list=obj.name.tmp)
		gc()
	}
	rm(x)
	gc()
	print(paste("Age densities saved for all locations. Compiling area burned and fire frequency data frame..."))
	ab <- as.integer(fire[, which(f.var=="AB")])
	fc <- as.integer(fire[, which(f.var=="FC")])
	fire <- data.frame(LocGroup=rep(grp.names.vec, each=length(years)), Location=rep(loc.names.vec, each=length(years)), Var=rep(c("Burn Area", "Fire Count"), each=length(ab)), Val=c(ab, fc), Year=years, Replicate=i, stringsAsFactors=FALSE)
	rownames(fire) <- NULL
	print("Returning list of (1) area burned and fire frequency data frame and (2) veg class areas data frame for all regions.")
	list(fire=fire, veg.area=veg.area)
}
