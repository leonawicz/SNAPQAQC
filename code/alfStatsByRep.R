# @knitr setup
library(raster)
rasterOptions(chunksize=10^12,maxmemory=10^11)

# @knitr getFireStats
getFireStats <- function(i, mainDir, years=NULL, cells.list, shp.names.list, n, ...){
	rep.lab <- paste0("_",c(0:199),"_")[i]
	pat <- gsub("expression","",paste(bquote(expression("^FireS.*.",.(rep.lab),".*.tif$")),collapse=""))
	files <- list.files(mainDir, pattern=pat, full=T)
	files.years <- as.numeric(substr(files,nchar(files)-7,nchar(files)-4))
	if(is.numeric(years)) { p <- files.years %in% years; if(any(p)) { files <- files[which(p)]; files.years <- files.years[which(p)] } }
	years <- unique(files.years)
	grp.names.vec <- rep(names(cells.list), times=sapply(cells.list, length))
	loc.names.vec <- as.character(unlist(lapply(cells.list, names)))
	
	getABFC <- function(x, values, index){
		x <- intersect(index, x)
		if(length(x)){
			AB <- length(which(values[x] > 0))
			FC <- length(unique(values[x]))
		} else AB <- FC <- 0
		c(AB, FC)
	}
	
	print("Beginning processing loop")
	for(j in 1:length(years)){
		rvals <- getValues(raster(files[j], band=2))
		f.ind <- which(!is.na(rvals))
		if(j==1) m <- matrix(NA,length(files),2*n) # create matrix for multiple subdomains
		v <- rapply(cells.list, f=getABFC, classes="integer", how="unlist", values=rvals, index=f.ind)
		m[j,1:(2*n)] <- v
		print(years[j])
	}
	varid <- rep(c("AB", "FC"), n)
	ab <- as.integer(m[, which(varid=="AB")])
	fc <- as.integer(m[, which(varid=="FC")])
	m <- data.frame(LocGroup=rep(grp.names.vec, each=length(years)), Location=rep(loc.names.vec, each=length(years)), Var=rep(c("Burn Area", "Fire Count"), each=length(ab)), Val=c(ab, fc), Year=years, Replicate=i, stringsAsFactors=FALSE)
	rownames(m) <- NULL
	print(paste("Returning area burned and fire frequency data frame."))
	m
}

# @knitr getAgeVegStats
getAgeVegStats <- function(i, mainDir, denDir, years=NULL, cells.list, shp.names.list, n, breaks, age.lab, veg.lab, n.samples, replace.deciduous=FALSE, ...){
	id.vals <- as.numeric(paste0(rep(1:(length(age.lab)), length(veg.lab)), rep(1:length(veg.lab), each=length(age.lab))))
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
	denFun <- function(x, n, min.zero=TRUE, diversify=FALSE){
		x <- x[!is.na(x)]
		lx <- length(x)
		if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
		if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
		dif <- diff(range(x))
		z <- density(x, adjust=2, n=n, from=min(x)-max(1, 0.05*dif), to=max(x)+max(1, 0.05*dif))
		if(min.zero && any(z$x < 0)) z <- density(x, adjust=2, n=n, from=0, to=max(x)+max(1, 0.05*dif))
		as.numeric(c(z$x, z$y))
	}
	
	for(j in 1:length(years)){
		a <- getValues(raster(filesA[j]))[data.ind]
		v <- getValues(raster(filesV[j]))[data.ind]
		#a[a<0] <- a[a<0] + 2147483647 # temporary hack
		for(l1 in 1:length(x)){
			for(l2 in 1:length(x[[l1]])){
				cells.tmp <- cells.list[[l1]][[l2]]
				a.tmp <- a[cells.tmp]
				v.tmp <- v[cells.tmp]
				#x[[l1]][[l2]][[j]] <- table(10*.bincode(a.tmp, breaks=breaks) + v.tmp)
				#uni.veg.tmp <- sort(unique(v.tmp))
				v.tmp.n <- tapply(cells.tmp, v.tmp, length)
				den.tmp <- tapply(a.tmp, v.tmp, denFun, n=n.samples)
				x[[l1]][[l2]][[j]] <- data.frame(LocGroup=names(x)[l1], Location=names(x[[l1]])[l2], VegID=rep(as.numeric(names(den.tmp)), each=2*n.samples), VegArea=rep(v.tmp.n, each=2*n.samples), Val=unlist(den.tmp), Year=years[j], stringsAsFactors=FALSE)
			}
		}
		print(paste("Age and veg maps:", years[j]))
	}
	
	#makeTable <- function(l, id, years, age, veg, loc, loc.grp){
	#	nr <- length(id)
	#	nc <- length(l)
	#	n.age <- length(age)
	#	n.veg <- length(veg)
	#	d <- matrix(0,nrow=nr,ncol=nc)
	#	x <- lapply(l, names)
	#	for(k in 1:nc) d[match(x[[k]], id), k] <- l[[k]]
	#	d <- data.frame(Veg=rep(veg, each=n.age), Age=rep(age, n.veg), LocGroup=loc.grp, Location=loc, Area=as.integer(d), Year=rep(years, each=n.age*n.veg), stringsAsFactors=FALSE)
	#	rownames(d) <- NULL
	#	d
	#}
	
	#for(l1 in 1:length(x)) for(l2 in 1:length(x[[l1]])) x[[l1]][[l2]] <- makeTable(l=x[[l1]][[l2]], id=id.vals, years=years, age=age.lab, veg=veg.lab, loc=names(x[[l1]])[l2], loc.grp=names(x)[l1])
	for(l1 in 1:length(x)) for(l2 in 1:length(x[[l1]])) x[[l1]][[l2]] <- do.call(rbind, x[[l1]][[l2]])
	x <- do.call(rbind, lapply(x, do.call, what=rbind))
	x$Replicate <- i
	veg.area <- x[seq(1, nrow(x), by=2*n.samples), -5]
	names(veg.area)[4] <- "Val"
	veg.area$Var <- "Vegetated Area"
	veg.area <- veg.area[,c(1,2,7,4,3,5,6)]
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
	print("Returning veg class areas data frame for subregions.")
	veg.area
}

# @knitr getAlfStats
getAlfStats <- function(i, mainDir, denDir, years=NULL, cells.list, shp.names.list, n, breaks, age.lab, veg.lab, n.samples, replace.deciduous=FALSE, ...){
	id.vals <- as.numeric(paste0(rep(1:(length(age.lab)), length(veg.lab)), rep(1:length(veg.lab), each=length(age.lab))))
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
	denFun <- function(x, n, min.zero=TRUE, diversify=FALSE){
		x <- x[!is.na(x)]
		lx <- length(x)
		if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
		if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
		dif <- diff(range(x))
		z <- density(x, adjust=2, n=n, from=min(x)-max(1, 0.05*dif), to=max(x)+max(1, 0.05*dif))
		if(min.zero && any(z$x < 0)) z <- density(x, adjust=2, n=n, from=0, to=max(x)+max(1, 0.05*dif))
		as.numeric(c(z$x, z$y))
	}
	print("Support functions defined and setup completed. Beginning annual prcoessing...")
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
