# @knitr setup
suppressMessages(library(rgdal))
suppressMessages(library(raster))
suppressMessages(library(data.table))
rasterOptions(chunksize=10^12,maxmemory=10^11)

# @knitr getFireStats
getFireStats <- function(i, mainDir, years=NULL, cells.list, shp.names.list, n, ...){
	rep.lab <- paste0("_",c(0:199),"_")[i]
	pat <- gsub("expression","",paste(bquote(expression("^FireS.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesF <- list.files(mainDir, pattern=pat, full=T)
	patV <- gsub("expression","",paste(bquote(expression("^V.*.",.(rep.lab),".*.tif$")),collapse=""))
	filesV <- list.files(mainDir, pattern=patV, full=T)
	files.years <- as.numeric(gsub("FireScar_\\d+_", "", gsub(".tif", "", basename(filesF))))
	ord <- order(files.years)
	filesF <- filesF[ord]
	filesV <- filesV[ord]
	if(is.numeric(years)) { p <- files.years %in% years; if(any(p)) { filesF <- filesF[which(p)]; filesV <- filesV[which(p)]; files.years <- files.years[which(p)] } } # test filesV addition
	years <- unique(files.years)
	grp.names.vec <- rep(names(cells.list), times=sapply(cells.list, length))
	loc.names.vec <- as.character(unlist(lapply(cells.list, names)))
	
	getABFC <- function(x, values, index){
		#x <- intersect(index, x)
		if(length(x)){
			AB <- length(which(values[match(x, index)] > 0))
			FC <- length(unique(values[match(x, index)]))
		} else AB <- FC <- 0
		c(AB, FC)
	}
	
	# fire size by vegetation class
	fsByVeg <- function(cells, i, v, v.fid, f.ind){
		if(!length(cells)) return(data.table())
		dlist <- vector("list", length(i))
		ind <- match(cells, f.ind)
        if(!length(ind)) return(data.table())
		if(length(ind) < length(f.ind)) v[-ind] <- NA
		for(k in 1:length(i)){
			x <- v
			x[x!=i[k]] <- NA
			x <- v.fid[!is.na(x) & !is.na(v.fid)]
			if(length(x)){
				x <- tapply(x, x, length)
				ord <- order(as.numeric(x))
				dlist[[k]] <- data.table(Vegetation=i[k], FS=as.integer(x)[ord], FID=as.integer(names(x))[ord])
			}
		}
		rbindlist(dlist)
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
		vid <- 1:5
		if(!all(is.na(v.fid))){
			system.time( dl <- rapply(cells.rmNA, f=fsByVeg, classes="integer", how="replace", i=vid, v=v.veg.rmNA, v.fid=v.fid.rmNA, f.ind=f.ind) )
            #for(jj in 1:1) dl[jj][which(!sapply(dl[[jj]], is.data.table))] <- data.frame(a=1) # set NULL to empty data table
			n.fires <- unlist(sapply(1:length(dl), function(x, l) sapply(l[[x]], function(a) if(is.data.table(a)) nrow(a) else 0), l=dl))
			d.j <- rbindlist(lapply(dl, rbindlist))
			dlist[[j]] <- data.table(LocGroup=rep(grp.names.vec, times=n.fires), Location=rep(loc.names.vec, times=n.fires), VegID=d.j$Vegetation, FID=d.j$FID, Val=d.j$FS, Year=as.integer(years[j]), Replicate=as.integer(i)) # end test lines
		}
		
		# basic aggregate burn area and fire frequency
		if(j==1) m <- matrix(NA,length(filesF),2*n) # create matrix for multiple subdomains
		v <- rapply(cells.rmNA, f=getABFC, classes="integer", how="replace", values=v.fid.rmNA, index=f.ind)
        v <- unlist(sapply(1:length(v), function(x, l) sapply(l[[x]], function(a) if(is.integer(a)) a else c(0, 0)), l=v))
		m[j,1:(2*n)] <- v
		
		print(paste0("Replicate ", i, ": ", years[j]))
	}
	d.fs <- rbindlist(dlist)

	#d.fs$Vegetation <- v.names[d.fs.veg$Vegetation]
	varid <- rep(c("AB", "FC"), n)
	ab <- as.integer(m[, which(varid=="AB")])
	fc <- as.integer(m[, which(varid=="FC")])
	m <- data.table(LocGroup=rep(grp.names.vec, each=length(years)), Location=rep(loc.names.vec, each=length(years)), Var=rep(c("Burn Area", "Fire Count"), each=length(ab)), Val=c(ab, fc), Year=as.integer(years), Replicate=as.integer(i))
	print(paste("Returning list of (1) area burned and fire frequency data table and (2) fire size by vegetation class data table."))
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
    
    # Bootstrapping
    btfun <- function(p, n.samples=length(p)/2, n.boot=10000, interp=FALSE, n.interp=1000, ...){
        if(!length(p)) return(p)
        if(all(is.na(p))) return(rep(NA, n.boot))
        p <- list(x=p[1:n.samples], y=p[(n.samples+1):(2*n.samples)])
        if(interp && length(unique(p[1:n.samples])) > 1) p <- approx(p$x, p$y, n=n.interp)
        p <- round(sample(p$x, n.boot, prob=p$y, rep=T), ...)
        p
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
                if(!length(cells.tmp)) next
				a.tmp <- a[cells.tmp]
				v.tmp <- v[cells.tmp]
				v.tmp.n <- tapply(cells.tmp, v.tmp, length)
				den <- unlist(tapply(a.tmp, v.tmp, denFun, n=n.samples))
                samp <- tapply(den, rep(1:(length(den)/(2*n.samples)), each=2*n.samples), btfun, n.samples=n.samples, n.boot=n.samples, interp=TRUE)
				dlist2[[l2]] <- data.table(LocGroup=names(x)[l1], Location=names(x[[l1]])[l2], VegID=rep(as.integer(names(samp)), each=n.samples), VegArea=rep(v.tmp.n, each=n.samples), Val=unlist(samp), Year=as.integer(years[j]))
			}
			dlist1[[l1]] <- rbindlist(dlist2)
		}
		dlist[[j]] <- rbindlist(dlist1)
        if(!is.data.table(dlist[[j]])){
            print(paste("PROBLEM WITH REPLICATE", i, "YEAR", years[j]))
            print(dlist[[j]])
            print(paste0("Class of dlist[[", j, "]]:", class(dlist[[j]])))
            stop(paste("PROBLEM WITH REPLICATE", i, "YEAR", years[j]))
        }
		print(paste("Replicate", i, "age and veg maps:", years[j]))
	}
	
	x <- rbindlist(dlist)
    
	x$Replicate <- as.integer(i)
	veg.area <- x[seq(1, nrow(x), by=2*n.samples)][, Val:=NULL]
	setnames(veg.area, "VegArea", "Val")
	veg.area[, Var := "Vegetated Area"]
	setcolorder(veg.area, c("LocGroup", "Location", "VegID", "Var", "Val", "Year", "Replicate"))
	locs <- unique(x$Location)
	x[, VegArea := NULL]
	x[, Replicate := NULL]
	setkey(x, Location)
	for(j in 1:length(locs)){
		obj.name.tmp <- paste0("age__", locs[j], "__rep", i)
		assign(obj.name.tmp, x[locs[j]])
		save(list=c("locs", obj.name.tmp), file=paste0(denDir, "/", obj.name.tmp, ".RData"))
		print(paste(obj.name.tmp, "object", j, "of", length(locs), "saved."))
		rm(list=obj.name.tmp)
		gc()
	}
	rm(x)
	gc()
	print("Returning veg class areas data frame for subregions.")
    if(!is.data.table(veg.area)){
        print(paste("PROBLEM WITH REPLICATE", i))
        print(veg.area)
        print(paste("Class of veg.area:", class(veg.area)))
        stop(paste("PROBLEM WITH REPLICATE", i))
    }
	veg.area
}
