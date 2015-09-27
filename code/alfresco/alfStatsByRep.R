# @knitr setup
suppressMessages(library(rgdal))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
rasterOptions(chunksize=10^12,maxmemory=10^11)

prep_alf_data <- function(dir, repid, years){
    rep.lab <- 
    p <- lapply(c("A", "V", "FireS"), function(p, r) gsub("expression","",paste(bquote(expression("^",.(p),".*.",.(paste0("_",c(0:199),"_")[r]),".*.tif$")),collapse="")), r=repid)
    files <- lapply(1:length(plist), function(i, dir, p) list.files(dir, pattern=p[[i]], full=TRUE), dir=dir, p=p)
    names(files) <- c("Age", "Veg", "FID")
    files.years <- as.numeric(gsub("FireScar_\\d+_", "", gsub(".tif", "", basename(files$FID))))
    ord <- order(files.years)
    files <- lapply(files, function(x, ord) x[ord], ord=ord)
    files <- lapply(files, function(x, file.yrs, keep.yrs) { k <- file.yrs %in% keep.yrs; if(any(k)) x[which(k)] else x }, file.yrs=files.years, keep.yrs=years)
    stopifnot(all.equal(length(files$Age), length(files$Veg), length(files$FID), length(years)))
    c(files, years=as.integer(years))
}

# @knitr getFireStats
getFireStats <- function(i, mainDir, years=NULL, cells, ...){
    if(is.null(list(...)$veg.labels)) veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra") else veg.labels <- list(...)$veg.labels
    x <- prep_alf_data(dir=mainDir, repid=i, years=years)
    cells <- group_by(cells, LocGroup, Location)
    d.fs <- vector("list", length(years))
    for(j in 1:length(x$years)){ # fire size by vegetation class
        v <- list(FID=getValues(raster(x$FID[j], band=2)), Veg=getValues(raster(x$Veg[j])))
        d.fs[[j]] <- filter(cells, Cell %in% which(!is.na(v$FID))) %>% mutate(Vegetation=veg.labels[v$Veg[Cell]], FID=v$FID[Cell]) %>% group_by(LocGroup, Location, Vegetation, FID) %>% summarise(FS=length(Cell), Year=x$years[j])
        print(paste0("Replicate ", i, ": ", x$years[j]))
    }
    d.fs <- rbindlist(d.fs)[, Replicate:=as.integer(i)]
    print(paste("Returning fire size by vegetation class data table."))
    d.fs
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
        a[a<0] <- a[a<0] + 2147483647 # temporary hack
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
    veg.area <- x[seq(1, nrow(x), by=n.samples)][, Val:=NULL]
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
