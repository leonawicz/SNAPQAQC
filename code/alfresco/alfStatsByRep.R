# @knitr setup
suppressMessages(library(rgdal))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
rasterOptions(chunksize=10^12,maxmemory=10^11)

prep_alf_data <- function(dir, repid, years){
    rep.lab <- 
    p <- lapply(c("A", "V", "FireS"), function(p, r) gsub("expression","",paste(bquote(expression("^",.(p),".*.",.(paste0("_",c(0:199),"_")[r]),".*.tif$")),collapse="")), r=repid)
    files <- lapply(1:length(p), function(i, dir, p) list.files(dir, pattern=p[[i]], full=TRUE), dir=dir, p=p)
    names(files) <- c("Age", "Veg", "FID")
    files.years <- as.numeric(gsub("FireScar_\\d+_", "", gsub(".tif", "", basename(files$FID))))
    ord <- order(files.years)
    files <- lapply(files, function(x, ord) x[ord], ord=ord)
    files <- lapply(files, function(x, file.yrs, keep.yrs) { k <- file.yrs %in% keep.yrs; if(any(k)) x[which(k)] else x }, file.yrs=files.years, keep.yrs=years)
    stopifnot(all.equal(length(files$Age), length(files$Veg), length(files$FID), length(years)))
    files$years <- as.integer(years)
    files
}

# @knitr getFireStats
getFireStats <- function(i, mainDir, years=NULL, cells, ...){
    if(is.null(list(...)$veg.labels)) veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra") else veg.labels <- list(...)$veg.labels
    x <- prep_alf_data(dir=mainDir, repid=i, years=years)
    cells <- group_by(cells, LocGroup, Location)
    d.fs <- vector("list", length(years))
    for(j in 1:length(x$years)){ # fire size by vegetation class
        v <- list(FID=getValues(raster(x$FID[j], band=2)), Veg=getValues(raster(x$Veg[j])))
        d.fs[[j]] <- filter(cells, Cell %in% which(!is.na(v$FID))) %>% mutate(Vegetation=veg.labels[v$Veg[Cell]], FID=v$FID[Cell]) %>%
            group_by(LocGroup, Location, Vegetation, FID) %>% summarise(FS=length(Cell), Year=x$years[j]) %>% setorder(LocGroup, Location, Vegetation, Year, FS)
        print(paste0("Replicate ", i, ": ", x$years[j]))
    }
    d.fs <- rbindlist(d.fs)[, Replicate:=as.integer(i)]
    print(paste("Returning fire size by vegetation class data table."))
    d.fs
}

# @knitr getAgeVegStats
getAgeVegStats <- function(i, mainDir, vegDir, years=NULL, cells, ...){
    if(is.null(list(...)$veg.labels)) {
        veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
    } else veg.labels <- list(...)$veg.labels
    if(!is.numeric(list(...)$age.samples)) age.samples <- 1000 else age.samples <- list(...)$age.samples
    x <- prep_alf_data(dir=mainDir, repid=i, years=years)
    cells <- group_by(cells, LocGroup, Location)
    r <- getValues(raster(x$Age[1])) # use as a template
    idx <- which(!is.na(r))
    idx.rmNA <- which(idx %in% 1:length(r))
    # Density estimation
    denFun <- function(x, n=1000, adj=0.1, out="vector", min.zero=TRUE, diversify=FALSE){
        x <- x[!is.na(x)]
        lx <- length(x)
        if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
        if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
        b <- max(1, 0.05*diff(range(x)))
        z <- density(x, adjust=adj, n=n, from=min(x)-b, to=max(x)+b)
        if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+b)
        if(out=="vector") return(as.numeric(c(z$x, z$y))) else if(out=="list") return(z)
    }
    
    d.age <- vector("list", length(x$years))
    for(j in 1:length(x$years)){
        v <- list(Age=getValues(raster(x$Age[j]))[idx], Veg=getValues(raster(x$Veg[j]))[idx])
        v$Age[v$Age < 0] <- v$Age[ v$Age < 0] + 2147483647 # temporary hack
        d.age[[j]] <- filter(cells, Cell_rmNA %in% idx.rmNA) %>% mutate(Vegetation=veg.labels[v$Veg[Cell_rmNA]], Age=v$Age[Cell_rmNA]) %>%
            group_by(LocGroup, Location, Vegetation, Age) %>% summarise(Area=length(Cell_rmNA), Year=x$years[j])
        print(paste0("Replicate ", i, ": ", x$years[j]))
    }
    d.age <- rbindlist(d.age)[, Replicate:=as.integer(i)]
    d.age <- group_by(d.age, Replicate, LocGroup, Location, Year, Vegetation) %>% setorder(Replicate, LocGroup, Location, Year, Vegetation, Age, Area)
    #d.age <- summarise(d.age, Val=denFun(rep(Age, times=Area), n=age.samples, out="list")$x, Prob=denFun(rep(Age, times=Area), n=age.samples, out="list")$y)
    #d.area <- summarise(d.age, Val=sum(Area))
    #d.area[, Var:= "Vegetated Area"]
    #setcolorder(d.area, c("LocGroup", "Location", "Vegetation", "Var", "Val", "Year", "Replicate"))
    locs <- unique(d.age$Location)
    setkey(d.age, Location)
    for(j in 1:length(locs)){
        obj.name.tmp <- paste0("age__", locs[j], "__rep", i)
        assign(obj.name.tmp, d.age[locs[j]])
        save(list=c("locs", obj.name.tmp), file=paste0(vegDir, "/", obj.name.tmp, ".RData"))
        print(paste(obj.name.tmp, "object", j, "of", length(locs), "saved."))
        rm(list=obj.name.tmp)
        gc()
    }
    rm(d.age)
    gc()
    #print("Returning veg class areas data table for subregions.")
    #if(!is.data.table(d.area)){
    #    print(paste("PROBLEM WITH REPLICATE", i))
    #    print(d.area)
    #    print(paste("Class of d.area:", class(d.area)))
    #    stop(paste("PROBLEM WITH REPLICATE", i))
    #}
    #d.area
    return()
}
