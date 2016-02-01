# @knitr setup
suppressMessages(library(rgdal))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
rasterOptions(chunksize=10^12,maxmemory=10^11)

prep_alf_files <- function(i, loopBy, mainDir, reps, years){
    if(is.null(years)) years <- 2008:2100
    if(loopBy=="rep"){
        iter <- reps
        keep <- reps - 1
        id <- paste0("_", years[i], ".tif$")
        mainDir <- file.path(mainDir, years[i])
    } else if(loopBy=="year") {
        keep <- iter <- years
        id <- paste0("_",c(0:199)[i],"_.*.tif$")
        mainDir <- file.path(mainDir, reps[i])
    }
    p <- lapply(c("A", "V", "FireS"), function(p, id) gsub("expression","",paste(bquote(expression("^",.(p),".*.",.(id))),collapse="")), id=id)
    files <- lapply(1:length(p), function(i, dir, p) list.files(dir, pattern=p[[i]], full=TRUE, recur=TRUE), dir=mainDir, p=p)
    names(files) <- c("Age", "Veg", "FID")
    if(loopBy=="rep") files.idx <- as.numeric(gsub("FireScar_", "", gsub("_\\d+\\.tif", "", basename(files$FID)))) + 1
    if(loopBy=="year") files.idx <- as.numeric(gsub("FireScar_\\d+_", "", gsub(".tif", "", basename(files$FID))))
    ord <- order(files.idx)
    files <- lapply(files, function(x, ord) x[ord], ord=ord)
    if(loopBy=="rep") files <- lapply(files, function(x, idx) x[idx], idx=reps)
    if(loopBy=="year") files <- lapply(files, function(x, file.idx, keep) { k <- file.idx %in% keep; if(any(k)) x[which(k)] else x }, file.idx=files.idx, keep=keep)
    files$iter <- if(is.null(iter)) files.idx[ord] else iter
    stopifnot(!any(unlist(sapply(files, is.na))))
    stopifnot(all(diff(unlist(sapply(files, length)))==0))
    files
}

# @knitr extract_data
extract_data <- function(i, type, loopBy, mainDir, ageDir=NULL, reps=NULL, years=NULL, cells, readMethod="loop", ...){
    stopifnot(length(type) > 0 && type %in% c("av", "fsv"))
    if(type=="av") return(extract_av(i, loopBy, mainDir, ageDir, reps, years, cells, readMethod, ...))
    if(type=="fsv") return(extract_fsv(i, loopBy, mainDir, reps, years, cells, readMethod, ...))
}

# @knitr extract_fsv
extract_fsv <- function(i, loopBy, mainDir, reps=NULL, years=NULL, cells, readMethod="loop", ...){
    if(is.null(list(...)$veg.labels)) {
        veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
    } else veg.labels <- list(...)$veg.labels
    x <- prep_alf_files(i=i, loopBy=loopBy, mainDir=mainDir, reps=reps, years=years)
    cells <- group_by(cells, LocGroup, Location)
    d.fs <- vector("list", length(x$iter))
    
    # fire size by vegetation class
    if(readMethod=="stack"){
        print("#### READING FULL REPLICATE STACK TO MEMORY ####")
        v <- list(FID=getValues(stack(x$FID, bands=2)), Veg=getValues(stack(x$Veg, quick=TRUE)))
        gc()
        print("#### READ COMPLETE ####")
        for(j in 1:length(x$iter)){ # fire size by vegetation class
            d <- filter(cells, Cell %in% which(!is.na(v$FID[,j]))) %>% mutate(Vegetation=factor(veg.labels[v$Veg[Cell,j]], levels=veg.labels), FID=v$FID[Cell,j]) %>%
                group_by(LocGroup, Location, Vegetation, FID) %>% summarise(Val=length(Cell), Var="Fire Size")
            d.fs[[j]] <- if(loopBy=="rep") mutate(d, Replicate=x$iter[j]) else if(loopBy=="year") mutate(d, Year=x$iter[j])
            print(switch(loopBy, "rep"=paste0("Year ", years[i], ": Replicate ", x$iter[j]), "year"=paste0("Replicate ", i, ": Year ", years[x$iter[j]])))
        }
        rm(v)
        gc()
    } else if(readMethod=="loop"){
        # fire size by vegetation class
        #fn <- function(j, x, v, cells, veg.labels, loopBy){
        #    d <- filter(cells, Cell %in% which(!is.na(v$FID))) %>% mutate(Vegetation=factor(veg.labels[v$Veg[Cell]], levels=veg.labels), FID=v$FID[Cell]) %>%
        #        group_by(LocGroup, Location, Vegetation, FID) %>% summarise(Val=length(Cell), Var="Fire Size")
        #    d <- if(loopBy=="rep") mutate(d, Replicate=x$iter[j]) else if(loopBy=="year") mutate(d, Year=x$iter[j])
        #    d
        #}
        #system.time({ d.fs <- lapply(1:10, fn, x, v, cells, veg.labels, loopBy) })
        
        # fire size by vegetation class
        #fn <- function(j, x, cells, veg.labels, loopBy){
        #    v <- list(FID=getValues(raster(x$FID[j], band=2)), Veg=getValues(raster(x$Veg[j])))
        #    d <- filter(cells, Cell %in% which(!is.na(v$FID))) %>% mutate(Vegetation=factor(veg.labels[v$Veg[Cell]], levels=veg.labels), FID=v$FID[Cell]) %>%
        #        group_by(LocGroup, Location, Vegetation, FID) %>% summarise(Val=length(Cell), Var="Fire Size")
        #    d <- if(loopBy=="rep") mutate(d, Replicate=x$iter[j]) else if(loopBy=="year") mutate(d, Year=x$iter[j])
        #    d
        #}
        #system.time({ d.fs <- lapply(1:length(x$iter), fn, x, cells, veg.labels, loopBy) })
        
        for(j in 1:length(x$iter)){ # fire size by vegetation class
            v <- list(FID=getValues(raster(x$FID[j], band=2)), Veg=getValues(raster(x$Veg[j])))
            d <- filter(cells, Cell %in% which(!is.na(v$FID))) %>% mutate(Vegetation=factor(veg.labels[v$Veg[Cell]], levels=veg.labels), FID=v$FID[Cell]) %>%
                group_by(LocGroup, Location, Vegetation, FID) %>% summarise(Val=length(Cell), Var="Fire Size")
            d.fs[[j]] <- if(loopBy=="rep") mutate(d, Replicate=x$iter[j]) else if(loopBy=="year") mutate(d, Year=x$iter[j])
            print(switch(loopBy, "rep"=paste0("Year ", years[i], ": Replicate ", x$iter[j]), "year"=paste0("Replicate ", i, ": Year ", years[x$iter[j]])))
        }
    }
    d.fs <- if(loopBy=="rep") rbindlist(d.fs)[, Year:=as.integer(years[i])] else if(loopBy=="year") rbindlist(d.fs)[, Replicate:=as.integer(i)]
    d.fs <- setcolorder(d.fs, c("LocGroup", "Location", "Var", "Vegetation", "Year", "Val", "FID", "Replicate")) %>%
        group_by(LocGroup, Location, Var, Vegetation, Year) %>% setorder(Replicate, LocGroup, Location, Var, Vegetation, Year, Val)
    print(paste("Returning fire size by vegetation class data table."))
    d.fs
}

# @knitr extract_av
extract_av <- function(i, loopBy, mainDir, ageDir=NULL, reps=NULL, years=NULL, cells, readMethod="loop", ...){
    if(is.null(list(...)$veg.labels)) {
        veg.labels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")
    } else veg.labels <- list(...)$veg.labels
    if(!is.numeric(list(...)$n.samples)) n.samples <- 1000 else n.samples <- list(...)$n.samples
    x <- prep_alf_files(i=i, loopBy=loopBy, mainDir=mainDir, reps=reps, years=years)
    cells <- group_by(cells, LocGroup, Location)
    r <- getValues(raster(x$Age[1])) # use as a template
    idx <- which(!is.na(r))
    idx.rmNA <- which(idx %in% 1:length(r))
    d.age <- vector("list", length(x$iter))
    for(j in 1:length(x$iter)){
        v <- list(Age=getValues(raster(x$Age[j]))[idx], Veg=getValues(raster(x$Veg[j]))[idx])
        v$Age[v$Age < 0] <- v$Age[ v$Age < 0] + 2147483647 # temporary hack
        d <- filter(cells, Cell_rmNA %in% idx.rmNA) %>% mutate(Vegetation=factor(veg.labels[v$Veg[Cell_rmNA]], levels=veg.labels), Age=v$Age[Cell_rmNA]) %>%
            group_by(LocGroup, Location, Vegetation, Age) %>% summarise(Freq=length(Cell_rmNA))
        d.age[[j]] <- if(loopBy=="rep") mutate(d, Replicate=x$iter[j]) else if(loopBy=="year") mutate(d, Year=x$iter[j])
        print(switch(loopBy, "rep"=paste0("Year ", years[i], ": Replicate ", x$iter[j]), "year"=paste0("Replicate ", i, ": Year ", years[x$iter[j]])))
    }
    d.age <- if(loopBy=="rep") rbindlist(d.age)[, Year:=as.integer(years[i])] else if(loopBy=="year") rbindlist(d.age)[, Replicate:=as.integer(i)]
    d.age <- group_by(d.age, LocGroup, Location, Year, Vegetation) %>% setorder(Replicate, LocGroup, Location, Year, Vegetation, Age, Freq)
    d.area <- group_by(d.age, Replicate, add=T) %>% summarise(Val=sum(Freq))
    d.area[, Var:= "Vegetated Area"]
    setcolorder(d.area, c("LocGroup", "Location", "Var", "Vegetation", "Year", "Val", "Replicate"))
    locs <- unique(d.age$Location)
    if(loopBy=="rep"){
        d.age <- group_by(d.age, Age, add=T) %>% summarise(Freq=sum(Freq))
        d.age[, Var:= "Vegetation Age"]
        setcolorder(d.age, c("LocGroup", "Location", "Var", "Vegetation", "Year", "Age", "Freq"))
        return(list(d.area=d.area, d.age=d.age))
    }
    if(loopBy=="year"){
        setkey(d.age, Location)
        for(j in 1:length(locs)){
            obj.name.tmp <- paste0("age__", locs[j], "__rep", i)
            assign(obj.name.tmp, d.age[locs[j]])
            save(list=c("locs", obj.name.tmp), file=paste0(ageDir, "/", obj.name.tmp, ".RData"))
            print(paste(obj.name.tmp, "object", j, "of", length(locs), "saved."))
            rm(list=obj.name.tmp)
            gc()
        }
        rm(d.age)
        gc()
        return(list(d.area=d.area))
    }
}
