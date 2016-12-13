##############################################################################
#### Stage two compilation of extracted Alfresco outputs into data tables ####
##############################################################################

# @knitr setup
comargs <- (commandArgs(TRUE))
if(!length(comargs)) q("no") else for(z in 1:length(comargs)) eval(parse(text=comargs[[z]]))

library(parallel)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)

if(!exists("period")) period <- stop("Must provide 'period' argument in escaped quotes. Options are 'historical' or 'projected'.")
stopifnot(length(period)==1 && period %in% c("historical", "projected"))
if(!exists("projectName")) projectName <- "Unnamed_Project_Run_Extractions"
if(!exists("mainDir")) mainDir <- "/atlas_scratch/mfleonawicz/alfresco"
mainDir <- file.path(mainDir, projectName, "extractions")
if(exists("gbm")) mainDir <- file.path(mainDir, gbm)
if(!exists("variable")) stop("Must provide 'variable' argument in escaped quotes. Options are 'age' (vegetation age), 'veg' (vegetated area), or 'fsv' (fire sizes by vegetation class).")
stopifnot(length(variable)==1 && variable %in% c("age", "veg", "fsv"))
if(!exists("reps")) reps <- NULL
inDir <- file.path(mainDir, variable)
if(exists("gbm")) projectName <- file.path(projectName, gbm)
dir.create(outDir <- file.path("/atlas_scratch/mfleonawicz/projects/SNAPQAQC/data/final/alfresco", projectName), showWarnings=F, recursive=T)
print(outDir)

# All regions for which stage-1 outputs currently exist on disk.
# Checking fsv files, but it is assumed stage one processing has been done on a common set of regions for all variables.
regions <- unique(sapply(strsplit(list.files(file.path(inDir)), "__"), "[", 2))
n.regions <- length(regions)
n.cores <- min(n.regions, 32)

# @knitr functions1
# Support functions
# Density estimation
dtDen <- function(x, n=1000, adj=0.1, out="vector", min.zero=TRUE, diversify=FALSE){
    x <- x[!is.na(x)]
    lx <- length(x)
    if(diversify && length(unique(x))==1) x <- rnorm(max(10, lx), mean=x[1]) # diversify constant values
    if(lx==1) x <- x + c(-1:1) #single pixel of veg type, add and subtract one age year to make procedure possible
    b <- max(1, 0.05*diff(range(x)))
    z <- density(x, adjust=adj, n=n, from=min(x)-b, to=max(x)+b)
    if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+b)
    if(out=="vector") return(as.numeric(c(z$x, z$y))) else if(out=="list") return(z)
}

# Bootstrapping
dtBoot <- function(p, p2=NULL, n.boot=10000, interp=TRUE, n.interp=100000, round.samples=FALSE){
    stopifnot(is.logical(round.samples) || is.na(as.integer(round.samples)))
    lp <- length(p)
	p <- if(is.null(p2)) list(x=p[1:(lp/2)], y=p[(lp/2+1):lp]) else list(x=p, y=p2)
	if(interp) p <- approx(p$x, p$y, n=n.interp)
    p <- sample(p$x, n.boot, prob=p$y, rep=T)
    if(round.samples==FALSE) return(p) else if(round.samples==TRUE) return(round(p)) else return(p, round.samples)
}

# @knitr prep_data
# Primary processing functions
prep_data <- function(j, inDir, outDir, n.samples=1000, n.boot=10000, period, ...){
  id <- basename(inDir)
  reps <- list(...)$reps
  exact <- list(...)$exact
  if(is.null(exact) || !is.logical(exact) || id!="veg") exact <- FALSE
	files.hist <- list.files(inDir, full=T, pattern=paste0("^", id, "__.*.CRU.*.RData$"))
  files.all <- list.files(inDir, full=T, pattern=paste0("^", id, "__.*.RData$"))
  if(period=="historical") files <- files.hist else if(period=="projected") files <- setdiff(files.all, files.hist)
	files.locs <- sapply(strsplit(files, "__"), "[", 2)
	locs <- unique(files.locs)
	if(j > length(locs)) return(NULL)
	loc <- locs[j]
	files <- files[which(files.locs %in% loc)]
  dat <- stat <- vector("list", length(files))
	for(i in 1:length(files)){
		load(files[i], envir=environment())
    d <- get(ls(pattern="^d\\."))
    rm(list=ls(pattern="^d\\."))
		loc.grp <- d$LocGroup[1]
		loc <- d$Location[1]
    if(id != "age"){ reps.all <- if(is.null(reps)) sort(unique(d$Replicate)) else reps }
    if(id=="fsv"){
      d2 <- group_by(d, Phase, Scenario, Model, Location, Var, Year, Replicate, FID) %>% summarise(Val=sum(Val)) %>% mutate(Vegetation="All") # agg-veg FS
      d <- data.table(bind_rows(d, d2)) %>% mutate(Vegetation=factor(Vegetation, levels=unique(Vegetation))) %>%
          group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) # individual and aggregate-veg fire sizes
      d2 <- group_by(d, Replicate, add=T) %>% summarise(BA=sum(Val), FC=length(Val)) # burn area and fire frequency
      d2 <- d2 %>% do(.,
          Expanded=suppressMessages(right_join(., data.table(Replicate=as.integer(reps.all)))) %>% # expand to include replicates with burn area and fire frequency of zero
              complete(nesting(Phase, Scenario, Model, Location, Var, Vegetation, Year), fill=list(BA=0L, FC=0L)) %>%
              fill(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% data.table
      ) %>% select(Expanded) %>% unnest(Expanded) %>% data.table %>% group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year)
      d <- dplyr::do(d, data.table::data.table(
                             Val=do.call(dtDen, list(.$Val, n=n.samples, out="list"))$x, Prob=do.call(dtDen, list(.$Val, n=n.samples, out="list"))$y))
      d2.ba <- dplyr::do(d2, data.table::data.table(
                             Val=do.call(dtDen, list(.$BA, n=n.samples, out="list"))$x, Prob=do.call(dtDen, list(.$BA, n=n.samples, out="list"))$y)) %>%
                             ungroup %>% mutate(Var="Burn Area")
      d2.fc <- dplyr::do(d2, data.table::data.table(
                             Val=do.call(dtDen, list(.$FC, n=n.samples, out="list"))$x, Prob=do.call(dtDen, list(.$FC, n=n.samples, out="list"))$y)) %>%
                             ungroup %>% mutate(Var="Fire Count")
      rm(d2)
      d <- data.table(bind_rows(d, d2.ba, d2.fc))
      rm(d2.ba, d2.fc)
      gc()
    }
		d <- group_by(d, Phase, Scenario, Model, Location, Var, Vegetation, Year)
    if(id=="age"){ # note that estimating age "distribution" is still allowed when all pixels for a given veg type have a contant age value across space in a given year
      d <- filter(d, !(Vegetation %in% c("Wetland Tundra", "Barren lichen-moss", "Temperate Rainforest")) & sum(Freq) > 30) # must be at least 30 pixels across all reps to estimate age distribution
      if(nrow(d)==0) return("Insufficient data.") # skip if any null data table for the location occurs for age, based on removal of irrelevant veg types or insufficient samples
      print(paste("j =", j, "| loc =", loc, "| loc.grp =", loc.grp, "| Age[1] = ", d$Age[1], "| Freq[1] =", d$Freq[1], "| nrow(d) =", nrow(d)))
      d <- mutate(d, N=n()) %>% group_by(Year, add=T)
      d <- dplyr::do(d, data.table::data.table( # constant age permitted (when length(Freq)==1)
                             Val=do.call(dtDen, list(ifelse(.$N==1, sample(.$Age, n.boot, T), sample(.$Age, n.boot, T, .$Freq)), n=n.samples, out="list"))$x,
                             Prob=do.call(dtDen, list(ifelse(.$N==1, sample(.$Age, n.boot, T), sample(.$Age, n.boot, T, .$Freq)), n=n.samples, out="list"))$y))
    }
    if(id=="veg"){
      dx <- d %>% do(.,
        Expanded=suppressMessages(right_join(., data.table(Replicate=as.integer(reps.all)))) %>% # expand to include replicates with veg area of zero
          complete(nesting(Phase, Scenario, Model, LocGroup, Location, Var, Vegetation, Year), fill=list(Val=0)) %>%
          fill(Phase, Scenario, Model, LocGroup, Location, Var, Vegetation, Year) %>% data.table
      ) %>% select(Expanded) %>% unnest(Expanded) %>% data.table %>% group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year)
      if(!exact){
        d <- dplyr::do(d, data.table::data.table(
                             Val=do.call(dtDen, list(.$Val, n=n.samples, out="list"))$x, Prob=do.call(dtDen, list(.$Val, n=n.samples, out="list"))$y))
      }
    }
    get_stats <- function(data, exact=FALSE){
      if(!exact) data <- dplyr::do(data, data.table::data.table(Val=do.call(dtBoot, list(.$Val,.$Prob, n.boot=n.boot))))
      summarise(data, Mean=round(mean(Val)), SD=round(sd(Val),1), Min=round(min(Val)),
        Pct_05=round(quantile(Val, 0.05)), Pct_10=round(quantile(Val, 0.10)), Pct_25=round(quantile(Val, 0.25)), Pct_50=round(quantile(Val, 0.50)),
        Pct_75=round(quantile(Val, 0.75)), Pct_90=round(quantile(Val, 0.90)), Pct_95=round(quantile(Val, 0.95)), Max=round(max(Val))) %>% group_by(Year, add=T)
    }
    s <- get_stats(d, exact=exact) # exact=TRUE only applies to veg area stats, for comparison with veg area stats computed after veg area density estimation
    if(id=="veg" & exact){
      d <- dplyr::do(d, data.table::data.table(
                             Val=do.call(dtDen, list(.$Val, n=n.samples, out="list"))$x, Prob=do.call(dtDen, list(.$Val, n=n.samples, out="list"))$y))
    }
    dat[[i]] <- d
    stat[[i]] <- s
    print(i)
	}
  dir.create(statsDir <- file.path(outDir, "stats", loc.grp, loc), recursive=T, showWarnings=F)
	dir.create(samplesDir <- file.path(outDir, "samples", loc.grp, loc), recursive=T, showWarnings=F)
  prefix <- if(stat[[1]]$Scenario[1]=="Historical") "historical" else "projected"
  if(id=="fsv"){
    dat <- rbindlist(dat)
    d.alf.fs <- filter(dat, Var=="Fire Size")
    d.alf.ba <- filter(dat, Var=="Burn Area")
    d.alf.fc <- filter(dat, Var=="Fire Count")
    stats.alf.fire <- rbindlist(stat)
    save(d.alf.fs, file=file.path(samplesDir, paste0(prefix, "_fsByVeg.RData")))
    save(d.alf.ba, file=file.path(samplesDir, paste0(prefix, "_baByVeg.RData")))
    save(d.alf.fc, file=file.path(samplesDir, paste0(prefix, "_fcByVeg.RData")))
    save(stats.alf.fire, file=file.path(statsDir, paste0(prefix, "_stats_fsbafcByVeg.RData")))
  } else {
    data.obj.name <- switch(id, age="d.alf.vegage", veg="d.alf.vegarea")
    stats.obj.name <- switch(id, age="stats.alf.vegage", veg="stats.alf.vegarea")
    assign(data.obj.name, rbindlist(dat))
    assign(stats.obj.name, rbindlist(stat))
    filename <- paste0(prefix, "_", tail(strsplit(data.obj.name, "\\.")[[1]], 1), ".RData")
    save(list=data.obj.name, file=file.path(samplesDir, filename))
    save(list=stats.obj.name, file=file.path(statsDir, filename))
  }
  return()
}

# @knitr run
mclapply(1:n.regions, prep_data, inDir, outDir, mc.cores=n.cores, period=period, reps=reps)
