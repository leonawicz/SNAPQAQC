# @knitr setup
set.seed(47)
source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/functions.R")
baseDir <- "/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final"
topDir <- "alfresco"
dataset <- "veg"
reg.grp <- "LCC Regions"
reg <- "NW Interior Forest N"
samples <- TRUE
dtype <- if(samples) "samples" else "stats"
yrs.lim <- c(2010,2100)

veg.given <- "Black Spruce"
mod.given <- "CCCMAcgcm31"
scen.given <- "SRES A1B"
year.given <- 2050

if(!exists("agg.veg")) agg.veg <- FALSE # aggregate to forest and tundra
if(!exists("lb") | !exists("ub")) { lb <- 0.025; ub <- 0.975 } # confidence limits

dir.create(plotDir <- file.path("/workspace/UA/mfleonawicz/projects/SNAPQAQC/plots/AlfTest", dataset), recur=T, showWarnings=F) # improve pathing specificity
setwd(wDir <- file.path(baseDir, topDir, dtype, reg.grp, reg))
lapply(c("reshape2", "dplyr", "data.table", "ggplot2"), library, character.only=T)

system.time( load(switch(dataset, ba="baByVeg.RData", fc="fcByVeg.RData", fs="fsByVeg.RData", veg="vegetatedArea.RData", age="vegetationAge.RData")) )# about 30 seconds on Atlas CPU
dataname <- switch(dataset, ba="Burn area", fc="Fire frequency", fs="Fire size", veg="Vegetated area", age="Vegetation age")
dataname2 <- tolower(dataname)
d <- get(ls(pattern="^d\\."))
rm(list=ls(pattern="^d\\."))
gc()

# @knitr data_prep
if(agg.veg & dataset=="veg"){
    veg.levels <- c("Forest", "Tundra", "All")
    agg.veg.lab <- "FT"
} else {
    veg.levels <- c("Black Spruce", "White Spruce", "Deciduous", "Shrub Tundra", "Graminoid Tundra", "All")
    agg.veg.lab <- ""
}

# primary data table (can subset rows in accordance with conditional distributions of the RV)
# For random variable X:
# d can be used to investigate the following probability distributions directly:
# f(x(R),S,M), f(x(R),M|S=s), f(x(R),S|M=m), f(x(R)|S=s,M=m)
# conditioning implies subsetting the table and joint refers to what remains

# probability mass functions of factor variables ("ID column" variables, e.g., Scenario and Model),
# conditional on values of the random variable of interest, X (derived from "Val" and "Prob" columns),
# can also be examined, but this is not done here.

system.time(d %>% filter(Vegetation!="Wetland Tundra" & Year >= yrs.lim[1] & Year <= yrs.lim[2]) %>%
    mutate(Vegetation=factor(Vegetation, levels=veg.levels)) %>% group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% disttable -> d) # about five seconds on Atlas CPU

if(agg.veg & dataset=="veg") d <- toForestTundra(d)
if(dataset %in% c("veg", "ba")){
    unit.scale <- 1000
    d[, Val:=Val/unit.scale]
} else unit.scale <- ""

# @knitr derived_datasets
# secondary data tables (can summarize/merge sections of rows in accordance with marginal distributions of the RV)
system.time({
d.R <- marginalize(d, c("Scenario", "Model")) # f(x(R))
d.RgS <- marginalize(d, "Model") # f(x(R)|S=s)
d.RgM <- marginalize(d, "Scenario") # f(x(R)|M=m)
d.msdist <- msdisttable(d, d.RgS, d.RgM, d.R) # smart bind
}) # about 520 seconds on Atlas CPU

# tertiary data tables (summarizing uncertainty bounds for conditional and marginal distributions of the RV)
# conditional uncertainty
system.time({
d.uc <- uc_table(d)
d.RgS.uc <- uc_table(d.RgS)
d.RgM.uc <- uc_table(d.RgM)
d.uc.cond <- uc_table(d.uc, d.RgS.uc, d.RgM.uc)  # smart bind
# average marginal uncertainty
d.mean.uc <- uc_table(d, condition.on.mean=c("Scenario", "Model"))
d.RgM.mean.uc <- uc_table(d.RgM, condition.on.mean="Model")
d.RgS.mean.uc <- uc_table(d.RgS, condition.on.mean="Scenario")
d.R.mean.uc <- uc_table(d.R)
# compound average marginal uncertainty levels
d.uc.mar.compound <- uc_table(d.mean.uc, d.RgM.mean.uc, d.RgS.mean.uc, d.R.mean.uc) # smart bind
# stepwise individual marginal uncertainty components
d.uc.mar.component <- uc_components(d.uc.mar.compound)
}) # about 230 seconds on Atlas CPU

# Skip this for now
do.stepwise <- FALSE
if(do.stepwise){
library(combinat)
library(parallel)
mod.seq <- permn(unique(d$Model))
system.time( d.step <- mclapply(1:length(mod.seq), function(i, data, models) uc_stepwise_gcm(data=data, models=models[[i]]), data=d, models=mod.seq, mc.cores=32) ) # 1.75 to 2 hours on one Atlas node, about 180 GB peak RAM
d.step <- rbindlist(d.step)
d.step %>% select(-Phase, -Location, -Var) %>% filter(Type=="GCM", Vegetation=="Black Spruce", nchar(GCMset)==3) %>%
    mutate(GCMset=sapply(lapply(strsplit(GCMset, " "), sort), paste, collapse="")) %>%
    group_by(GCMset) %>% summarise(Magnitude=mean(Magnitude)) %>% setorder(-Magnitude) %>% print
save(d.step, file="/workspace/UA/mfleonawicz/projects/SNAPQAQC/workspaces/test_dt_stepwise_uncertainty.RData")
}

tables()

# @knitr plots
# plot_setup
setwd(plotDir)
w <- 4800
h <- 2800
res <- 300

if(dataset %in% c("veg", "ba", "fs")){
    lb_ts <- bquote(.(dataname)~(.(unit.scale)~km^2)~"")
    lb_tsu <- bquote(.(dataname)~"uncertainty"~(.(unit.scale)~km^2)~"")
} else if(dataset=="age") {
    lb_ts <- bquote(.(dataname)~(years)~"")
    lb_tsu <- bquote(.(dataname)~"uncertainty"~(years)~"")
} else if(dataset=="fc") {
    lb_ts <- dataname
    lb_tsu <- paste(dataname, "uncertainty")
}

# example histograms of probability distributions of a RV
png(paste0(dataset, "_hist", agg.veg.lab, ".png"), height=h, width=w, res=res)
distplot(d, facet.formula="Decade ~ Model", facet.scales="free_y", decade.start.years=seq(2010, 2090, by=20), Scenario=scen.given, Vegetation=veg.given, xlab=lb_ts)
dev.off()
png(paste0(dataset, "_hist2", agg.veg.lab, ".png"), height=h, width=w, res=res)
distplot(d, facet.formula="Model ~ Vegetation", decade.start.years=seq(2010, 2090, by=20), Scenario=scen.given, colour="Decade", xlab=lb_ts)
dev.off()

# example time series of conditional uncertainty and average marginal uncertainty for a RV with respect to GCMs and scenarios
# combined total
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucTotal.png"), height=h, width=w, res=res)
distplot(d.uc.mar.compound, type="total", facet.formula="~ Vegetation", facet.ncol=3, xlab="", ylab=lb_ts)
dev.off()
# compound
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucCompound.png"), height=h, width=w, res=res)
distplot(d.uc.mar.compound, type="compound", facet.formula="~ Vegetation", facet.ncol=3, xlab="", ylab=lb_tsu)
dev.off()
# individual components: total
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucComponentStack.png"), height=h, width=w, res=res)
distplot(d.uc.mar.component, type="stack", facet.formula="~ Vegetation", facet.ncol=3, xlab="")
dev.off()
# individual components: proportions
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucComponentProp.png"), height=h, width=w, res=res)
distplot(d.uc.mar.component, type="proportion", facet.formula="~ Vegetation", facet.ncol=3, xlab="")
dev.off()
# conditional uncertainty
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucConditional.png"), height=h, width=w, res=res)
distplot(d.uc.cond, type="conditional", facet.formula="~ Vegetation", facet.ncol=3, xlab="", ylab=lb_tsu)
dev.off()

# example cycling of density estimation and sampling of a RV
png(paste0(dataset, "_bootstrapDensityCycling_10k.png"), height=h, width=w, res=res)
distplot(d, dist.cycle=TRUE, Scenario=scen.given, Year=year.given, Vegetation=veg.given, Model=mod.given, xlab=lb_ts)
dev.off()
png(paste0(dataset, "_bootstrapDensityCycling2_10k.png"), height=h, width=w, res=res)
distplot(d, dist.cycle=TRUE, Scenario=scen.given, Year=year.given, facet.formula="Model ~ Vegetation", xlab=lb_ts)
dev.off()
png(paste0(dataset, "_bootstrapDensityCycling3_10k.png"), height=h, width=w, res=res)
distplot(d, dist.cycle=TRUE, Scenario=scen.given, Vegetation=veg.given, facet.formula="Model ~ Decade", decade.start.years=seq(2010, 2090, by=20), xlab=lb_ts)
dev.off()
png(paste0(dataset, "_bootstrapDensityCycling4_10k.png"), height=h, width=w, res=res)
distplot(d, dist.cycle=TRUE, Scenario=scen.given, group.vars=c("Model", "Location", "Var", "Vegetation"), facet.formula="Model ~ Vegetation", xlab=lb_ts)
dev.off()

# example comparisons of marginal and conditional distributions of a RV with respect to GCMs and scenarios
png(paste0(dataset, "_marCondDensity.png"), height=h, width=w, res=res)
distplot(d.msdist, facet.formula="~ DistType", Year=year.given, Vegetation=veg.given, xlab=lb_ts)
dev.off()
png(paste0(dataset, "_marCondDensity2.png"), height=h, width=w, res=res)
distplot(d.msdist, facet.formula="DistType ~ Vegetation", Year=year.given, xlab=lb_ts, strip.text.size=8)
dev.off()
png(paste0(dataset, "_marCondDensity3.png"), height=h, width=w, res=res)
distplot(d.msdist, facet.formula="DistType ~ Decade", Vegetation=veg.given, decade.start.years=seq(2010, 2090, by=20), xlab=lb_ts, , strip.text.size=8)
dev.off()
png(paste0(dataset, "_marCondDensity4.png"), height=h, width=w, res=res)
distplot(d.msdist, facet.formula="~ DistType", Vegetation=veg.given, xlab=lb_ts)
dev.off()
