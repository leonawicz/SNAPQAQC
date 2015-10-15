# @knitr setup
set.seed(47)
#source("C:/github/SNAPQAQC/code/alfresco/functions.R")
#baseDir <- "C:/github/SNAPQAQC/data"
source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/functions.R")
baseDir <- "/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final"
topDir <- "alfresco"
dataset <- "fs"
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

#dir.create(plotDir <- file.path("C:/github/SNAPQAQC/plots/AlfTest", dataset), recur=T, showWarnings=F) # improve pathing specificity
#setwd(wDir <- baseDir)
dir.create(plotDir <- file.path("/workspace/UA/mfleonawicz/projects/SNAPQAQC/plots/AlfTest", dataset), recur=T, showWarnings=F) # improve pathing specificity
setwd(wDir <- file.path(baseDir, topDir, dtype, reg.grp, reg))
lapply(c("reshape2", "dplyr", "data.table", "ggplot2"), library, character.only=T)

system.time( load(switch(dataset, ba="baByVeg.RData", fc="fcByVeg.RData", fs="fsByVeg.RData", veg="vegarea.RData", age="vegage.RData")) )# about 15 seconds on Atlas CPU
dataname <- switch(dataset, ba="Burn area", fc="Fire frequency", fs="Fire size", veg="Vegetated area", age="Vegetation age")
dataname2 <- tolower(dataname)
d <- get(ls(pattern="^d\\."))
rm(list=ls(pattern="^d\\."))
gc()

# @knitr data_prep

# primary data table (can subset rows in accordance with conditional distributions of the RV)
# For random variable X:
# d can be used to investigate the following probability distributions directly:
# f(x(R),S,M), f(x(R),M|S=s), f(x(R),S|M=m), f(x(R)|S=s,M=m)
# conditioning implies subsetting the table and joint refers to what remains

# probability mass functions of factor variables ("ID column" variables, e.g., Scenario and Model),
# conditional on values of the random variable of interest, X (derived from "Val" and "Prob" columns),
# can also be examined, but this is not done here.

d <- filter(d, !(Vegetation %in% c("Barren lichen-moss", "Wetland Tundra", "All")) & Year >= yrs.lim[1] & Year <= yrs.lim[2]) %>% group_by(Phase, Scenario, Model, Location, Var, Vegetation, Year) %>% disttable

if(agg.veg & dataset=="veg") { d <- toForestTundra(d); agg.veg.lab <- "FT" } else agg.veg.lab <- ""
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
}) # about 260 seconds on Atlas CPU

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
d.uc.mar.compound2 <- uc_table(d.mean.uc, d.RgS.mean.uc, d.RgM.mean.uc, d.R.mean.uc) # smart bind
# stepwise individual marginal uncertainty components
d.uc.mar.component <- uc_components(d.uc.mar.compound)
d.uc.mar.component2 <- uc_components(d.uc.mar.compound2)
d.uc.mar.component2 <- bind_rows(
    mutate(d.uc.mar.component, Stepwise="Stepwise variable addition: Scenario, then Model"),
    mutate(d.uc.mar.component2, Stepwise="Stepwise variable addition: Model, then Scenario")
) %>% data.table
class(d.uc.mar.component2) <- c("uc_comp_table", class(d.uc.mar.component2))
}) # about 300 seconds on Atlas CPU

# stepwise uncertainty increase with addition of models
do.stepwise <- TRUE
if(do.stepwise){
library(parallel)
library(purrr)
id <- expand.grid(1:5, 1:5) # limit stepwise comparison to all pairwise combinations of models
id <- id[id[,1]!=id[,2],]
mod.seq <- unique(lapply(strsplit(paste(id[,1], id[,2]), " "), function(x,m) m[as.numeric(sort(x))], m=unique(d$Model)))
system.time( d.step <- mclapply(1:length(mod.seq), function(i, data, models) uc_stepwise(data=data, models=models[[i]]), data=d, models=mod.seq, mc.cores=10) ) # about 320 seconds on Atlas CPU
d.step <- rbindlist(d.step)
class(d.step) <- unique(c("ucsteptable", "uccomptable", class(d.step)))
d.step %>% select(-Phase, -Location, -Var) %>% filter(Type=="Model" & Vegetation==veg.given) -> d.step.tmp
d.step.tmp <- d.step.tmp[grep("\n", d.step.tmp$GCMset),]
d.step2 <- group_by(d.step.tmp, GCMset) %>% summarise(Magnitude=mean(Magnitude)) %>% setorder(-Magnitude) %>% mutate(GCMset=factor(GCMset, levels=rev(GCMset)))
class(d.step2) <- unique(c("ucsteptable", "uccomptable", class(d.step2)))
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
    lb_comp_stack <- bquote(.(dataname)~"uncertainty"~(.(unit.scale)~km^2)~"")
    lb_comp_prop <- paste(dataname, "proportional uncertainty")
} else if(dataset=="age") {
    lb_ts <- bquote(.(dataname)~(years)~"")
    lb_tsu <- bquote(.(dataname)~"uncertainty"~(years)~"")
    lb_comp_stack <- bquote(.(dataname)~"uncertainty"~(years)~"")
    lb_comp_prop <- paste(dataname, "proportional uncertainty")
} else if(dataset=="fc") {
    lb_ts <- dataname
    lb_tsu <- paste(dataname, "uncertainty")
    lb_comp_stack <- paste(dataname, "uncertainty")
    lb_comp_prop <- paste(dataname, "proportional uncertainty")
}

# example of RV uncertainty by model pairs
png(paste0(dataset, "_ucModelPairs", agg.veg.lab, ".png"), height=h, width=w, res=res)
distplot(d.step2, facet.formula=NULL, facet.scales="free_y", Scenario=scen.given, Vegetation=veg.given, ylab=lb_tsu)
dev.off()

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
distplot(d.uc.mar.component, type="stack", facet.formula="~ Vegetation", facet.ncol=3, xlab="", ylab=lb_comp_stack)
dev.off()
# individual components: proportions
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucComponentProp.png"), height=h, width=w, res=res)
distplot(d.uc.mar.component, type="proportion", facet.formula="~ Vegetation", facet.ncol=3, xlab="", ylab=lb_comp_prop)
dev.off()
# stepwise variable addition comparison, individual components: total
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucComponentStack_stepwise.png"), height=h, width=w, res=res)
distplot(d.uc.mar.component2 %>% filter(Vegetation=="Black Spruce"), type="stack", facet.formula="~ Stepwise", facet.ncol=2, xlab="", ylab=lb_comp_stack)
dev.off()
# stepwise variable addition comparison, individual components: proportions
png(paste0(dataset, "_tsByVeg", agg.veg.lab, "_ucComponentProp_stepwise.png"), height=h, width=w, res=res)
distplot(d.uc.mar.component2 %>% filter(Vegetation=="Black Spruce"), type="proportion", facet.formula="~ Stepwise", facet.ncol=2, xlab="", ylab=lb_comp_prop)
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
