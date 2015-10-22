# @knitr setup
set.seed(47)
#source("C:/github/SNAPQAQC/code/alfresco/functions.R")
#baseDir <- "C:/github/SNAPQAQC/data"
source("/workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/functions.R")
baseDir <- "/workspace/UA/mfleonawicz/projects/SNAPQAQC/data/final"
topDir <- "alfresco"
dataset <- c("fs", "ba")
reg.grp <- "LCC Regions"
reg <- c("NW Interior Forest S", "NW Interior Forest N")
samples <- TRUE
dtype <- if(samples) "samples" else "stats"
yrs.lim <- c(2010,2100)

veg.given <- "Black Spruce"
mod.given <- "CCCMAcgcm31"
scen.given <- "SRES A1B"
dec.given <- "2050s"

if(!exists("agg.veg")) agg.veg <- FALSE # aggregate to forest and tundra
if(!exists("lb") | !exists("ub")) { lb <- 0.025; ub <- 0.975 } # confidence limits

dir.create(plotDir <- file.path("/workspace/UA/mfleonawicz/projects/SNAPQAQC/plots/AlfTest_multiVar"), recur=T, showWarnings=F) # improve pathing specificity
lapply(c("reshape2", "dplyr", "data.table", "ggplot2"), library, character.only=T)
datname <- c()
d <- vector("list", length(reg))

for(i in 1:length(reg)){
    setwd(wDir <- file.path(baseDir, topDir, dtype, reg.grp, reg[i]))
    dlist <- vector("list", length(dataset))
    for(j in 1:length(dataset)){
        system.time( load(switch(dataset[j], ba="baByVeg.RData", fc="fcByVeg.RData", fs="fsByVeg.RData", veg="vegarea.RData", age="vegage.RData")) )# about 15 seconds on Atlas CPU
        if(i==1) dataname <- c(datname, switch(dataset[j], ba="Burn area", fc="Fire frequency", fs="Fire size", veg="Vegetated area", age="Vegetation age"))
        dlist[[j]] <- get(ls(pattern="^d\\.")) %>%
            filter(Year %in% yrs.lim[1]:yrs.lim[2] & !(Vegetation %in% c("Barren lichen-moss", "Wetland Tundra"))) %>%
            byDecade(decade.start.years=seq(2010,2090,by=10))
        rm(list=ls(pattern="^d\\."))
        gc()
        print(paste("Region", i, "of", length(reg), ": Variable", j, "of", length(dataset), "loaded and prepped."))
    }
    d[[i]] <- rbindlist(dlist)
}

rm(dlist)
d <- rbindlist(d)
gc()
dataname2 <- tolower(dataname)

# @knitr data_prep

# primary data table (can subset rows in accordance with conditional distributions of the RV)
# For random variable X:
# d can be used to investigate the following probability distributions directly:
# f(x(R),S,M), f(x(R),M|S=s), f(x(R),S|M=m), f(x(R)|S=s,M=m)
# conditioning implies subsetting the table and joint refers to what remains

# probability mass functions of factor variables ("ID column" variables, e.g., Scenario and Model),
# conditional on values of the random variable of interest, X (derived from "Val" and "Prob" columns),
# can also be examined, but this is not done here.

d <- group_by(d, Phase, Scenario, Model, Location, Var, Vegetation, Decade) %>% disttable

#if(agg.veg & dataset=="veg") { d <- toForestTundra(d); agg.veg.lab <- "FT" } else agg.veg.lab <- ""
#if(dataset %in% c("veg", "ba")){
#    unit.scale <- 1000
#    d[, Val:=Val/unit.scale]
#} else 
unit.scale <- ""

# @knitr derived_datasets
# secondary data tables (can summarize/merge sections of rows in accordance with marginal distributions of the RV)
system.time({
d.R <- marginalize(d, c("Location", "Scenario", "Model")) # f(x(R))
d.RgM <- marginalize(d, c("Location", "Scenario")) # f(x(R))
d.RgLS <- marginalize(d, "Model") # f(x(R)|S=s)
d.RgLM <- marginalize(d, "Scenario") # f(x(R)|M=m)
d.RgSM <- marginalize(d, "Location") # f(x(R)|L=l)
#d.msdist <- msdisttable(d, d.RgLS, d.RgLM, d.RgSM, d.R) # smart bind
}) # about 260 seconds on Atlas CPU

# tertiary data tables (summarizing uncertainty bounds for conditional and marginal distributions of the RV)
# conditional uncertainty
drop.vars=c("Var", "Vegetation", "Decade", "Year")
system.time({
d.uc <- uc_table(d, drop.vars=drop.vars)
d.RgLS.uc <- uc_table(d.RgLS, drop.vars=drop.vars)
d.RgLM.uc <- uc_table(d.RgLM, drop.vars=drop.vars)
d.RgSM.uc <- uc_table(d.RgSM, drop.vars=drop.vars)
d.RgM.uc <- uc_table(d.RgM, drop.vars=drop.vars)
d.uc.cond <- uc_table(d.uc, d.RgLS.uc, d.RgLM.uc, d.RgSM.uc, d.RgM.uc, drop.vars=drop.vars)  # smart bind
# average marginal uncertainty
d.mean.uc <- uc_table(d, condition.on.mean=c("Location", "Scenario", "Model"), drop.vars=drop.vars)
d.RgLM.mean.uc <- uc_table(d.RgLM, condition.on.mean=c("Location", "Model"), drop.vars=drop.vars)
d.RgLS.mean.uc <- uc_table(d.RgLS, condition.on.mean=c("Location", "Scenario"), drop.vars=drop.vars)
d.RgSM.mean.uc <- uc_table(d.RgSM, condition.on.mean=c("Scenario", "Model"), drop.vars=drop.vars)
d.RgM.mean.uc <- uc_table(d.RgM, condition.on.mean="Model", drop.vars=drop.vars)
d.R.mean.uc <- uc_table(d.R, drop.vars=drop.vars)
# compound average marginal uncertainty levels
d.uc.mar.compound <- uc_table(d.mean.uc, d.RgLM.mean.uc, d.RgSM.mean.uc, d.RgM.mean.uc, d.R.mean.uc) # smart bind
d.uc.mar.compound2 <- uc_table(d.mean.uc, d.RgSM.mean.uc, d.RgLM.mean.uc, d.RgM.mean.uc, d.R.mean.uc) # smart bind
# stepwise individual marginal uncertainty components
d.uc.mar.component <- uc_components(d.uc.mar.compound)
d.uc.mar.component2 <- uc_components(d.uc.mar.compound2)
d.uc.mar.component2 <- bind_rows(
    mutate(d.uc.mar.component, Stepwise="Stepwise variable addition: Location, Scenario, Model"),
    mutate(d.uc.mar.component2, Stepwise="Stepwise variable addition: Scenario, Location, Model")
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
system.time( d.step <- mclapply(1:length(mod.seq), function(i, data, models) uc_stepwise(data=data, models=models[[i]], drop.vars=drop.vars), data=d, models=mod.seq, mc.cores=10) ) # about 320 seconds on Atlas CPU
d.step <- rbindlist(d.step)
class(d.step) <- unique(c("ucsteptable", "uccomptable", class(d.step)))
d.step %>% select(-Phase) %>% filter(Type=="Model") -> d.step.tmp
d.step.tmp <- d.step.tmp[grep("\n", d.step.tmp$GCMset),]
d.step2 <- group_by(d.step.tmp, Location, Var, Vegetation, Decade, GCMset) %>% summarise(Magnitude=mean(Magnitude)) %>% setorder(Location, Var, Vegetation, Decade, -Magnitude) #%>% mutate(GCMset=factor(GCMset, levels=rev(GCMset)))
class(d.step2) <- unique(c("ucsteptable", "uccomptable", class(d.step2)))
}

tables()

# @knitr plots
# plot_setup
setwd(plotDir)
w <- 4800
h <- 2800
res <- 300

#if(dataset %in% c("veg", "ba", "fs")){
#    lb_ts <- bquote(.(dataname)~(.(unit.scale)~km^2)~"")
#    lb_tsu <- bquote(.(dataname)~"uncertainty"~(.(unit.scale)~km^2)~"")
#    lb_comp_stack <- bquote(.(dataname)~"uncertainty"~(.(unit.scale)~km^2)~"")
#    lb_comp_prop <- paste(dataname, "proportional uncertainty")
#} else if(dataset=="age") {
#    lb_ts <- bquote(.(dataname)~(years)~"")
#    lb_tsu <- bquote(.(dataname)~"uncertainty"~(years)~"")
#    lb_comp_stack <- bquote(.(dataname)~"uncertainty"~(years)~"")
#    lb_comp_prop <- paste(dataname, "proportional uncertainty")
#} else if(dataset=="fc") {
#    lb_ts <- dataname
#    lb_tsu <- paste(dataname, "uncertainty")
#    lb_comp_stack <- paste(dataname, "uncertainty")
#    lb_comp_prop <- paste(dataname, "proportional uncertainty")
#}

# example of RV uncertainty by model pairs
png(paste0("ucModelPairs.png"), height=h, width=w, res=res)
distplot(d.step2 %>% filter(Var=="Burn Area" & Vegetation %in% c("Black Spruce", "Graminoid Tundra", "All")), facet.formula="Location ~ Vegetation", facet.scales="free_y", facet.ncol=3, ylab=expression("Burn area uncertainty"~(~km^2)~""), flip=FALSE)
dev.off()

# example time series of conditional uncertainty and average marginal uncertainty for a RV with respect to GCMs and scenarios
# combined total
png(paste0("tsByVeg", agg.veg.lab, "_ucTotal.png"), height=h, width=w, res=res)
distplot(d.uc.mar.compound %>% filter(Vegetation %in% c("Black Spruce", "Graminoid Tundra", "All")), type="total", facet.formula="Var ~ Vegetation", facet.ncol=3, xlab="", ylab=expression("Uncertainty"~(~km^2)~""))
dev.off()
# compound
png(paste0("tsByVeg", agg.veg.lab, "_ucCompound.png"), height=h, width=w, res=res)
distplot(d.uc.mar.compound %>% filter(Vegetation %in% c("Black Spruce", "Graminoid Tundra", "All")), type="compound", facet.formula="Var ~ Vegetation", facet.ncol=3, xlab="", ylab=expression("Uncertainty"~(~km^2)~""))
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
