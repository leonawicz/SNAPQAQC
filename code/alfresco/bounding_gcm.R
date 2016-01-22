# @knitr setup
set.seed(47)
wDir <- "/atlas_scratch/mfleonawicz/projects/SNAPQAQC"
setwd(wDir)
source("code/alfresco/functions.R")
projectName <- "IEM"
baseDir <- file.path("data/final/alfresco", projectName)
dataset <- c("ba")
reg.grp <- "Alaska L2 Ecoregions"
reg <- c("Intermontane Boreal")
samples <- TRUE
dtype <- if(samples) "samples" else "stats"
yrs.lim <- c(2008,2100)
period <- if(yrs.lim[1] >= 2008) "projected" else if(yrs.lim[2] <= 2007) "historical" else c("historical", "projected")
if(!exists("agg.veg")) agg.veg <- FALSE # aggregate to forest and tundra
if(!exists("lb") | !exists("ub")) { lb <- 0.025; ub <- 0.975 } # confidence limits

dir.create(plotDir <- file.path("plots", projectName, "bounding_models"), recur=T, showWarnings=F) # improve pathing specificity
lapply(c("dplyr", "data.table", "ggplot2"), library, character.only=T)
datname <- c()
d <- vector("list", length(reg))

for(i in 1:length(reg)){
    setwd(file.path(baseDir, dtype, reg.grp, reg[i]))
    dlist <- vector("list", length(dataset))
    for(j in 1:length(dataset)){
        infile <- switch(dataset[j], ba="baByVeg.RData", fc="fcByVeg.RData", fs="fsByVeg.RData", veg="vegarea.RData", age="vegage.RData")
        infile <- paste(period, infile, sep="_")
        d_tmp <- vector("list", length(infile))
        for(k in 1:length(infile)){
            system.time(load(infile[k])) # about 15 seconds on Atlas CPU
            if(i==1) dataname <- c(datname, switch(dataset[j], ba="Burn area", fc="Fire frequency", fs="Fire size", veg="Vegetated area", age="Vegetation age"))
            d_tmp[[k]] <- get(ls(pattern="^d\\.")) %>%
                filter(Year %in% yrs.lim[1]:yrs.lim[2] & (Vegetation %in% c("All"))) #%>% byDecade(decade.start.years=seq(2010,2090,by=10))
        }
        d_tmp <- rbindlist(d_tmp)
        dlist[[j]] <- d_tmp
        rm(list=ls(pattern="^d\\."))
        gc()
        print(paste("Region", i, "of", length(reg), ": Variable", j, "of", length(dataset), "loaded and prepped."))
    }
    d[[i]] <- rbindlist(dlist)
}
setwd(wDir)
rm(dlist, d_tmp)
d <- rbindlist(d)
gc()
dataname2 <- tolower(dataname)

# @knitr data_prep
d <- group_by(d, Phase, Scenario, Model, Location, Var, Vegetation) %>% disttable %>% marginalize(margin="Year")

# @knitr derived_datasets
# secondary data tables (can summarize/merge sections of rows in accordance with marginal distributions of the RV)
system.time({
d.RgM <- marginalize(d, "Scenario") # f(x(R)|M=m)
d.RgS <- marginalize(d, "Model") # f(x(R)|S=s)
}) # about 27 seconds on Atlas CPU

# inverse probability mass functions
system.time( d.pmf.M <- inverse_pmf(d, val.range=c(15000,1e6), var.new="Model") )
d.pmf.S <- inverse_pmf(d, val.range=c(1000,1e6), var.new="Scenario")
d.pmf.M.aggS <- inverse_pmf(d.RgM, val.range=c(1000,1e6), var.new="Model")
d.pmf.S.aggM <- inverse_pmf(d.RgS, val.range=c(1000,1e6), var.new="Scenario")

tables()

# @knitr plots
# plot_setup
setwd(plotDir)
w <- 3200
h <- 1600
res <- 300

# examples of inverse pmf
png("pmf_model_ts.png", height=h, width=w, res=res)
#distplot(d.pmf1, facet.formula="Location ~ Vegetation", fill="Model", xlab="")
ggplot(d.pmf.M, aes(Scenario, Prob, fill=Model)) + geom_bar(stat="identity", position="dodge") #+ facet_wrap(~Scenario, ncol=1)
dev.off()


