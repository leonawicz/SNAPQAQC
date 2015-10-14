# @knitr functions
# Distribution data table class constructor
disttable <- function(data, check.values=TRUE){
    stopifnot(any(class(data) %in% c("grouped_dt", "tbl_dt", "tbl", "data.table", "data.frame")))
    stopifnot(all(c("Val", "Prob") %in% names(data)))
    if(check.values){
        stopifnot(is.numeric(data$Val) && is.numeric(data$Prob))
        stopifnot(!any(is.na(data$Val)) && !any(is.na(data$Prob)))
        stopifnot(min(data$Prob) >= 0)
    }
    class(data) <- unique(c("disttable", class(data)))
    data
}

# bootstrap sample from estimated pdf
dtBoot <- function(p, p2=NULL, n.boot=10000, interp=TRUE, n.interp=100000, round.samples=FALSE){
    stopifnot(is.logical(round.samples) || is.na(as.integer(round.samples)))
    lp <- length(p)
	p <- if(is.null(p2)) list(x=p[1:(lp/2)], y=p[(lp/2+1):lp]) else list(x=p, y=p2)
	if(interp) p <- approx(p$x, p$y, n=n.interp)
    p <- sample(p$x, n.boot, prob=p$y, rep=T)
    if(round.samples==FALSE) return(p) else if(round.samples==TRUE) return(round(p)) else return(p, round.samples)
}

# sample from probability densities data table and replace with bootstrap samples data table
sample_densities <- function(data, check.class=TRUE){
    if(check.class && class(data)[1]!="disttable") class(data) <- unique(c("disttable", class(data)))
    if("Prob" %in% names(data)) return(summarise(data, Val=dtBoot(Val, Prob)))
    summarise(data, Val=dtBoot(Val))
}

# summarize distributions of a RV with specified uncertainty bounds around mean, conditioned on or marginalized over other variables
# constructs a uc_table class object from specifically disttable class objects or row binds multiple uc_table objects
uc_table <- function(..., lb=0.025, ub=0.975, condition.on.mean=NULL, margin=NULL, bindType.as.factor=TRUE, ignore.constants=TRUE, var.set=c("Phase", "Scenario", "Model", "Location", "Var", "Vegetation", "Year"), use.var.set=TRUE, check.class=TRUE){
    dots <- list(...)
    if(!length(dots)) stop("No data provided.")
    if(length(dots) > 1){
        stopifnot(all(sapply(dots, function(x) class(x)[1]=="uc_table")))
        stopifnot(sum(diff(sapply(dots, function(x) attributes(x)$lb)))==0)
        stopifnot(sum(diff(sapply(dots, function(x) attributes(x)$ub)))==0)
        data <- rbindlist(dots, fill=T)
        if(bindType.as.factor) data <- mutate(data, Type=factor(Type, levels=unique(Type)))
        cl <- class(data)
        cl <- unique(c("uc_table", cl[cl!="disttable"]))
        class(data) <- cl
        attr(data, "lb") <- lb
        attr(data, "ub") <- ub
        return(data)
    }
    stopifnot(is.null(condition.on.mean) || is.character(condition.on.mean))
    if(any(margin %in% c("Val", "Prob"))) stop("Invalid marginalization.")
    if(any(var.set %in% c("Val", "Prob"))) stop("Invalid variable set.")
    data <- dots[[1]]
    if(check.class && class(data)[1]!="disttable") data <- disttable(data)
    id <- names(data)
    stopifnot(all(margin %in% id))
    
    get_conditioning_vars <- function(x, vars, ignore.constants){
        if(!(ignore.constants && length(vars))) return(vars)
        idx <- c()
        for(i in 1:length(vars)) if(nrow(unique(x[, vars[i], with=F]))==1) idx <- c(idx, i)
        if(length(idx)) vars <- vars[-idx]
        if(!length(vars)) vars <- ""
        vars
    }
    
    get_marginalized_vars <- function(x, vars, var.set, use.var.set){
        if(!use.var.set) return(c("Sim", vars))
        stopifnot(length(var.set) > 0)
        idx <- which(!(var.set %in% names(x)))
        if(length(idx)) vars <- c(var.set[idx], vars)
        return(c("Sim", vars))
    }
    
    conditional.vars <- get_conditioning_vars(x=data, vars=id[!(id %in% c("Val", "Prob", margin))], ignore.constants=ignore.constants)
    stopifnot(all(condition.on.mean %in% conditional.vars))
    marginalized.vars <- get_marginalized_vars(x=data, vars=margin, var.set=var.set, use.var.set=use.var.set)
    if(length(margin)) data <- marginalize(data, margin=margin)
    id <- names(data)
    label.mar <- paste0(marginalized.vars, collapse=" + ")
    label.con <- if(conditional.vars[1]=="") "" else paste("|", paste0(conditional.vars, collapse=", "))
    label <- paste(label.mar, label.con) 
    dots <- lapply(id[!(id %in% c("Val", "Prob"))], as.symbol)
    sample_densities(data) %>% group_by_(.dots=dots) %>% summarise(LB=quantile(Val, lb), Mean=mean(Val), UB=quantile(Val, ub)) %>%
        mutate(Magnitude=UB-LB, Type=label) %>% group_by_(.dots=dots) -> data
    if(!is.null(condition.on.mean)){
        id <- names(data)
        dots2 <- lapply(id[!(id %in% c("LB", "Mean", "UB", "Magnitude", condition.on.mean))], as.symbol)
        group_by_(data, .dots=dots2) %>% summarise(LB=mean(LB), Mean=mean(Mean), UB=mean(UB), Magnitude=mean(Magnitude)) %>% group_by_(.dots=dots2) -> data
        data <- setcolorder(data, id[id %in% names(data)])
    }
    cl <- class(data) 
    cl <- unique(c("uc_table", cl[cl!="disttable"]))
    class(data) <- cl
    attr(data, "lb") <- lb
    attr(data, "ub") <- ub
    data
}

# estimated pdf from bootstrap sample
dtDen <- function(x, n=1000, adj=0.1, out="vector", min.zero=TRUE, diversify=FALSE){
    b <- max(1, 0.05*diff(range(x)))
    z <- density(x, adjust=adj, n=n, from=min(x)-b, to=max(x)+b)
    if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+b)
    if(out=="vector") return(as.numeric(c(z$x, z$y))) else if(out=="list") return(z)
}

# reclass Vegetation column to aggregate Forest and Tundra entries
toForestTundra <- function(data, group.vars=as.character(groups(data)), keep.prob.col=TRUE, check.class=TRUE){
    if(check.class && class(data)[1]!="disttable") data <- disttable(data)
    nam <- names(data)
    tundra <- c("Graminoid Tundra", "Shrub Tundra", "Wetland Tundra")
    f1 <- function(x) summarise(x, Val=dtBoot(Val))
    f2 <- function(x) summarise(x, Val=dtDen(Val))
    if("Prob" %in% nam){
        f1 <- function(x) summarise(x, Val=dtBoot(Val, Prob))
        if(keep.prob.col) f2 <- function(x) summarise(x, Val=dtDen(Val, out="list")$x, Prob=dtDen(Val, out="list")$y)
    }
    f1(data) %>% mutate(Vegetation2=ifelse(as.character(Vegetation) %in% tundra, "Tundra", "Forest"), Obs=1:nboot) %>%
        select(-Vegetation) %>% mutate(Vegetation=Vegetation2) %>% select(-Vegetation2) %>%
        group_by_(.dots=lapply(c(group.vars, "Obs"), as.symbol)) %>% summarise(Val=sum(Val)) %>%
        f2 %>% group_by_(.dots=lapply(group.vars, as.symbol)) %>% setcolorder(nam) %>% setkey -> data
    class(data) <- unique(c("disttable", class(data)))
    data
}

# merge distributions using a cycle of bootstrap resampling followed by density re-estimation
merge_densities <- function(data, keep.prob.col=TRUE, check.class=TRUE){
    if(check.class && class(data)[1]!="disttable") data <- disttable(data)
    g <- as.character(groups(data))
    if("Prob" %in% names(data)){
        if(keep.prob.col){
            data <- summarise(data, Val=dtDen(dtBoot(Val, Prob), out="list")$x, Prob=dtDen(dtBoot(Val, Prob), out="list")$y)
        } else data <- summarise(data, Val=dtDen(dtBoot(Val, Prob)))
    } else data <- summarise(data, Val=dtDen(dtBoot(Val)))
    if(all(g %in% names(data))) data <- group_by_(data, .dots=lapply(g, as.symbol))
    class(data) <- unique(c("disttable", class(data)))
    data
}

# Repeat cycle of bootstrap resampling followed by density re-estimation n-1 times, assumes Prob column present
bootDenCycle <- function(data, n, start=NULL, group.vars=as.character(groups(data)), check.class=TRUE){
    if(check.class && class(data)[1]!="disttable") data <- disttable(data)
    if(!("Cycle" %in% names(data))) data$Cycle <- 1
    if(is.null(start)) start <- max(data$Cycle)
    stopifnot(is.null(group.vars) || all(group.vars %in% names(data)))
    if(n<=1){
        class(data) <- unique(c("disttable", class(data)))
        return(select_(data, .dots=lapply(c(group.vars, "Val", "Prob", "Cycle"), as.symbol)) %>% group_by(Cycle, add=T))
    }
    data %>% bind_rows( filter(data, Cycle==start) %>% merge_densities %>% mutate(Cycle=start+1) ) %>%
        data.table %>% group_by_(.dots=lapply(group.vars, as.symbol)) %>% disttable %>% bootDenCycle(n-1, start+1, group.vars)
}

# marginalize distribution of RV (Val column) over selected categorical variables (Scenario, Model)
marginalize <- function(data, margin, keep.prob.col=TRUE, check.class=TRUE){
    if(check.class && class(data)[1]!="disttable") data <- disttable(data)
    id <- names(data)
    if(length(margin) && any(!(margin %in% id))) stop("Marginalizing variable(s) not found.")
    dots <- lapply(id[!(id %in% c("Val", "Prob", margin))], as.symbol)
    group_by_(data, .dots=dots) %>% merge_densities(keep.prob.col, check.class) %>% group_by_(.dots=dots) -> data
    class(data) <- unique(c("disttable", class(data)))
    data
}

# make data table of uncertainty components contributed to a total uncertainty by various factors
# constructs a uc_com_table class object from specifically uc_table class objects
# stepwise individual uncertainty components, given uncertainty from variables already marginalized over
uc_components <- function(data){
    stopifnot(class(data)[1]=="uc_table")
    require(reshape2)
    data <- mutate(data, Type2=sapply(strsplit(data[, as.character(Type)], " \\| "), "[", 1), LB=NULL, Mean=NULL, UB=NULL)
    x.len <- sapply(strsplit(data[, Type2], " \\+ "), length)
    stepwise.vars <- data[, Type2][match(unique(x.len), x.len)]
    variable.order <- c()
    for(i in 1:length(stepwise.vars)){
        v <- strsplit(stepwise.vars[i], " \\+ ")[[1]]
        variable.order[i] <- v[!(v %in% variable.order)][1]
    }
    stopifnot(all(stepwise.vars %in% data[, Type2]))
    data <- filter(data, Type2 %in% stepwise.vars) %>% mutate(Type2=factor(Type2, levels=stepwise.vars))
    string <- paste(paste(names(data)[-which(names(data) %in% c("Magnitude", "Type", "Type2"))], collapse=" + "), "~ Type2")
    data <- data.table(dcast(data, as.formula(string), value.var="Magnitude"))
    for(i in 2:length(stepwise.vars)){
        expr <- lazyeval::interp(~x-y, x=as.name(stepwise.vars[i]), y=as.name(stepwise.vars[i-1]))
        data <- mutate_(data, .dots=setNames(list(expr), variable.order[i]))
    }
    for(i in 2:length(stepwise.vars)){
        expr <- lazyeval::interp(~list(NULL), x=as.name(stepwise.vars[i]))
        data <- mutate_(data, .dots=setNames(list(expr), stepwise.vars[i]))
    }
    data <- data.table(melt(data, measure.vars=variable.order, variable.name="Type", value.name="Magnitude")) %>% 
        mutate(Type=factor(Type, levels=variable.order))
    data[Magnitude < 0, Magnitude:=0]
    class(data) <- unique(c("uc_comp_table", class(data)))
    data
}

# make stepwise GCM data table of uncertainty components (simulation, scenario, and model)
uc_stepwise_gcm <- function(data, models=NULL){
    uc.mar <- c("Sim", "Sim + Scenario", "Sim + Model", "Sim + Scenario + Model")
    uc.mar <- factor(uc.mar, levels=uc.mar)
    if(is.null(models)) models <- unique(data$Model)
    for(i in 1:length(models)){
        gcmset <- models[1:i]
        abb <- sapply(gcmset, function(x) switch(x, "CCCMAcgcm31"="C", "GFDLcm21"="G", "MIROC32m"="M", "MPIecham5"="E", "ukmoHADcm3"="H", "unknown"))
        data %>% filter(Model %in% gcmset) %>%
        (function(d) {
            d %>% uncertainty(drop=c("Scenario", "Model"), label=uc.mar[1]) -> d.sim
            d %>% marginalize("Scenario") %>% uncertainty(drop="Model", label=uc.mar[2]) -> d.simScen
            d %>% marginalize("Model") %>% uncertainty(drop="Scenario", label=uc.mar[3]) -> d.simMod
            d %>% marginalize(c("Scenario", "Model")) %>% uncertainty(label=uc.mar[4]) -> d.simScenMod
            uc_components(sim=d.sim, simScen=d.simScen, simMod=d.simMod, simScenMod=d.simScenMod)
        }) %>% mutate(GCMset=paste(abb, collapse=" ")) -> d.tmp
        if(i==1) d <- d.tmp else d <- bind_rows(d, d.tmp)
    }
    class(d) <- unique(c("ucsteptable", "uccomptable", "disttable", class(d)))
    d
}

# make stepwise GCM data table of uncertainty components (simulation, scenario, and model)
uc_stepwise <- function(data, models=NULL, use.abb=FALSE){
    if(is.null(models)) models <- unique(data$Model)
    stopifnot(length(models) > 1)
    for(i in 1:length(models)){
        gcmset <- models[1:i]
        abb <- sapply(gcmset, function(x) switch(x, "CCCMAcgcm31"="C", "GFDLcm21"="G", "MIROC32m"="M", "MPIecham5"="E", "ukmoHADcm3"="H", "unknown"))
        abb.collapse <- " "
        if(!use.abb) {abb <- gcmset; abb.collapse <- "\n"}
        data %>% filter(Model %in% gcmset) %>%
        (function(d) {
            d %>% uc_table(condition.on.mean=c("Scenario", "Model"), ignore.constants=FALSE) -> d.sim
            print(1)
            d %>% uc_table(condition.on.mean="Model", margin="Scenario", ignore.constants=FALSE) -> d.simScen
            print(2)
            d %>% uc_table(condition.on.mean="Scenario", margin="Model", ignore.constants=FALSE) -> d.simMod
            print(3)
            d %>% uc_table(margin=c("Scenario", "Model"), ignore.constants=FALSE) -> d.simScenMod
            print(4)
            uc_components(uc_table(d.sim, d.simScen, d.simMod, d.simScenMod))
        }) %>% mutate(GCMset=paste(abb, collapse=abb.collapse)) -> d.tmp
        if(i==1) d <- d.tmp else d <- bind_rows(d, d.tmp)
    }
    d <- data.table(d)
    class(d) <- unique(c("ucsteptable", "uccomptable", class(d)))
    d
}

# add Decade column based on Year column and optionally filter decades using vector of decade leading years
byDecade <- function(data, decade.start.years=NULL) {
    g <- as.character(groups(data))
    data <- mutate(data, Decade=paste0(10*Year%/%10, "s"))
    if(!is.null(decade.start.years)) data <- filter(data, Decade %in% paste0(decade.start.years, "s"))
    data %>% group_by_(.dots=lapply(c(g, "Decade"), as.symbol)) %>% marginalize("Year")
}

# compare marginal and conditional distributions of a RV with respect to models and scenarios
# constructs a msdisttable class object from specifically disttable class objects
msdisttable <- function(data, data.RgS, data.RgM, data.R){
    stopifnot(all(sapply(list(data, data.RgS, data.RgM, data.R), function(x) class(x)[1]=="disttable")))
    x <- c("Model", "Scenario")
    stopifnot(all(x %in% names(data)) && !(x[1] %in% names(data.RgS)) && x[2] %in% names(data.RgS) &&
        x[1] %in% names(data.RgM) && !(x[2] %in% names(data.RgM)) && !any(x %in% data.R))
    g <- groups(data)
    dist.types <- c("given scenario and model", "given scenario", "given model", "marginalized over both")
    data <- mutate(data, DistType=factor(dist.types[1], levels=dist.types))
    data.RgS <- mutate(data.RgS, DistType=factor(dist.types[2], levels=dist.types))
    data.RgM <- mutate(data.RgM, DistType=factor(dist.types[3], levels=dist.types))
    data.R <- mutate(data.R, DistType=factor(dist.types[4], levels=dist.types))
    bind_rows(data, data.RgS, data.RgM, data.R) %>% data.table %>%
        group_by_(.dots=lapply(c(as.character(g), "DistType"), as.symbol)) -> data
    class(data) <- unique(c("msdisttable", "disttable", class(data)))
    data
}

distplot <- function(x, ...) UseMethod("distplot")

# helper function for filtering factor variable levels in multiple plotting functions
.filter_plot_data <- function(data, dots){
    if(!length(dots)) return(data)
    dots <- dots[names(dots) %in% names(data)]
    if(!length(dots)) return(data)
    dots <- lapply(1:length(dots), function(i, dots){
            if(!is.numeric(dots[[i]][1])) dots[[i]] <- paste0("\'", dots[[i]], "\'")
            paste0(names(dots)[i], " %in% c(", paste(dots[[i]], collapse=","), ")")
        }, dots=dots)
    data <- filter_(data, .dots=dots)
    data
}

# helper function aggregating annual distributions of a RV to decades by marginalizing over years by decade
.aggToDecades <- function(data, group.vars=as.character(groups(data)), decades=NULL){
    if(is.null(decades)) return(data)
    data <- byDecade(data, decades)
    group.vars[group.vars=="Year"] <- "Decade"
    group_by_(data, .dots=lapply(group.vars, as.symbol))
}

# helper function for making a default plot title when title is not specified as a ... argument
.plottitle_label <- function(d, id.vars=c("Scenario", "Model", "Year", "Vegetation", "Var", "Location"), suffix=""){
    f <- function(x, d) if(!(x %in% names(d)) || nrow(unique(d[,x, with=F])) > 1) return("") else return(as.character(unlist(unique(d[,x, with=F]))))
    paste(gsub(" +", " ", paste(sapply(id.vars, f, d=d), collapse=" ")), suffix)
}

# Plot distributions of a RV after optionally conditioning (filtering) and/or marginalizing over (merging) factor levels of associated variables
# Plot repeated cycle of bootstrap resampling followed by density re-estimation of a RV with dist.cycle=TRUE and n > 1
distplot.disttable <- function(data, n=10, dist.cycle=FALSE, group.vars=as.character(groups(data)), facet.formula=NULL, facet.scales="free", show.plot=TRUE, return.data=TRUE, ...){
    stopifnot(is.null(group.vars) || all(group.vars %in% names(data)))
    merge.factors <- if(length(group.vars) < length(groups(data))) TRUE else FALSE
    dots <- list(...)
    .filter_plot_data(data, dots) %>% group_by_(.dots=lapply(group.vars, as.symbol)) %>%
        .aggToDecades(group.vars, dots$decade.start.years) -> data
    if(merge.factors) data <- merge_densities(data)
    if(dist.cycle) data <- bootDenCycle(data, n)
    data <-sample_densities(data)
    if(show.plot){
        colour <- if(!is.null(dots$colour)) dots$colour else NULL
        clrs <- if(!is.null(dots$color.vec)) dots$color.vec else c("#E69F00", "#0072B2", "#CC79A7", "#D55E00", "#009E73")
        xlb <- if(!is.null(dots$xlab)) dots$xlab else "X"
        ylb <- if(!is.null(dots$ylab)) dots$ylab else "Density"
        suffix <- if(dist.cycle) "pdf cycling" else "distribution"
        title <- if(!is.null(dots$title)) dots$title else .plottitle_label(data, suffix=suffix)
        if(dist.cycle){
            g <- ggplot(data=data %>% filter(Cycle==1), aes(x=Val)) +
                geom_histogram(aes(y=..density..), colour="black", fill="gray") +
                geom_line(data=data, aes(colour=factor(Cycle)), stat="density")
        } else {
            g <- ggplot(data=data, aes_string(x="Val", colour=colour)) +
                geom_histogram(aes(y=..density..), colour="black", fill="gray") + geom_line(stat="density", fill="black")
        }
        g <- g + guides(colour=guide_legend(override.aes=list(alpha=1))) + labs(x=xlb, y=ylb, title=title) + theme_bw(base_size=16) +
            theme(legend.position="bottom", legend.box="horizontal", axis.text.y=element_blank(), axis.ticks.y=element_blank())
        if(!is.null(colour) && nrow(unique(data[, colour, with=F])) <= length(clrs)) g <- g + scale_colour_manual(name="", values=clrs)
        if(!is.null(facet.formula)) g <- g + facet_wrap(as.formula(facet.formula), scales=facet.scales)
        print(g)
    }
    if(return.data) return(data)
}

# Plot comparison of marginal and conditional distributions of a RV with respect to GCMs and scenarios
distplot.msdisttable <- function(data, group.vars=as.character(groups(data)), facet.formula=NULL, facet.scales="free", show.plot=TRUE, return.data=TRUE, ...){
    dots <- list(...)
    group_by_(data, .dots=lapply(group.vars, as.symbol)) %>% .filter_plot_data(dots) %>%
        .aggToDecades(decades=dots$decade.start.years) %>% sample_densities -> data
    if(show.plot){
        xlb <- if(!is.null(dots$xlab)) dots$xlab else "X"
        ylb <- if(!is.null(dots$ylab)) dots$ylab else "Density"
        title <- if(!is.null(dots$title)) dots$title else .plottitle_label(data, c("Year", "Vegetation", "Var", "Location"), "marginal and conditional distributions")
        clrs <- if(!is.null(dots$color.vec)) dots$color.vec else c("#E69F00", "#0072B2", "#CC79A7", "#D55E00", "#009E73")
        strip.text.size <- if(!is.null(dots$strip.text.size)) dots$strip.text.size else 16
        g <- ggplot(data=data %>% filter(DistType=="given scenario and model"), aes(x=Val, colour=Model, linetype=Scenario, group=interaction(Scenario, Model))) +
            scale_colour_manual(name="GCM", values=clrs) +
            geom_line(stat="density") +
            geom_line(data=data %>% filter(DistType=="given model"), aes(linetype=NULL, group=Model), stat="density") +
            geom_line(data=data %>% filter(DistType=="given scenario"), aes(colour=NULL, linetype=Scenario, group=Scenario), stat="density") +
            geom_histogram(data=data %>% filter(DistType=="marginalized over both"),
                aes(y=..density.., fill=NULL, colour=NULL, linetype=NULL, group=NULL), colour="black", fill="gray", show_guide=F) +
            geom_line(data=data %>% filter(DistType=="marginalized over both"), aes(colour=NULL, linetype=NULL, group=NULL), stat="density") +
            theme_bw(base_size=16) + 
            theme(legend.position="bottom", legend.box="horizontal", axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                strip.text=element_text(size=strip.text.size)) +
            guides(colour=guide_legend(override.aes=list(alpha=1))) + labs(x=xlb, y=ylb, title=title)
        if(!is.null(facet.formula)) g <- g + facet_wrap(as.formula(facet.formula), scales=facet.scales)
        print(g)
    }
    if(return.data) return(data)
}

# Plot total or compound (by component) average marginal uncertainty in a RV over time with respect to underlying simulation uncertainty and other stepwise added variables
# optionally condition (filter) on factor levels of associated with other variables not marginalized over
distplot.uc_table <- function(data, type="total", facet.formula=NULL, facet.scales="free_y", facet.ncol=1, show.plot=TRUE, return.data=TRUE, ...){
    dots <- list(...)
    data <- .filter_plot_data(data, dots)
    x <- as.character(unique(data$Type))
    idx <- which.max(sapply(strsplit(x, " \\+ "), length))
    label.total <- x[idx]
    stopifnot(length(idx)==1)
    if(show.plot){
        if(!(type %in% c("total", "compound", "conditional"))) stop("type must be 'total', 'compound', or 'conditional'.")
        yrs.brks <- seq(min(data$Year) - min(data$Year) %% 10, max(data$Year) + (10 - max(data$Year) %% 10), by=10)
        clrs2 <- c("#00000030", "#00000030")
        names(clrs2) <- c("Uncertainty", paste("Combined uncertainty:", label.total))
        colour <- if(!is.null(dots$colour)) dots$colour else NULL
        clrs <- if(!is.null(dots$color.vec)) dots$color.vec else c("#E69F00", "#0072B2", "#CC79A7", "#D55E00", "#009E73")
        xlb <- if(!is.null(dots$xlab)) dots$xlab else "X"
        ylb <- if(!is.null(dots$ylab)) dots$ylab else "Y"
        prefix <- if(type=="conditional") "Conditional uncertainty in projected annual" else ""
        if(type=="total") suffix <- "with uncertainty" else if(type=="compound") suffix <- "compound uncertainty" else suffix <- ""
        title <- if(!is.null(dots$title)) dots$title else paste0(prefix, .plottitle_label(data, id.vars=c("Year", "Vegetation", "Var", "Location"), suffix=suffix))
        if(type=="total"){
            g <- ggplot(data=data %>% filter(Type==label.total), aes_string(x="Year", y="Mean", colour=colour), environment=environment()) +
                geom_ribbon(aes(ymin=LB, ymax=UB, fill=paste("Combined uncertainty:", label.total))) +
                geom_line(size=1) + geom_point(size=2)
        } else if(type=="compound"){
            g <- ggplot(data=data, aes(x=Year, y=Magnitude, colour=Type)) + geom_line(size=1) + expand_limits(y=0)
        } else if(type=="conditional"){
            g <- ggplot(data=data, aes_string(x="Year", y="Magnitude", colour="Type")) + geom_point() + expand_limits(y=0)
        }
        g <- g + scale_colour_manual(name="", values=clrs) + scale_fill_manual(name="",values=clrs2) + scale_x_continuous(breaks=yrs.brks) +
            theme_bw(base_size=16) + theme(legend.position="bottom", legend.box="horizontal") +
            guides(colour=guide_legend(override.aes=list(alpha=1))) + labs(x=xlb, y=ylb, title=title)
        if(!is.null(facet.formula)) g <- g + facet_wrap(as.formula(facet.formula), scales=facet.scales, ncol=facet.ncol)
        print(g)
    }
    if(return.data) return(data)
}

# Plot stacks or proprtions of average simulation, scenario, and model component marginal uncertainty in a RV over time
# optionally condition (filter) on other (not model or scenario) factor levels of associated variables
distplot.uc_comp_table <- function(data, type="stack", facet.formula=NULL, facet.scales="free_y", facet.ncol=1, show.plot=TRUE, return.data=TRUE, ...){
    dots <- list(...)
    data <- .filter_plot_data(data, dots)
    if(show.plot){
        if(!(type %in% c("stack", "proportion"))) stop("type must be 'stack' or 'proportion'.")
        yrs.brks <- seq(min(data$Year) - min(data$Year) %% 10, max(data$Year) + (10 - max(data$Year) %% 10), by=10)
        clrs <- if(!is.null(dots$color.vec)) dots$color.vec else c("#E69F00", "#0072B2", "#CC79A7", "#D55E00", "#009E73")
        xlb <- if(!is.null(dots$xlab)) dots$xlab else "X"
        if(!is.null(dots$ylab)) ylb <- dots$ylab else if(type=="stack") ylb <- "Combined uncertainty by source" else ylb <- "Proportion of combined uncertainty"
        prefix <- if(type=="stack") "Combined uncertainty in projected annual " else "Proportional uncertainty in projected annual "
        suffix <- "\nEstimated uncertainty component by source"
        title <- if(!is.null(dots$title)) dots$title else paste0(prefix, .plottitle_label(data, id.vars=c("Year", "Vegetation", "Var", "Location"), suffix=suffix))
        g <- ggplot(data=data, aes(x=Year, y=Magnitude, fill=Type))
        g <- if(type=="stack") g + geom_bar(stat="identity") else if(type=="proportion") g + geom_bar(stat="identity", position="fill")
        g <- g + scale_fill_manual(name="",values=clrs) + scale_x_continuous(breaks=yrs.brks) +
            theme_bw(base_size=16) + theme(legend.position="bottom", legend.box="horizontal") +
            guides(colour=guide_legend(override.aes=list(alpha=1))) + labs(x=xlb, y=ylb, title=title)
        if(!is.null(facet.formula)) g <- g + facet_wrap(as.formula(facet.formula), scales=facet.scales, ncol=facet.ncol)
        print(g)
    }
    if(return.data) return(data)
}

# Plot stacks or proprtions of average simulation, scenario, and model component marginal uncertainty in a RV over time
# optionally condition (filter) on other (not model or scenario) factor levels of associated variables
distplot.ucsteptable <- function(data, facet.formula=NULL, facet.scales="free_y", facet.ncol=1, show.plot=TRUE, return.data=TRUE, ...){
    dots <- list(...)
    data <- .filter_plot_data(data, dots)
    if(show.plot){
        clrs <- if(!is.null(dots$color.vec)) dots$color.vec else c("#E69F00", "#0072B2", "#CC79A7", "#D55E00", "#009E73")
        xlb <- if(!is.null(dots$xlab)) dots$xlab else "Model pairs"
        if(!is.null(dots$ylab)) ylb <- dots$ylab else ylb <- "Model pair uncertainty"
        prefix <- "Average uncertainty by model pair "
        title <- if(!is.null(dots$title)) dots$title else paste0(prefix, .plottitle_label(data, id.vars=c("Year", "Vegetation", "Var", "Location")))
        g <- ggplot(data=data, aes(x=GCMset, y=Magnitude))
        g <- g + geom_bar(stat="identity") + coord_flip()
        g <- g + scale_fill_manual(name="",values=clrs) +
            theme_bw(base_size=16) + theme(legend.position="bottom", legend.box="horizontal") +
            guides(colour=guide_legend(override.aes=list(alpha=1))) + labs(x=xlb, y=ylb, title=title)
        if(!is.null(facet.formula)) g <- g + facet_wrap(as.formula(facet.formula), scales=facet.scales, ncol=facet.ncol)
        print(g)
    }
    if(return.data) return(data)
}

