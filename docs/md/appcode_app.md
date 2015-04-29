


##
##
## SNAP data QA/QC Shiny app R code

The `app.R` script can be viewed as the central hub of the server-side **R** code.
It is sourced by `server.R`.
I keep it externalized because `server.R` also sources `serverHead.R`,
the initial server-side code external to and preceding the `shinyServer` call, which is also a sizeable script on its own.
In turn, `app.R` sources several other **R** scripts as shown below, covering the remainder of the server-side code that flows into `server.R`.
This includes all strict server-side reactive expressions, reactive outputs (assigned specifically to `output` via the `$` operator), and custom `ggplot2`-based plotting functions.

### app.R

#### Source R scripts


```r
# Source reactive expressions and other code
source("external/reactives.R", local = T)  # source reactive expressions
source("external/io_sidebar.R", local = T)  # source input/output objects associated with sidebar
source("external/io_main.R", local = T)  # source input/output objects associated with mainPanel

tsPlot <- source("external/plot_ts.R", local = T)$value
scatterPlot <- source("external/plot_scatter.R", local = T)$value
varPlot <- source("external/plot_variability.R", local = T)$value
heatPlot <- source("external/plot_heatmap.R", local = T)$value
spatialPlot <- source("external/plot_spatial.R", local = T)$value
```

#### Define plot wrapper functions

Each plot function sourced by `app.R` is wrapped in a parent function, which is similarly non-reactive but references known reactive expressions and reactive inputs.

##### Time series wrapper function


```r
# Specific plot function setup
doPlot_ts <- function(...) {
    if (permitPlot() & !is.null(input$group)) {
        if (!(input$group != "None" & !length(input$colorpalettes_ts))) {
            tsPlot(d = dat(), x = input$xtime, y = aggStatsID(), y.name = input$aggStats, 
                Log = input$log_ts, d.grp = datCollapseGroups(), d.pool = datCollapsePooled(), 
                grp = input$group, n.grp = n.groups(), ingroup.subjects = subjectChoices(), 
                panels = facet.panels(), facet.by = input$facet, vert.facet = input$vert.facet, 
                fontsize = input$plotFontSize, colpal = input$colorpalettes_ts, 
                linePlot = input$linePlot, barPlot = input$barPlot, pts.alpha = input$alpha1, 
                bartype = input$bartype, bardirection = input$bardirection, 
                show.points = input$ts_showpts, show.lines = input$ts_showlines, 
                show.overlay = input$ts_showCRU, overlay = CRU(), jit = input$ts_jitterXY, 
                plot.title = plot_ts_title(), plot.subtitle = plot_ts_subtitle(), 
                show.panel.text = input$ts_showPanelText, show.title = input$ts_showTitle, 
                lgd.pos = input$legendPos1, units = currentUnits(), yrange = input$yrange, 
                clbootbar = input$clbootbar, clbootsmooth = input$clbootsmooth, 
                pooled.var = pooled.var(), plot.theme.dark = input$ts_plotThemeDark, 
                logo.mat = logo.mat, ...)
        } else NULL
    } else NULL
}
```

##### Scatter plot wrapper function


```r
doPlot_scatter <- function(...) {
    if (permitPlot() & !is.null(input$group2)) {
        if (!(input$group2 != "None" & !length(input$colorpalettes_sc))) {
            scatterPlot(d = dat2(), x = input$vars, y = input$vars2, x.name = input$aggStats, 
                y.name = input$aggStats2, Logx = input$log_sc_x, Logy = input$log_sc_y, 
                flip.axes = sc_flip_xy(), grp = input$group2, n.grp = n.groups2(), 
                panels = facet.panels2(), facet.by = input$facet2, vert.facet = input$vert.facet2, 
                fontsize = input$plotFontSize2, colpal = input$colorpalettes_sc, 
                show.points = input$sc_showpts, contourlines = input$sc_showlines, 
                hexbin = input$hexbin, pts.alpha = input$alpha2, show.overlay = input$sc_showCRU, 
                overlay = CRU2(), jit = input$sc_jitterXY, plot.title = plot_sp_title(), 
                plot.subtitle = plot_sp_subtitle(), show.panel.text = input$sc_showPanelText, 
                show.title = input$sc_showTitle, lgd.pos = input$legendPos2, 
                units = currentUnits(), pooled.var = pooled.var2(), plot.theme.dark = input$sc_plotThemeDark, 
                logo.mat = logo.mat, ...)
        } else NULL
    } else NULL
}
```

##### Heat map wrapper function


```r
doPlot_heatmap <- function(...) {
    if (permitPlot() & !is.null(input$heatmap_x) & !is.null(input$heatmap_y) & 
        length(input$colorpalettes_hm)) {
        heatPlot(d = dat(), d.stat = aggStatsID(), d2 = dat_heatmap(), x = input$heatmap_x, 
            y = input$heatmap_y, z = input$statHeatmap, Log = input$log_hm, 
            panels = facetPanelsHeatmap(), facet.by = input$facetHeatmap, fontsize = input$plotFontSizeHeatmap, 
            colpal = input$colorpalettes_hm, reverse.colors = input$revHeatmapColors, 
            aspect_1to1 = input$aspect1to1, show.values = input$showHeatmapVals, 
            show.overlay = input$hm_showCRU, overlay = CRU(), plot.title = plot_hm_title(), 
            plot.subtitle = plot_hm_subtitle(), show.panel.text = input$hm_showPanelText, 
            show.title = input$hm_showTitle, lgd.pos = input$legendPosHeatmap, 
            units = currentUnits(), pooled.var = pooledVarHeatmap(), plot.theme.dark = input$hm_plotThemeDark, 
            logo.mat = logo.mat, ...)
    } else NULL
}
```

##### Variability plots wrapper function


```r
doPlot_var <- function(...) {
    if (permitPlot() & !is.null(pooled.var3()) & !is.null(input$group3)) {
        if (!(input$group3 != "None" & !length(input$colorpalettes_vr))) {
            varPlot(d = dat(), x = input$xvar, y = aggStatsID(), y.name = input$aggStats, 
                stat = stat(), around.mean = input$variability, d.grp = datCollapseGroups(), 
                d.pool = datCollapsePooled(), grp = input$group3, n.grp = n.groups3(), 
                ingroup.subjects = subjectChoices3(), panels = facet.panels3(), 
                facet.by = input$facet3, vert.facet = input$vert.facet3, fontsize = input$plotFontSize3, 
                colpal = input$colorpalettes_vr, boxplots = input$boxplots, 
                pts.alpha = input$alpha3, bartype = input$bartype3, bardirection = input$bardirection3, 
                show.points = input$vr_showpts, show.lines = input$vr_showlines, 
                show.overlay = input$vr_showCRU, overlay = CRU(), jit = input$vr_jitterXY, 
                plot.title = plot_var_title(), plot.subtitle = plot_var_subtitle(), 
                show.panel.text = input$vr_showPanelText, show.title = input$vr_showTitle, 
                lgd.pos = input$legendPos3, units = currentUnits(), yrange = input$yrange, 
                clbootbar = input$clbootbar, clbootsmooth = input$clbootsmooth, 
                plot.theme.dark = input$vr_plotThemeDark, logo.mat = logo.mat, 
                ...)
        } else NULL
    } else NULL
}
```

##### Distribution plots wrapper function


```r
doPlot_spatial <- function(...) {
    if (permitPlot() & !is.null(pooledVarSpatial()) & !is.null(input$groupSpatial)) {
        if (!(input$groupSpatial != "None" & !length(input$colorpalettes_sp))) {
            spatialPlot(d = dat_spatial(), x = input$spatial_x, y = "Val", grp = input$groupSpatial, 
                n.grp = nGroupsSpatial(), ingroup.subjects = subjectChoicesSpatial(), 
                plottype = input$plotTypeSpatial, thin.sample = as.numeric(input$thinSpatialSample), 
                panels = facetPanelsSpatial(), facet.by = input$facetSpatial, 
                vert.facet = input$vertFacetSpatial, fontsize = input$plotFontSizeSpatial, 
                colpal = input$colorpalettes_sp, linePlot = input$linePlotSpatial, 
                boxplots = input$boxplotsSpatial, pts.alpha = input$alphaSpatial, 
                density.type = input$densityTypeSpatial, strip.direction = input$stripDirectionSpatial, 
                show.points = input$sp_showpts, show.lines = input$sp_showlines, 
                show.overlay = input$sp_showCRU, overlay = CRU_spatial(), jit = input$sp_jitterXY, 
                plot.title = plot_spatial_title(), plot.subtitle = plot_spatial_subtitle(), 
                show.panel.text = input$sp_showPanelText, show.title = input$sp_showTitle, 
                lgd.pos = input$legendPosSpatial, units = currentUnits(), plot.theme.dark = input$sp_plotThemeDark, 
                logo.mat = logo.mat, ...)
        } else NULL
    } else NULL
}
```

Outputs below include browser-rendered plots, plot downloads, and data downloads.

#### Time series outputs


```r
# Primary outputs Time series plot
output$PlotTS <- renderPlot({
    if (twoBtnNullOrZero_ts()) 
        return()
    isolate({
        progress <- Progress$new(session, min = 1, max = 10)
        on.exit(progress$close())
        progress$set(message = "Generating plot...", value = 10)
        doPlot_ts(show.logo = F)
    })
}, height = function() {
    w <- if (twoBtnNullOrZero_ts()) 
        0 else session$clientData$output_PlotTS_width
    round((7/12) * w)
}, width = "auto")

output$dlCurPlotTS <- downloadHandler(filename = "timeseries.pdf", content = function(file) {
    if (input$ts_plotThemeDark) 
        bg <- "black" else bg <- "white"
    pdf(file = file, width = 1.5 * 12, height = 1.5 * 7, pointsize = 12, onefile = FALSE, 
        bg = bg)
    doPlot_ts(show.logo = T)
    dev.off()
})

output$dlCurTableTS <- downloadHandler(filename = function() {
    "timeseries_data.csv"
}, content = function(file) {
    write.csv(dat(), file)
})
```

#### Scatter plot outputs


```r
# Scatterplot
plot_scatter_ht <- function() {
    if (twoBtnNullOrZero_sc()) 
        return(0)
    ht <- 700
    if (!is.null(facet.panels2())) {
        cols <- ceiling(sqrt(facet.panels2()))
        rows <- ceiling(facet.panels2()/cols)
        ht <- ht * (rows/cols)
    }
    ht
}

output$PlotScatter <- renderPlot({
    if (twoBtnNullOrZero_sc()) 
        return()
    isolate({
        progress <- Progress$new(session, min = 1, max = 10)
        on.exit(progress$close())
        progress$set(message = "Generating plot...", value = 10)
        doPlot_scatter(show.logo = F)
    })
}, height = function() {
    w <- if (twoBtnNullOrZero_sc()) 
        0 else session$clientData$output_PlotScatter_width
    w
}, width = "auto")

output$dlCurPlotScatter <- downloadHandler(filename = "scatterplot.pdf", content = function(file) {
    if (input$sc_plotThemeDark) 
        bg <- "black" else bg <- "white"
    pdf(file = file, width = 1.5 * 12, height = 1.5 * 12, pointsize = 12, onefile = FALSE, 
        bg = bg)
    doPlot_scatter(show.logo = T)
    dev.off()
})

output$dlCurTableScatter <- downloadHandler(filename = function() {
    "scatterplot_data.csv"
}, content = function(file) {
    write.csv(dat2(), file)
})
```

#### Heat map outputs


```r
# Variability plot
output$PlotVariability <- renderPlot({
    if (twoBtnNullOrZero_vr()) 
        return()
    isolate({
        progress <- Progress$new(session, min = 1, max = 10)
        on.exit(progress$close())
        progress$set(message = "Generating plot...", value = 10)
        doPlot_var(show.logo = F)
    })
}, height = function() {
    w <- if (twoBtnNullOrZero_vr()) 
        0 else session$clientData$output_PlotVariability_width
    round((7/12) * w)
}, width = "auto")

output$dlCurPlotVariability <- downloadHandler(filename = "variability.pdf", 
    content = function(file) {
        if (input$vr_plotThemeDark) 
            bg <- "black" else bg <- "white"
        pdf(file = file, width = 1.5 * 12, height = 1.5 * 7, pointsize = 12, 
            onefile = FALSE, bg = bg)
        doPlot_var(show.logo = T)
        dev.off()
    })

output$dlCurTableVariability <- downloadHandler(filename = function() {
    "variability_data.csv"
}, content = function(file) {
    write.csv(dat(), file)
})
```

#### Variability plots outputs


```r
# Heatmap plot
output$PlotHeatmap <- renderPlot({
    if (twoBtnNullOrZero_hm()) 
        return()
    isolate({
        progress <- Progress$new(session, min = 1, max = 10)
        on.exit(progress$close())
        progress$set(message = "Generating plot...", value = 10)
        doPlot_heatmap(show.logo = F)
    })
}, height = function() {
    w <- if (twoBtnNullOrZero_hm()) 
        0 else session$clientData$output_PlotHeatmap_width
    round((7/12) * w)
}, width = "auto")

output$dlCurPlotHeatmap <- downloadHandler(filename = "heatmap.pdf", content = function(file) {
    if (input$hm_plotThemeDark) 
        bg <- "black" else bg <- "white"
    pdf(file = file, width = 1.5 * 12, height = 1.5 * 7, pointsize = 12, onefile = FALSE, 
        bg = bg)
    doPlot_heatmap(show.logo = T)
    dev.off()
})

output$dlCurTableHeatmap <- downloadHandler(filename = function() {
    "heatmap_data.csv"
}, content = function(file) {
    write.csv(dat(), file)
})
```

#### Distribution plots outputs


```r
# Spatial plot
output$PlotSpatial <- renderPlot({
    if (twoBtnNullOrZero_sp()) 
        return()
    isolate({
        progress <- Progress$new(session, min = 1, max = 10)
        on.exit(progress$close())
        progress$set(message = "Generating plot...", value = 10)
        doPlot_spatial(show.logo = F)
    })
}, height = function() {
    w <- if (twoBtnNullOrZero_sp()) 
        0 else session$clientData$output_PlotSpatial_width
    round((7/12) * w)
}, width = "auto")

output$dlCurPlotSpatial <- downloadHandler(filename = "spatial.pdf", content = function(file) {
    if (input$sp_plotThemeDark) 
        bg <- "black" else bg <- "white"
    pdf(file = file, width = 1.5 * 12, height = 1.5 * 7, pointsize = 12, onefile = FALSE, 
        bg = bg)
    doPlot_spatial(show.logo = T)
    dev.off()
})

output$dlCurTableSpatial <- downloadHandler(filename = function() {
    "spatial_data.csv"
}, content = function(file) {
    write.csv(dat_spatial(), file)
})
```
