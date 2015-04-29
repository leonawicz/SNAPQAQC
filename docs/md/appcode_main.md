


##
##
## SNAP data QA/QC Shiny app R code

The `main.R` script contains browser-side **R** code and is sourced by `ui.R`.
It employs conditional panels to break up the outputs associated with the different plot types available in the app.

### main.R


```r
column(8, help_tabpanel_conditional, conditionalPanel(condition = "input.tsp == 'plot_ts' && input.goButton !== null && input.goButton > 0", 
    plotOutput("PlotTS", width = "100%", height = "auto"), uiOutput("tsTextSub"), 
    uiOutput("TableTS")), conditionalPanel(condition = "input.tsp == 'plot_scatter' && input.goButton !== null && input.goButton > 0", 
    plotOutput("PlotScatter", width = "100%", height = "auto"), uiOutput("spTextSub"), 
    uiOutput("TableScatter")), conditionalPanel(condition = "input.tsp == 'plot_heatmap' && input.goButton !== null && input.goButton > 0", 
    plotOutput("PlotHeatmap", width = "100%", height = "auto"), uiOutput("hmTextSub"), 
    uiOutput("TableHeatmap")), conditionalPanel(condition = "input.tsp == 'plot_variability' && input.goButton !== null && input.goButton > 0", 
    plotOutput("PlotVariability", width = "100%", height = "auto"), uiOutput("varTextSub"), 
    uiOutput("TableVariability")), conditionalPanel(condition = "input.tsp == 'plot_spatial' && input.goButton !== null && input.goButton > 0", 
    plotOutput("PlotSpatial", width = "100%", height = "auto"), uiOutput("spatialTextSub"), 
    uiOutput("TableSpatial")))
```
