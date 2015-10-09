


##
##
## SNAP data QA/QC Shiny app R code

`io_main.R` is a server-side **R** script sourced by `app.R`, which is in turn sourced by `server.R`.
It contains a number of reactive outputs rendered in the browser, specifically outputs which are sent to the main panel area of the app.

### io_main.R

#### Plot-specific sub-titling

These outputs are text which appear below a plot output and are not part of the plot itself.
They are used to draw attention to which variables selected by the user are pooled in a plot.

For example, the data informing a time series may contain multiple years, months, locations, and climate models.
If years are along the x-axis and the user has elected to color points or lines in the plot by locations and facet by month,
this still leaves climate models pooled together.
As a result, multiple points or lines of the same color can be displayed across the same x-axis values repeatedly.
These text displays below the plot simply help to remind the user what they have chosen to display.

Due to the complexity of the app and the myriad ways in which variables can be selected, combined, aggregated, or broken apart,
it is difficult to make the text outputs perfectly "smart", and not worth the effort compared to more important app developments.
However, they work decently enough, and it is because of their limitations that they appear below a plot rather than in the plot itself.


```r
output$tsTextSub <- renderUI({
    if (twoBtnNullOrZero_ts()) 
        return()
    isolate(pooledVarsCaption(pv = pooled.var(), permit = permitPlot(), ingrp = input$group))
})

output$spTextSub <- renderUI({
    if (twoBtnNullOrZero_sc()) 
        return()
    isolate(pooledVarsCaption(pv = pooled.var2(), permit = permitPlot(), ingrp = input$group2))
})

output$hmTextSub <- renderUI({
    if (twoBtnNullOrZero_hm()) 
        return()
    isolate(pooledVarsCaption(pv = pooledVarHeatmap(), permit = permitPlot()))
})

output$varTextSub <- renderUI({
    if (twoBtnNullOrZero_vr()) 
        return()
    isolate(pooledVarsCaption(pv = pooled.var3(), permit = permitPlot(), ingrp = input$group3))
})

output$spatialTextSub <- renderUI({
    if (twoBtnNullOrZero_sp()) 
        return()
    isolate(pooledVarsCaption(pv = pooledVarSpatial(), permit = permitPlot(), 
        ingrp = input$groupSpatial))
})
```

The data table outputs below are displayed in the app main panel area after a user has completed data selection, but prior to plotting.
Once a plot is displayed on the screen, the data table is bumped down to remain below the plot.

Note that the data table displayed to the screen for the variability plots  is the same as that shown for time series.
This is because they use the same source data, but conditional tab panel reactivity I have built into the app requires that they be defined by their own reactive outputs.

Also worth noting is that these tables are built from GCM data only.
The user option to overlay CRU data on GCM data in a plot is a plot option, not a data selection option.
As a result, tables are already compiled before CRU is considered.
Realistically, the tables are not worth literally reading in the browser, so it is not a big deal that CRU data, when selected under plot options. are not included in the tables.
Displaying the tables best serves the purpose of letting the user know something has not gone clearly wrong with their intended data selection and subsetting,
so that they may have confidence that their plot shows what they think it shows.
The data download option will mirror this absence of selected CRU data in the table.

#### Time series data table


```r
output$SubsetTableTS <- renderDataTable({
    if (!is.null(dat())) 
        dat()[, !"Decade", with = FALSE]
}, options = list(orderClasses = TRUE, lengthMenu = c(5, 10, 25, 50), pageLength = 5), 
    style = "bootstrap", rownames = F, filter = "bottom", caption = "Table 1: GCM data selection for time series plots.")

output$TableTS <- renderUI({
    if (goBtnNullOrZero()) 
        return()
    isolate(if (permitPlot()) 
        fluidRow(column(12, dataTableOutput("SubsetTableTS"))))
})
```

#### Scatter plot data table


```r
output$SubsetTableScatter <- renderDataTable({
    if (!is.null(dat2())) 
        dat2()[, !"Decade", with = FALSE]
}, options = list(orderClasses = TRUE, lengthMenu = c(5, 10, 25, 50), pageLength = 5), 
    style = "bootstrap", rownames = F, filter = "bottom", caption = "Table 2: GCM data selection for scatter plots.")

output$TableScatter <- renderUI({
    if (goBtnNullOrZero()) 
        return()
    isolate(if (permitPlot()) 
        fluidRow(column(12, dataTableOutput("SubsetTableScatter"))))
})
```

#### Heat map data table


```r
output$SubsetTableHeatmap <- renderDataTable({
    if (!is.null(dat_heatmap()) && nrow(dat_heatmap() > 0)) {
        if (ncol(dat_heatmap()) >= 9) 
            dat_heatmap()[, !"Decade", with = FALSE] else dat_heatmap()
    }
}, options = list(orderClasses = TRUE, lengthMenu = c(5, 10, 25, 50), pageLength = 5), 
    style = "bootstrap", rownames = F, filter = "bottom", caption = "Table 3: GCM data selection for heat maps.")

output$TableHeatmap <- renderUI({
    if (goBtnNullOrZero()) 
        return()
    isolate(if (permitPlot() & !is.null(dat_heatmap())) 
        if (nrow(dat_heatmap() > 0)) 
            fluidRow(column(12, dataTableOutput("SubsetTableHeatmap"))))
})
```

#### Variability plots data table


```r
output$SubsetTableVariability <- renderDataTable({
    if (!is.null(dat())) 
        dat()[, !"Decade", with = FALSE]
}, options = list(orderClasses = TRUE, lengthMenu = c(5, 10, 25, 50), pageLength = 5), 
    style = "bootstrap", rownames = F, filter = "bottom", caption = "Table 4: GCM data selection for variability assessment.")  # same as table 1

output$TableVariability <- renderUI({
    if (goBtnNullOrZero()) 
        return()
    isolate(if (permitPlot()) 
        fluidRow(column(12, dataTableOutput("SubsetTableVariability"))))
})
```

#### Distribution plots data table


```r
output$SubsetTableSpatial <- renderDataTable({
    if (!is.null(dat_spatial())) 
        dat_spatial()[, !"Decade", with = FALSE]
}, options = list(orderClasses = TRUE, lengthMenu = c(5, 10, 25, 50), pageLength = 5), 
    style = "bootstrap", rownames = F, filter = "bottom", caption = "Table 5: GCM data selection for spatial distributions.")

output$TableSpatial <- renderUI({
    if (goBtnNullOrZero()) 
        return()
    isolate(if (permitPlot()) 
        fluidRow(column(12, dataTableOutput("SubsetTableSpatial"))))
})
```
