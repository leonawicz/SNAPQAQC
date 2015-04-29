


##
##
## SNAP data QA/QC Shiny app R code

`server.R` is very short because it sources longer **R** scripts, `serverHead.R` and `app.R`, the latter of which sources additional **R** code.

### server.R


```r
source("external/serverHead.R", local = TRUE)
shinyServer(function(input, output, session) source("external/app.R", local = TRUE))
```
