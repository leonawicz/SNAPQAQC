# @knitr sankey_packages
require(igraph)
require(rCharts)

# @knitr files
out <- "X:/leonawicz/projects/AR4_AR5_comparisons/docs/diagrams/codeflow.html"

c0 <- c("AR4_AR5_extract.slurm", "CRU_extract.slurm")
c1 <- c("AR4_AR5_extract.R", "CRU_extract.R")
c2a <- c("cities_setup.R", "stats_setup.R", "samples_setup.R")
c2b <- c("cities_setup_CRU31.R", "stats_setup_CRU31.R", "samples_setup_CRU31.R")

d1 <- c("shapes2cells_AKCAN2km_5pct.RData", "Community location data: locs.csv", "2-km AK-CAN downscaled climate data")
d2 <- "Regional shapefiles"

c3 <- "shapes2cells.R"

e1 <- "qaqc_app_metadata.R"
e2 <- "Master QAQC app: meta.RData"
e3 <- "cities_meta_akcan2km.RData"

f1 <- c("alfStatsByRep_Rmpi.slurm", "alfStatsByRep_Rmpi.R", "alfStatsByRep.R", "Intermediary data", "getAlfStatsAndDensities.R")
f2 <- "1-km ALFRESCO outputs"
f3 <- c("shapes2cells_AKCAN1km.RData", "shapes2cells_AKCAN1km_rmNA.RData")

# @knitr links
from <- c(
	c0,
	rep(c1, times=c(length(c2a), length(c2b))),
	rep(d1, length(c1)),
	d2,
	rep(c3, 3),
	"cities_setup.slurm",
	c(c2a[2:3], c2b[3], e1),
	e1, e3,
	c(f1[1], f3),
	f1[2], f1[3], f1[4],
	rep(f1[5], 3),
	rep(f2, 2)
)

to <- c(
	c1,
	c2a, c2b,
	rep(c1, each=length(d1)),
	c3,
	c(d1[1], f3),
	c2a[1],
	rep(e2, 4),
	e3, e2,
	rep(f1[2], 3),
	f1[3], f1[4], f1[5],
	c(e2, "App data", "Other QAQC analyses"),
	f1[2:3]
)

val <- rep(1, length(to))

# @knitr igraph
relations <- data.frame(from=from, to=to)
g <- graph.data.frame(relations, directed=T, vertices=data.frame(c(c0, c1, c2a, c2b, c3, d1, d2)))

gw <- get.data.frame(g)
colnames(gw) <- c("source","target","value")
gw$value <- 1
gw$source <- as.character(gw$source)
gw$target <- as.character(gw$target)

# @knitr rcharts
p <- rCharts$new()
p$setLib('http://timelyportfolio.github.io/rCharts_d3_sankey/libraries/widgets/d3_sankey')
p$setTemplate(script = "http://timelyportfolio.github.io/rCharts_d3_sankey/libraries/widgets/d3_sankey/layouts/chart.html")
p$set(data=gw, nodeWidth=15, nodePadding=10, layout=32, width=900, height=800, margin=list(right=20, left=20, bottom=50, top=50), title="Code Flow")

p$setTemplate(
  afterScript="
<script>
  var cscale = d3.scale.category20b();
  d3.selectAll('#{{ chartId }} svg path.link')
    .style('stroke', function(d){
      return cscale(d.source.name);
    })
  d3.selectAll('#{{ chartId }} svg .node rect')
    .style('fill', function(d){
      return cscale(d.name)
    })
    .style('stroke', 'none')
</script>
")

# @knitr code_sankey_save
p$save(out)

# @knitr code_sankey_embed
p$show("iframesrc", cdn=T)
