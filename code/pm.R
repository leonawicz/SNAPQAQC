# @knitr create_project
source("C:/github/ProjectManagement/code/rpm.R") # eventually load a package instead of source script
proj.name <- "SNAPQAQC" # Project name
proj.location <- matt.proj.path # Use default file location

docDir <- c("Rmd/include", "md", "html", "Rnw", "pdf", "timeline")
newProject(proj.name, proj.location, docs.dirs=docDir, overwrite=T) # create a new project

appcode.files <- list.files("C:/github/shiny-apps/ar4ar5", pattern="\\.R$", full=TRUE, recursive=TRUE)
rfile.path <- file.path(proj.location, proj.name, "code") # path to R scripts
file.copy(appcode.files, paste0(rfile.path, "/climate/ar4ar5_appcode/appcode_", basename(appcode.files)), overwrite=TRUE)
docs.path <- file.path(proj.location, proj.name, "docs")
rmd.path <- file.path(docs.path, "Rmd")

# generate Rmd files from existing R scripts using default yaml front-matter
genRmd(path=rfile.path) # specify header.args list argument if necessary

# @knitr update_project
# update yaml front-matter only
genRmd(path=rfile.path, update.header=TRUE)

# obtain knitr code chunk names in existing R scripts
chunkNames(path=file.path(proj.location, proj.name, "code"))

# append new knitr code chunk names found in existing R scripts to any Rmd files which are outdated
chunkNames(path=file.path(proj.location, proj.name, "code"), append.new=TRUE)

# @knitr website
# Setup for generating a project website
user <- "leonawicz"
proj.github <- file.path("https://github.com", user, proj.name)
index.url <- "index.html"
#file.copy(index.url, "index.html")

proj.title <- "SNAP Data QA/QC"
proj.menu <- c("Overview", "Climate", "Fire & Veg", "Analyses", "App docs", "App code", "All Projects")

proj.submenu <- list(
	c("empty"),
	c("Extraction", "Downscaled GCMs", "Downscaled CRU 3.x", "2-km PRISM", "divider", "GCM organization", "Regions: stats", "Regions: densities", "Cities: stats", "CRU organization", "Regions: stats", "Regions: densities", "Cities: stats"),
	c("Extraction", "ALFRESCO: Rmpi", "ALFRESCO: functions", "divider", "Organization", "ALFRESCO: stats and densities"),
	c("Autocorrelation", "Intro/example"),
	c("About", "divider", "Help documents",
		"Getting started", "Working with data", "Graphical options", "Updating settings",
			"Graphing: time series", "Graphing: scatter plots", "Graphing: heat maps", "Graphing: variability", "Graphing: distributions",
		"divider", "Prep R code", c("QA/QC App Metadata")),
	c("App R code",
			c("global.R", "ui.R", "server.R", "about.R"), "divider",
			c("sidebar.R", "main.R", "serverHead.R", "app.R"), "divider",
			c("io_sidebar.R", "io_main.R", "reactives.R"), "divider",
			c("plot_ts.R", "plot_scatter.R", "plot_heatmap.R", "plot_variability.R", "plot_spatial.R")),
	c("empty")
)

proj.files <- list(
	c("index.html"),
	c("header", "AR4_AR5_extract.html", "CRU_extract.html", "PRISM_extract.html", "divider", "header", "stats_setup.html", "samples_setup.html", "cities_setup.html", "header", "stats_setup_CRU.html", "samples_setup_CRU.html", "cities_setup_CRU.html"),
	c("header", "alfStatsByRep_Rmpi.html", "alfStatsByRep.html", "divider", "header", "getAlfStatsAndDensities.html"),
	c("header", "sp_ac_clump.html"),
	c("aboutapp.html", "divider", "header",
		paste0("help0", c("1_start", "2_data", "3_plotOptions", "4_updating", "5_01_graphTS", "5_02_graphScatter", "5_03_graphHeat", "5_04_graphVar", "5_05_graphSpatial"), ".html"),
		"divider", "header", c("qaqc_app_metadata.html")),
	c("header",
			c("appcode_global.html", "appcode_ui.html", "appcode_server.html", "appcode_about.html"), "divider",
			c("appcode_sidebar.html", "appcode_main.html", "appcode_serverHead.html", "appcode_app.html"), "divider",
			c("appcode_io_sidebar.html", "appcode_io_main.html", "appcode_reactives.html"), "divider",
			c("appcode_plot_ts.html", "appcode_plot_scatter.html", "appcode_plot_heatmap.html", "appcode_plot_variability.html", "appcode_plot_spatial.html")),
	c("http://leonawicz.github.io")
)
			
# generate navigation bar html file common to all pages
genNavbar(htmlfile=file.path(proj.location, proj.name, "docs/Rmd/include/navbar.html"), title=proj.title, menu=proj.menu, submenus=proj.submenu, files=proj.files, title.url="index.html", home.url="index.html", site.url=proj.github, include.home=FALSE)

# generate _output.yaml file
# Note that external libraries are expected, stored in the "libs" below
yaml.out <- file.path(proj.location, proj.name, "docs/Rmd/_output.yaml")
libs <- "libs"
common.header <- "include/in_header.html"
genOutyaml(file=yaml.out, lib=libs, header=common.header, before_body="include/navbar.html")

# @knitr knit_setup
library(rmarkdown)
library(knitr)
setwd(rmd.path)

# Rmd files
files.Rmd <- list.files(pattern=".Rmd$", full=T)

# @knitr save
# write all yaml front-matter-specified outputs to Rmd directory for all Rmd files
lapply(files.Rmd, render, output_format="all")
insert_gatc(list.files(pattern=".html$"))
moveDocs(path.docs=docs.path)

# if also making PDFs for a project, speed up the Rmd to Rnw file conversion/duplication
rnw.path <- file.path(docs.path, "Rnw")
setwd(rnw.path)
#themes <- knit_theme$get()
highlight <- "solarized-dark"
convertDocs(path=rmd.path, emphasis="replace", overwrite=TRUE, highlight=highlight) # Be careful
lapply(list.files(pattern=".Rnw$"), knit2pdf)
moveDocs(path.docs=docs.path, type="pdf", remove.latex=FALSE)
