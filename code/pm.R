# @knitr create_project
source("C:/github/ProjectManagement/code/rpm.R") # eventually load a package instead of source script
proj.name <- "SNAPQAQC" # Project name
proj.location <- matt.proj.path # Use default file location

docDir <- c("Rmd/include", "md", "html", "Rnw", "pdf", "timeline")
newProject(proj.name, proj.location, docs.dirs=docDir, overwrite=T) # create a new project

rfile.path <- file.path(proj.location, proj.name, "code") # path to R scripts
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

proj.title <- "SNAP Data QAQC"
proj.menu <- c("Overview", "Code flow", "GCM climate", "CRU 3.x climate", "Fire & Vegetation", "Other R Code", "All Projects")

proj.submenu <- list(
	c("empty"),
	c("empty"),
	c("Extraction", "Downscaled GCMs", "divider", "Organization", "Regions: stats", "Regions: densities", "Cities: stats"),
	c("Extraction", "Downscaled CRU 3.x", "divider", "Organization", "Regions: stats", "Regions: densities", "Cities: stats"),
	c("Extraction", "ALFRESCO: Rmpi", "ALFRESCO: functions", "divider", "Organization", "ALFRESCO: stats and densities"),
	c("QAQC App Metadata"),
	c("empty")
)

proj.files <- list(
	c("index.html"),
	c("code_sankey.html"),
	c("header", "AR4_AR5_extract.html", "divider", "header", "stats_setup.html", "samples_setup.html", "cities_setup.html"),
	c("header", "CRU_extract.html", "divider", "header", "stats_setup_CRU.html", "samples_setup_CRU.html", "cities_setup_CRU.html"),
	c("header", "alfStatsByRep_Rmpi.html", "alfStatsByRep.html", "divider", "header", "getAlfStatsAndDensities.html"),
	c("qaqc_app_metadata.html"),
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