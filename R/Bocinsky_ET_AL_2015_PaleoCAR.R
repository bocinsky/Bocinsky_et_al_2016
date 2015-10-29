### Script for PRISM climate retrodiction using PaleoCAR reported in
### R. Kyle Bocinsky, Johnathan Rush, Keith W. Kintigh, and Timothy A. Kohler, 
### Macrohistory of the Prehispanic Pueblo Southwest: Cycles of Exploration and Exploitation

### WARNING: This will take a long time to run on large areas.

# Set the working directory to the location of R-scripts
setwd("/Volumes/BOCINSKY_DATA/WORKING/BOCINSKY_ET_AL_2015/R/")

# Load the functions for all analyses below
# install.packages("devtools", dependencies=T, repos = "http://cran.rstudio.com")
# update.packages(ask=F, repos = "http://cran.rstudio.com")
devtools::install_github("bocinsky/FedData@2.0.0")
devtools::install_github("bocinsky/PaleoCAR@2.0")
library(FedData)
library(PaleoCAR)
pkg_test("parallel")
pkg_test("rgdal")
pkg_test("raster")
pkg_test("rgeos")

# Suppress use of scientific notation
options(scipen=999)

# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+07,maxmemory=2e+08)

## Set the calibration period
# Here, I use a 60 year period ending at 1983 
# to maximize the number of dendro series.
calibration.years <- 1924:1983

## Set the retrodiction years
# Here, we are setting it for 1--2000,
# for consistency with Bocinsky & Kohler 2014
prediction.years <- 1:2000

## Set the directory of the signal geotiffs on which you want to run PaleoCAR
signal.dir <- "../DATA/PRISM/GDD_may_sept_demosaic"

## Set the output directory
out.dir <- "../DATA/PaleoCAR"

# The Four Corners states
if(!dir.exists("../DATA/NATIONAL_ATLAS/statesp010g")){
  dir.create("../DATA/NATIONAL_ATLAS/", showWarnings = F, recursive = T)
  download.file("http://dds.cr.usgs.gov/pub/data/nationalatlas/statesp010g.shp_nt00938.tar.gz", destfile="../DATA/NATIONAL_ATLAS/statesp010g.shp_nt00938.tar.gz", mode='wb')
  untar("../DATA/NATIONAL_ATLAS/statesp010g.shp_nt00938.tar.gz", exdir="../DATA/NATIONAL_ATLAS/statesp010g")
}
states <- readOGR("../DATA/NATIONAL_ATLAS/statesp010g", layer='statesp010g')
states <- states[states$NAME %in% c("Colorado","Utah","New Mexico","Arizona"),]
states <- rgeos::gUnaryUnion(states)

# Create a 10-degree buffer around the 4C states
treePoly <- suppressWarnings(rgeos::gBuffer(states, width=10, quadsegs=1000))

# Extract the Four Corners standardized tree-ring chronologies
ITRDB <- get_itrdb(template=treePoly, label="SKOPE_4CORNERS_PLUS_10DEG", raw.dir = "../DATA/ITRDB/RAW/ITRDB/", extraction.dir = "../DATA/ITRDB/EXTRACTIONS/ITRDB/", recon.years=prediction.years, calib.years=calibration.years, measurement.type="Ring Width", chronology.type="Standard")

# Load the annual chunked raster bricks
chunks.files <- list.files(signal.dir, full.names=T)

## BEGIN PARALLELIZATION!
process.brick <- function(brick.file){
  dir.create(out.dir, recursive=T, showWarnings=F)
  the.brick <- raster::brick(brick.file)
  if(all(is.na(the.brick[]))) return()
  # These bricks are for the whole PRISM time period (1895--2013)
  # Get only calibration years
  the.brick <- raster::subset(the.brick,which(1896:2013 %in% calibration.years))
  names(the.brick) <- calibration.years
  
  junk <- PaleoCAR::paleoCAR.batch(predictands=the.brick, label=basename(tools::file_path_sans_ext(brick.file)), out.dir=out.dir, calibration.years=calibration.years, prediction.years=prediction.years, chronologies=ITRDB, meanVar="chained", asInt=T, min.width=5, floor=0, verbose=F, force.redo=F, generate.reconstruction=T,return.objects=F)

  rm(junk)
  gc();gc()
  return()
}

## PARALLEL RUN
# Set the number of cores on which you want this to run!
this.working <- getwd()
cl <- makeCluster(20)
clusterExport(cl, varlist=c("signal", "out.dir", "calibration.years", "prediction.years", "ITRDB", "this.working"))
clusterEvalQ(cl, {
  # Set the working directory to the location of R-scripts
  setwd(this.working)
  # And libraries
  library(PaleoCAR)
  pkg_test("rgdal")
  pkg_test("raster")
  pkg_test("rgeos")
})
clusterApplyLB(cl, chunks.files, process.brick)

stopCluster(cl)
