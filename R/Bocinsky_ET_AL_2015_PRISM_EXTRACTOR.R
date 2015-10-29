### Script for PRISM climate extraction reported in
### R. Kyle Bocinsky, Johnathan Rush, Keith W. Kintigh, and Timothy A. Kohler, 
### Macrohistory of the Prehispanic Pueblo Southwest: Cycles of Exploration and Exploitation

## This script extracts data from the ~800 m monthly PRISM dataset
## for the SKOPE study area: roughly the total extent of the four-corners states.
## It first defines the extent, chops thattest extent into 800 m (30 arc-second) PRISM cells,
## and subdivides the extent into 14,440 (120x120) cell chunks for computation.
## (the chunks are 1x1 degree).
## These chunks are saved for later computation.

## Set the working directory to the directory of this file!
setwd("/Volumes/BOCINSKY_DATA/WORKING/BOCINSKY_ET_AL_2015/R/")

# Load the functions for all analyses below
# install.packages("devtools", dependencies=T, repos = "http://cran.rstudio.com")
# update.packages(ask=F, repos = "http://cran.rstudio.com")
devtools::install_github("bocinsky/FedData@2.0.0")
library(FedData)
pkg_test("parallel")
pkg_test("sp")
pkg_test("raster")
pkg_test("Hmisc")

# Create an output directory
dir.create("../DATA", showWarnings = F, recursive = T)

# Load all the auxillary functions
all.functions <- lapply(list.files("./src",full.names=T),source)

# Suppress scientific notation
options(scipen=999)

# Force Raster to load large rasters into memory
raster::rasterOptions(chunksize=2e+08,maxmemory=2e+09)

# This MUST point at an original LT81 dataset available from the PRISM climate group
# http://www.prism.oregonstate.edu
PRISM800.DIR <- "/Volumes/DATA/PRISM/LT81_800M/"

# Specify a directory for extraction
EXTRACTION.DIR <- "../DATA/PRISM/"

# The climate parameters to be extracted
types <- c("ppt","tmin","tmax")

##### BEGIN RAW DATA EXTRACTION #####
# Create data output directory if it doesn't already exist
dir.create(EXTRACTION.DIR, showWarnings = F, recursive = T)

# (Down)Load the states shapefile form the National Atlas
if(!dir.exists("../DATA/NATIONAL_ATLAS/statesp010g")){
  dir.create("../DATA/NATIONAL_ATLAS/", showWarnings = F, recursive = T)
  download.file("http://dds.cr.usgs.gov/pub/data/nationalatlas/statesp010g.shp_nt00938.tar.gz", destfile="../DATA/NATIONAL_ATLAS/statesp010g.shp_nt00938.tar.gz", mode='wb')
  untar("../DATA/NATIONAL_ATLAS/statesp010g.shp_nt00938.tar.gz", exdir="../DATA/NATIONAL_ATLAS/statesp010g")
}
states <- readOGR("../DATA/NATIONAL_ATLAS/statesp010g", layer='statesp010g')
states <- states[states$NAME %in% c("Colorado","Utah","New Mexico","Arizona"),]
# Transform to the CRS of the PRISM data
states <- sp::spTransform(states, sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# Get the extent
extent.states <- raster::extent(states)
# Floor the minimums, ceiling the maximums
extent.states@xmin <- floor(extent.states@xmin)
extent.states@ymin <- floor(extent.states@ymin)
extent.states@xmax <- ceiling(extent.states@xmax)
extent.states@ymax <- ceiling(extent.states@ymax)

# Extract and crop all types
type.stacks <- lapply(types, function(type){
  dir.create(paste0(EXTRACTION.DIR,type,'/'), showWarnings = F, recursive = T)
  
  # Get all file names
  monthly.files <- list.files(paste0(PRISM800.DIR,type), recursive=T, full.names=T)
  
  # Trim to only file names that are rasters
  monthly.files <- grep("*\\.bil$", monthly.files, value=TRUE)
  monthly.files <- grep("spqc", monthly.files, value=TRUE, invert=T)
  monthly.files <- grep("/cai", monthly.files, value=TRUE)
  
  # Generate the raster stack
  type.list <- raster::unstack(raster::stack(monthly.files,native=F,quick=T))
  
  cl <- parallel::makeCluster(8)
  system.time(type.list.cropped <- parallel::parLapply(cl,type.list,function(type,PRISM800.DIR,EXTRACTION.DIR,template,rast){
    layer.name <- basename(raster::filename(rast))
    yearmonth <- gsub(".*_","",layer.name)
    yearmonth <- gsub(".bil","",yearmonth)
    year <- as.numeric(substr(yearmonth,1,4))
    month <- as.numeric(substr(yearmonth,5,6))
    if(file.exists(paste0(EXTRACTION.DIR,type,'/','Y',sprintf("%04d", year),'M',sprintf("%02d", month),'.tif'))) return()
    out.rast <- round(raster::crop(rast,template)*ifelse(type=='ppt',1,10))
    return(raster::writeRaster(out.rast,file=paste0(EXTRACTION.DIR,type,'/','Y',sprintf("%04d", year),'M',sprintf("%02d", month),'.tif'), datatype=ifelse(type=='ppt',"INT2U","INT2S"), options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE))
  },template=extent.states, type=type, PRISM800.DIR=PRISM800.DIR, EXTRACTION.DIR=EXTRACTION.DIR))
  parallel::stopCluster(cl)
  
  type.list.cropped <- raster::stack(list.files(paste0(EXTRACTION.DIR,type,'/'), full.names=T), native=F, quick=T)
  
  return(type.list.cropped)
})

names(type.stacks) <- types


# Water-year precipitation
ppt.water_year <- annualize_prism_monthly(prism.brick=type.stacks[['ppt']], months=c(-2:9), fun='sum', out_dir=paste0(EXTRACTION.DIR,"PPT_water_year/"))
ppt.water_year.chunks <- demosaic(raster_brick=ppt.water_year, corners=extent.states.SW.corners, out_dir=paste0(EXTRACTION.DIR,"PPT_water_year_demosaic/"))
rm(ppt.water_year); rm(ppt.water_year.chunks); gc(); gc()

# May--Sept GDD
dir.create(paste0(EXTRACTION.DIR,"gdd/"), recursive=T, showWarnings=F)
gdd.monthly <- calc_gdd_monthly(tmin_brick=type.stacks[['tmin']], tmax_brick=type.stacks[['tmax']], t.base=10, t.cap=30, multiplier=10, to_fahrenheit=T, output.dir=paste0(EXTRACTION.DIR,"gdd/"))
gdd.may_sept <- annualize_prism_monthly(prism.brick=gdd.monthly, months=c(5:9), fun='sum', out_dir=paste0(EXTRACTION.DIR,"GDD_may_sept/"))
gdd.may_sept.chunks <- demosaic(raster_brick=gdd.may_sept, corners=extent.states.SW.corners, out_dir=paste0(EXTRACTION.DIR,"GDD_may_sept_demosaic/"))
rm(gdd.may_sept); rm(gdd.may_sept.chunks); gc(); gc()
