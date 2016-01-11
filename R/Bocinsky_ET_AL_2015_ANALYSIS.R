### Script for all analyses reported in
### R. Kyle Bocinsky, Johnathan Rush, Keith W. Kintigh, and Timothy A. Kohler, 
### Exploration and Exploitation in the Macrohistory of the Prehispanic Pueblo Southwest

### *** DUE TO THEIR SENSITIVE NATURE, TREE-RING DATA ARE AVAILABLE ONLY WITH PERMISSION OF THE CORRESPONDING AUTHORS *** ###
### Contact Timothy A. Kohler (tako@wsu.edu) for access through the Digital Archaeological Record (www.tdar.org).

## This script requires the installation of several command line tools in order to function
## correctly, and has only been tested on a Mac running OS 10.10.
## REQUIRED COMMAND-LINE TOOLS:
## ghostscript ("gs" command)
## imagemagick ("convert" command)
## ffmpeg

## Set the working directory to the directory of this file!
setwd("/Volumes/BOCINSKY_DATA/PUBLICATIONS/SUBMITTED/BOCINSKY_ET_AL_2015/R/")

## Use the checkpoint package to ensure successful execution.
## This package installs temporary versions of the necessary packages
## that conform to the versions available on the run-date of these
## analyses, September 1, 2015.
## Initial installation of the packages might take a few minutes.
# install.packages("checkpoint")
# checkpoint::checkpoint("2015-10-09", checkpointLocation = "./")

## Load the required packages
library("FedData")
pkg_test("bocinsky/PaleoCAR")
pkg_test("data.table")
pkg_test("gtools")
pkg_test("RSQLite")
pkg_test("parallel")
pkg_test("snow")
pkg_test("raster")
pkg_test("rgdal")
pkg_test("rgeos")
pkg_test("png")
pkg_test("rasterVis")
pkg_test("rgl")
pkg_test("psych")
pkg_test("spatstat")
pkg_test("zoo")
pkg_test("scales")
pkg_test("maptools")
pkg_test("RColorBrewer")


# Create an output directory
dir.create("../OUTPUT", showWarnings = F, recursive = T)

# Load all the auxillary functions
all.functions <- lapply(list.files("./src",full.names=T),source)

# Suppress scientific notation
options(scipen=999999)

# Force Raster to load large rasters into memory
raster::rasterOptions(chunksize=2e+08,maxmemory=2e+09)

# Define the final study area
SWUS.extent <- extent(-113,-105,32,38)
SWUS.poly <- polygon_from_extent(SWUS.extent,proj4string="+proj=longlat +datum=WGS84")

# Colors for plots: blue, purple, red
level.colors <- c("#e41a1c","#984ea3","#377eb8")
# and a color ramp interpolation between them
level.colors.ramp <- colorRampPalette(level.colors)(101)

# A smoothing distribution to use throughout the analyses
# This is a 21-year wide Gaussian distribution with a 
# mean of 0 and standard deviation of 5 years
dist <- dnorm(seq(-10,10,1), sd=5)

# A function that counts occurances of a number in a vector of expected range,
# and fills the rest of the range with zeros
my.counter <- function(v,fill.range){
  v.dt <- data.table(v=v)
  v.dt.counts <- v.dt[,.N,by=v]
  setorder(v.dt.counts,v)
  v.dt.counts <- merge(v.dt.counts,data.table(v=fill.range),all=T,by="v")
  v.dt.counts$N[is.na(v.dt.counts$N)] <- 0
  return(v.dt.counts[['N']])
}

## A function to calculate the standard error of the mean
se <- function(x){
  x <- na.omit(x)
  return(sd(x)/sqrt(length(x)))
}

##### TREE-RING ANALYSIS #####
# connect to tree-ring database
# These data are available through the Digital Archaeological Record (www.tdar.org) with the permission
# of the corresponding authors.
sites <- data.table(read.csv("../DATA/BOCINSKY_ET_AL_2015_TREE_RINGS.csv", stringsAsFactors = F))
sites.sp <- SpatialPointsDataFrame(coords=sites[,.(LON,LAT)], data=sites, proj4string=CRS("+proj=longlat +ellps=WGS84"))
# Only include sites in the SWUS study area
sites.sp <- crop(sites.sp,SWUS.poly, snap="out")
sites <- data.table(sites.sp@data)

### Extract 30 arc-second cell numbers from the PRISM Four Corners master raster
## Generate master raster if not already created
if(!file.exists("../OUTPUT/FOUR_CORNERS.rast.tif")){
  all.files <- list.files("../DATA/PaleoCAR/GDD_may_sept_demosaic", pattern="recon", full.names=T)
  all.rasts <- lapply(all.files,raster)
  FOUR_CORNERS.rast <- do.call(raster::merge, all.rasts)
  # Switch to binary representation
  FOUR_CORNERS.rast <- !is.na(FOUR_CORNERS.rast)
  raster::writeRaster(FOUR_CORNERS.rast, "../OUTPUT/FOUR_CORNERS.rast.tif", datatype="INT1U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE) 
}
FOUR_CORNERS.rast <- raster("../OUTPUT/FOUR_CORNERS.rast.tif")

# How many non-NA cells are there?
sum(FOUR_CORNERS.rast[], na.rm=T) ## 2164084 non-NA cells

# Extract site locations, cells first
sites[,CELL:=raster::extract(FOUR_CORNERS.rast,sites[,.(LON,LAT)],cellnumbers=T)[,'cells']]
# Then cell coordinates
sites[,c("CELL.X","CELL.Y"):=list(xFromCell(FOUR_CORNERS.rast,CELL),yFromCell(FOUR_CORNERS.rast,CELL))]

# Create a SpatialPointsDataFrame of just the unique cell locations
sites.raster.points <- unique(sites[,c("CELL","CELL.X","CELL.Y"),with=F])
sites.raster.points <- SpatialPointsDataFrame(as.matrix(sites.raster.points[,c("CELL.X","CELL.Y"),with=F]),data=sites.raster.points[,"CELL",with=F],proj4string=CRS(projection(FOUR_CORNERS.rast)))

## Print some stats about the tree-ring database
nrow(sites) # 32,761 dates
length(unique(sites[['Seq_Number']])) # 1,176 unique sites
nrow(sites[Outer_Date_AD %in% 500:1400]) # 29,311 dates in 500--1400
nrow(sites[Outer_Date_AD %in% 500:1400 & C_Level %in% 3]) # 10,489 cutting dates in 500--1400
nrow(sites[Outer_Date_AD %in% 500:1400 & C_Level %in% 2]) # 3,661 near-cutting dates in 500--1400
nrow(sites[Outer_Date_AD %in% 500:1400 & C_Level %in% 1]) # 15,161 cutting or near-cutting dates in 500--1400
length(unique(sites[Outer_Date_AD %in% 500:1400][['Seq_Number']])) # 1,002 unique sites in 500:1400
100*length(unique(sites[Outer_Date_AD %in% 500:1400 & C_Level %in% 2:3][['Seq_Number']]))/length(unique(sites[Outer_Date_AD %in% 500:1400][['Seq_Number']])) # 70 percent of sites have a cutting date

## Get the counts of cutting versus non-cutting dates in each cell
cutting.cells <- sites[C_Level %in% c(2,3),.N,by=.(Outer_Date_AD,CELL)]
non_cutting.cells <- sites[C_Level  == 1,.N,by=.(Outer_Date_AD,CELL)]
setkey(cutting.cells,Outer_Date_AD,CELL)
setkey(non_cutting.cells,Outer_Date_AD,CELL)
# Merge the cell-level data
all.cells <- merge(cutting.cells,non_cutting.cells, all=T)
setnames(all.cells,c("Outer_Date_AD","CELL","cutting","non_cutting"))
all.cells[["cutting"]][is.na(all.cells[["cutting"]])] <- 0
all.cells[["non_cutting"]][is.na(all.cells[["non_cutting"]])] <- 0
# Calculate the total dates
all.cells[,TOTAL:=cutting+non_cutting]
# create size classes for movie
all.cells[,SIZE:=scales::rescale(log(TOTAL),to=c(0.02,0.1))]
# calculate the proportion of cutting dates
all.cells[,P_TOTAL:=cutting/TOTAL]
# assign to color ramp
all.cells[,COLOR:=colorRampPalette(rev(level.colors))(101)[round(100*P_TOTAL)+1]]

# Create smoothed, stacked probablity distributions of counts of cutting and non-cutting dates
dates.stack <- do.call(rbind,list(my.counter(sites[C_Level %in% c(2,3)][["Outer_Date_AD"]], 1:2000),
                                  my.counter(sites[C_Level == 1][["Outer_Date_AD"]], 1:2000)))
dates.stack.smooth <- t(apply(dates.stack,1,function(x){filter(x=x, filter=dist)}))
dates.stack.smooth[is.na(dates.stack.smooth)] <- 0
dates.stack <- apply(dates.stack,2,cumsum)
dates.stack.smooth <- apply(dates.stack.smooth,2,cumsum)

# Create smoothed, stacked probability distributions of counts of cells containing varying 
# proportions of cutting dates
cells.stack <- do.call(rbind,
                       lapply(rev(seq(0,1,by=0.01)), function(prob){
                         temp <- all.cells[P_TOTAL >= prob & P_TOTAL < (prob+0.01)][["Outer_Date_AD"]]
                         return(my.counter(temp, 1:2000))
                       }))

cells.stack.smooth <- t(apply(cells.stack,1,function(x){filter(x=x, filter=dist)}))
cells.stack.smooth[is.na(cells.stack.smooth)] <- 0
cells.stack <- apply(cells.stack,2,cumsum)
cells.stack.smooth <- apply(cells.stack.smooth,2,cumsum)

## Generate Pecos periods and subperiods based on number of cells with cutting dates
if(!file.exists("../OUTPUT/pueblo.subperiod.breaks.Rds")){
  ## Calculate inflection points of probability density distributions with bandwidths
  ## ranging between 0.01 and 50, by 0.01
  test <- unlist(lapply(seq(0.01,50,by=0.01),function(i){
    cutting.seq <- density(all.cells[cutting > 0 & non_cutting==0 & Outer_Date_AD %in% 500:1400][["Outer_Date_AD"]], from=500, to=1400, n=901, bw=i, kernel='g')$y
    cutting.seq.inflections <- which(diff(sign(diff(diff(cutting.seq)))) %in% c(-2,2))+501
    cutting.seq.inflections[c(1,length(cutting.seq.inflections))] <- c(500,1400)
    return(cutting.seq.inflections)
  }))
  # Calculate new density distribution of inflection points
  test.density <- density(test,from=500, to=1400, n=901, bw=5)$y
  # Get the peaks of the distribution
  test.peaks <- which(diff(sign(diff(test.density))) %in% c(-2))
  # Extract the top ten peaks
  tester <- sort((test.peaks[order(test.density[test.peaks],decreasing = T)]+501)[1:10]) # 602  700  787  892 1036 1144 1199 1285 1365 1400
  # A function to round to the nearest "base" level
  mround <- function(x,base){ 
    base*round(x/base) 
  } 
  # Round to the nearest 5 years, and drop the 1365 inflection
  pueblo.subperiod.breaks <- c(500,mround(tester,5))[-10] # 500  600  700  790  890 1035 1145 1200 1285 1400
  # Save!
  saveRDS(pueblo.subperiod.breaks,"../OUTPUT/pueblo.subperiod.breaks.Rds")  
}
pueblo.subperiod.breaks <- readRDS("../OUTPUT/pueblo.subperiod.breaks.Rds")
# Name the subperiods for later
period.labels <- paste0( c("BM III Exploration, AD ","BM III Exploitation, AD ","P I Exploration, AD ","P I Exploitation, AD ","P II Exploration, AD ","P II Exploitation, AD ","P III Exploration, AD ","P III Exploitation, AD ","P IV, AD "),pueblo.subperiod.breaks[0-(length(pueblo.subperiod.breaks))]+1,"–",pueblo.subperiod.breaks[-1])

### Tree-Ring Robusticity Analysis ###
## Here, we perform a resampling exercise in which we
## restrict each site to one tree-ring date and see if the two
## distributions hold.
## The `Seq_Number` variable is a site-level indicator.
# Split tree-ring data by site
sites.slim <- sites[,.(Seq_Number,Outer_Date_AD,C_Level,CELL)]
if(!file.exists("../OUTPUT/sites_montecarlo.Rds")){
  sites.mc.1 <- rbindlist(mclapply(1:999,function(i){
    out <- sites.slim[,.SD[sample(.N,1)],by = Seq_Number]
    out[,Run:=i]
  }, mc.cores=8))
  saveRDS(sites.mc.1,"../OUTPUT/sites_montecarlo.Rds",compress="xz")
  rm(sites.mc.1); gc(); gc()
}
sites.montecarlo <- readRDS("../OUTPUT/sites_montecarlo.Rds")

### Tree-ring clustering analysis ###
## Here, we test whether spatial clustering of cells containing tree-ring
## dates is significant through time, and we estimate the domain radius,
## or the radius of maximal clustering through time.
## We do so by using Ripley's K (and the derived L and H statistics).
# First, get an annual level cell data, and convert to 
# Lambert's Conic Conformal projection for accurate distance calculations (in meters).
cells <- unique(sites[,.(Outer_Date_AD,CELL.X,CELL.Y)])
cells <- SpatialPointsDataFrame(coords=as.matrix(cells[,.(CELL.X,CELL.Y)]),data=cells[,.(Outer_Date_AD)],proj4string = CRS("+proj=longlat +ellps=WGS84"))
cells <- spTransform(cells,CRS("+proj=lcc +lat_1=37 +lon_0=-109.045225"))
# Convert to a point-pattern object, as used by package spatstat
cells <- as(cells, "ppp")
# Create a vecot of distances over which to calculate the K statistic
distances <- seq(0,dist(cbind(cells$window$xrange,cells$window$yrange)),1000)
# Calculate the K statistic over all points, and use 199 random replicates
# to estimate the L statistic and confidence envelope.
if(!file.exists("../OUTPUT/hats.out.Rds")){
  hats.out <- mclapply(1:2000,function(year){
    this.sites <- subset(cells,marks==year)
    test <- suppressWarnings(envelope(this.sites,Lest,nsim=199,rank=5,global=T,verbose=F))
    return(list(YEAR=year, NPOINTS=npoints(this.sites), TEST=test))
  }, mc.cores=8)
  saveRDS(hats.out,"../OUTPUT/hats.out.Rds")
}
hats.out <- readRDS("../OUTPUT/hats.out.Rds")
# Then, calculate the H statistic [ H(r) = L(r) - r ],
# and estimate the domain radius using the method defined in
# Kiskowski, Hancock, and Kenworthy 2009 (see Supplementary Materials and Methods).
# These are graphed later.
hats.results <- lapply(hats.out,function(hats){
  
  dists.sig <- hats$TEST$obs - hats$TEST$hi>0
  dists.h <- hats$TEST$obs - hats$TEST$r
  rad.max.h <- hats$TEST$r[min(which(dists.h==max(dists.h)))]
  
  if(is.na(dists.sig[1]) | hats$NPOINTS<=2){
    dists.h.smooth.prime.min <- NA
    rad.max.agg.sig <- F
  }else{
    dists.h.smooth.prime <- predict(smooth.spline(x=hats$TEST$r, y=dists.h, spar=0.5),deriv=1)$y
    min.idxs <- which((-1-dists.h.smooth.prime)>0)
    if(!is.empty(min.idxs)){
      dists.h.smooth.prime.min <- hats$TEST$r[min(min.idxs)-1]/2
    }else{
      dists.h.smooth.prime.min <- hats$TEST$r[min(which((-1-dists.h.smooth.prime)==max((-1-dists.h.smooth.prime))))-1]/2
    }
    
    if(is.empty(dists.h.smooth.prime.min)) dists.h.smooth.prime.min <- NA
    rad.max.agg.sig <- dists.sig[min(which(dists.h==max(dists.h)))]
  }
  
  return(list(rad=dists.h.smooth.prime.min,sig=rad.max.agg.sig))
})


##### CONCLUDE PRIMARY TREE-RING ANALYSIS #####




##### NICHE ANALYSIS #####
## By this point, the PaleoCAR reconstruction has been 
## performed for the entire study area.
## Here, we calculate the niche for each year, taken to be 
## water-year precipitation >= 30cm and 
## growing-season GDD >= 1800 (1000 GDD in centigrade)
# List tiles in data directory for PPT; we'll use these names to align the corresponding tiles
all.tiles <- list.files("../DATA/PaleoCAR/PPT_water_year_demosaic",pattern="recon",include.dirs=F)
all.tiles <- all.tiles[grep("nc4",all.tiles,invert=T)]
# Create an output directory for niche tiles
dir.create("../OUTPUT/NICHE_demosaic/", showWarnings = F, recursive = T)
# Calculate niche for each tile
for(tile in all.tiles){
  # Only recalculate if not already calculated!
  if(file.exists(paste0("../OUTPUT/NICHE_demosaic/",tile))) next
  tile.PPT <- brick(paste0("../DATA/PaleoCAR/PPT_water_year_demosaic/",tile)) >= 300
  tile.GDD <- brick(paste0("../DATA/PaleoCAR/GDD_may_sept_demosaic/",tile)) >= 1800
  tile.niche <- tile.PPT*tile.GDD
  writeRaster(tile.niche,paste0("../OUTPUT/NICHE_demosaic/",tile), datatype="INT1U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
}

# These are pretty small, so we can mosaic all of them together
if(!file.exists("../OUTPUT/NICHE_FINAL.tif")){
  all.tiles <- list.files("../OUTPUT/NICHE_demosaic",pattern="recon",full.names = T,include.dirs=F)
  system(paste0("gdalbuildvrt -overwrite ../OUTPUT/NICHE_FINAL.vrt ",paste(all.tiles,collapse=" ")))
  system("gdal_translate -ot Byte -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND ../OUTPUT/NICHE_FINAL.vrt ../OUTPUT/NICHE_FINAL.tif")
}
NICHE <- brick("../OUTPUT/NICHE_FINAL.tif")

# Extract just the study area
if(!file.exists("../OUTPUT/NICHE_FINAL_SWUS.tif")){
  NICHE.SWUS <- crop(NICHE, SWUS.extent, snap='out')
  writeRaster(NICHE.SWUS,"../OUTPUT/NICHE_FINAL_SWUS.tif", datatype="INT1U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
}
NICHE.SWUS <- brick("../OUTPUT/NICHE_FINAL_SWUS.tif")

## Calculate the percent of cells in the niche, through time, for the entire landscape
if(!file.exists("../OUTPUT/all.niche.percent.Rds")){
  all.niche.percent <- cellStats(NICHE.SWUS, mean)
  saveRDS(all.niche.percent, "../OUTPUT/all.niche.percent.Rds", compress = "xz")
}
all.niche.percent <- readRDS("../OUTPUT/all.niche.percent.Rds")
# smooth it for eventual plotting
all.niche.percent.smooth <- filter(all.niche.percent, filter=dist)

## Now, we calculate the count of years (of the last four) each cell was in the niche.
## We have to do this by tile or else it's too much to handle.
dir.create("../OUTPUT/NICHE_demosaic_SMOOTH/", showWarnings = F, recursive = T)
if(!file.exists("../OUTPUT/NICHE_FINAL_SMOOTH.tif")){
  all.rasts <- list.files("../OUTPUT/NICHE_demosaic",pattern="recon",full.names=T,include.dirs=F)
  all.rasts <- all.rasts[grepl(paste(106:113,collapse="|"),all.rasts)]
  all.rasts <- all.rasts[grepl(paste(32:37,collapse="|"),all.rasts)]
  
  for(rast in all.rasts){
    if(file.exists(paste0("../OUTPUT/NICHE_demosaic_SMOOTH/", basename(rast)))) next
    rast.brick <- brick(rast)==1
    gc();gc()
    
    # This runs the calculation on multiple cores
    beginCluster(8)
    f1 <- function(x){calc(x,function(v){as.integer(filter(as.logical(v), filter=rep(1,4), sides=1))})}
    system.time(clusterR(rast.brick, f1, filename=paste0("../OUTPUT/NICHE_demosaic_SMOOTH/", basename(rast)), datatype = "INT1U", options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite = T, setStatistics = FALSE, rast=rast))
    endCluster()
    gc();gc()
    
  }
  gc();gc()
  
  all.rasts <- list.files("../OUTPUT/NICHE_demosaic_SMOOTH/",pattern="recon",full.names=T,include.dirs=F)
  all.rasts <- do.call(raster::merge,lapply(all.rasts,brick))
  writeRaster(all.rasts,"../OUTPUT/NICHE_FINAL_SMOOTH.tif", datatype="INT1U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
  rm(all.rasts)
  gc();gc()
}
NICHE.SMOOTH <- brick("../OUTPUT/NICHE_FINAL_SMOOTH.tif")

## Calculate the percent of years in the niche, AD 500--1400
if(!file.exists("../OUTPUT/NICHE_PERCENT.tif")){
  NICHE.PERCENT <- calc(NICHE[[500:1400]], mean)
  writeRaster(NICHE.PERCENT, filename="../OUTPUT/NICHE_PERCENT.tif", datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
  gc();gc()
}
NICHE.PERCENT <- raster("../OUTPUT/NICHE_PERCENT.tif")

## Calculate the percent of years in the niche for each 
## Pecos subperiod, as defined above
if(!file.exists("../OUTPUT/subperiod.niches.Rds")){
  EXP2.NICHES <- lapply(2:length(pueblo.subperiod.breaks),function(i){
    calc(NICHE[[(pueblo.subperiod.breaks[i-1]+1):(pueblo.subperiod.breaks[i])]], mean)
  })
  names(EXP2.NICHES) <- c("BM III Exploration","BM III Exploitation","P I Exploration","P I Exploitation","P II Exploration","P II Exploitation","P III Exploration","P III Exploitation","P IV")
  saveRDS(EXP2.NICHES,"../OUTPUT/subperiod.niches.Rds",compress="xz")
}
EXP2.NICHES <- readRDS("../OUTPUT/subperiod.niches.Rds")

## Extract the prediction error information from the reconstruction.
## We do this for both water-year precipitation and growing-season GDD.
## We have to do it by-tile, or else it is too much to handle.
if(!file.exists("../OUTPUT/PPT_ERROR_500-1400.tif")){
  all.models <- list.files("../DATA/PaleoCAR/PPT_water_year_demosaic",pattern="models",full.names=T,include.dirs=F)
  
  all.models <- all.models[grepl(paste(106:113,collapse="|"),all.models)]
  all.models <- all.models[grepl(paste(32:37,collapse="|"),all.models)]
  cl <- parallel::makeCluster(8)
  # set up each worker.
  clusterEvalQ(cl, {
    library(PaleoCAR)
    library(raster)
    # Force Raster to load large rasters into memory
    raster::rasterOptions(chunksize=2e+09,maxmemory=2e+10)
  })
  # Calculate space-time averages
  out <- clusterApplyLB(cl, all.models, function(model){
    the.models <- readRDS(model)
    rast <- PaleoCAR::errors.paleocar.models.batch(the.models)[[500:1400]]
    return(list(time.avg = mean(rast),space.avg = cellStats(rast,mean)))
  })
  stopCluster(cl)
  
  space.errors <- do.call(raster::merge,lapply(out,"[[","time.avg"))
  time.errors <- colMeans(do.call(rbind,lapply(out,"[[","space.avg")))
  
  writeRaster(space.errors,"../OUTPUT/PPT_ERROR_500-1400.tif", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
  saveRDS(time.errors,"../OUTPUT/PPT_ERROR_500-1400.Rds", compress="xz")
  rm(out,all.models)
  gc();gc()
}
PPT.ERROR.SPACE <- raster("../OUTPUT/PPT_ERROR_500-1400.tif")
PPT.ERROR.TIME <- readRDS("../OUTPUT/PPT_ERROR_500-1400.Rds")

# Then do it all for the GDD data
if(!file.exists("../OUTPUT/GDD_ERROR_500-1400.tif")){
  all.models <- list.files("../DATA/PaleoCAR/GDD_may_sept_demosaic",pattern="models",full.names=T,include.dirs=F)
  
  all.models <- all.models[grepl(paste(106:113,collapse="|"),all.models)]
  all.models <- all.models[grepl(paste(32:37,collapse="|"),all.models)]

  cl <- parallel::makeCluster(8)
  # set up each worker.
  clusterEvalQ(cl, {
    library(PaleoCAR)
    library(raster)
    # Force Raster to load large rasters into memory
    raster::rasterOptions(chunksize=2e+09,maxmemory=2e+10)
  })
  # Calculate space-time averages
  out <- clusterApplyLB(cl, all.models, function(model){
    the.models <- readRDS(model)
    rast <- PaleoCAR::errors.paleocar.models.batch(the.models)[[500:1400]]
    return(list(time.avg = mean(rast),space.avg = cellStats(rast,mean)))
  })
  stopCluster(cl)
  
  space.errors <- do.call(raster::merge,lapply(out,"[[","time.avg"))
  time.errors <- colMeans(do.call(rbind,lapply(out,"[[","space.avg")))
  
  writeRaster(space.errors,"../OUTPUT/GDD_ERROR_500-1400.tif", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
  saveRDS(time.errors,"../OUTPUT/GDD_ERROR_500-1400.Rds", compress="xz")
  rm(out,all.models)
  gc();gc()
}
GDD.ERROR.SPACE <- raster("../OUTPUT/GDD_ERROR_500-1400.tif")
GDD.ERROR.TIME <- readRDS("../OUTPUT/GDD_ERROR_500-1400.Rds")

##### CONCLUDE PRIMARY NICHE ANALYSIS #####



##### BEGIN JOINT NICHE/TREE-RING ANALYSIS #####

## First, we extract the niche record under
## each cell that ever has a tree-ring date
if(!file.exists("../OUTPUT/NICHE.extraction.Rds")){
#   cl <- parallel::makeCluster(8)
#   clusterExport(cl,"NICHE")
#   NICHE.extraction <- parallel::parRapply(cl,sites.raster.points@coords, raster::extract, x=NICHE)
#   stopCluster(cl)
#   NICHE.TEMP <- readAll(NICHE.SWUS)
  NICHE.extraction <- raster::extract(NICHE,sites.raster.points)
  colnames(NICHE.extraction) <- 1:2000
  NICHE.extraction <- data.table::data.table(NICHE.extraction)
  NICHE.extraction <- SpatialPointsDataFrame(as.matrix(unique(sites[,c("CELL.X","CELL.Y"),with=F])),data=NICHE.extraction,proj4string=CRS(projection(FOUR_CORNERS.rast)))
  saveRDS(NICHE.extraction,"../OUTPUT/NICHE.extraction.Rds",compress="xz")
  rm(NICHE.extraction); gc(); gc()
}
all.niche.points <- readRDS("../OUTPUT/NICHE.extraction.Rds")

## Then, create a masking table to identify only those years
## that each cell has a date.
sites.raster.dates <- do.call(rbind,lapply(sites.raster.points$CELL, function(cell){
  mat.out <- matrix(data=0, nrow=1,ncol=length(1:2000))
  mat.out[1,all.cells[CELL==cell]$Outer_Date_AD] <- 1
  return(mat.out)
}))

## Smooth the masking table so that the current year a previous 3 years are included
sites.raster.dates.run <- t(apply(sites.raster.dates, 1, function(x){
  as.numeric(running(x, fun=max, width=4, align="left", pad=T, allow.fewer=T, na.rm=T))
})) 
sites.raster.dates.run[sites.raster.dates.run==0] <- NA

## Finally, generate the "all" and "local" datasets (see Fig. 2 and Supplementary Materials and Methods)
all.data <- all.niche.points@data
local.data <- all.data * sites.raster.dates.run

## What is the niche percent for occupies versus all cells?
niche.percent.swus <- crop(NICHE.PERCENT, SWUS.extent, snap='out')
niche.percent.occupied <- extract(niche.percent.swus,sites.raster.points)

## Calculate the percent of cells in the niche each year, 
## and smooth using a 21-year Gaussian filter with a 5-year standard deviation,
## as defined at the top of this script.
all.data.mean <- filter(apply(all.data,2,mean, na.rm=T),filter=dist)

## Do the same with the local data, 
## and calculate upper and lower confidence intervals
local.data.mean <- filter(apply(local.data,2,mean, na.rm=T),filter=dist)
local.data.se <-  filter(apply(local.data,2,se)*1.96,filter=dist)
local.data.se.lower <- local.data.mean-local.data.se; local.data.se.lower[local.data.se.lower<0] <- 0; local.data.se.lower[is.na(local.data.se.lower) | local.data.se.lower==1] <- 0
local.data.se.upper <- local.data.mean+local.data.se; local.data.se.upper[local.data.se.upper>1] <- 1; local.data.se.upper[is.na(local.data.se.upper) | local.data.se.upper==0] <- 1

## Create a table summarizing niche distribution by subperiod
## Include numbers/proportions of t-r dates
EXP2.NICHES.swus <- lapply(EXP2.NICHES,crop, y=SWUS.extent, snap='out')
all.means <- lapply(EXP2.NICHES.swus,cellStats, mean)
occupied.means <- lapply(lapply(EXP2.NICHES.swus, extract, y=sites.raster.points),mean)
local.means <- lapply(2:length(pueblo.subperiod.breaks), function(i){
  return(mean(local.data.mean[(pueblo.subperiod.breaks[i-1]+1):pueblo.subperiod.breaks[i]]))
})
names(local.means) <- names(occupied.means)
treering.periods <- lapply(2:length(pueblo.subperiod.breaks), function(i){
  cutting <- nrow(sites[Outer_Date_AD %in% (pueblo.subperiod.breaks[i-1]+1):pueblo.subperiod.breaks[i] & C_Level %in% c(2,3)])
  noncutting <- nrow(sites[Outer_Date_AD %in% (pueblo.subperiod.breaks[i-1]+1):pueblo.subperiod.breaks[i] & C_Level %in% 1])
  p.cutting <- cutting/(cutting+noncutting)
  return(list(cutting=cutting, noncutting=noncutting, p.cutting=p.cutting))
})
names(treering.periods) <- names(occupied.means)

period.summaries <- as.data.frame(rbindlist(list(lapply(treering.periods,'[[','cutting'),lapply(treering.periods,'[[','noncutting'),lapply(treering.periods,'[[','p.cutting'),all.means,occupied.means,local.means)))
rownames(period.summaries) <- c('Cutting','Non-cutting','Percent Cutting','MFN - All', 'MFN - Occupied','MFN - Local')
write.csv(period.summaries,"../OUTPUT/period.summaries.csv")


## Summarize the niche distribution by period
## And taking the average over the last four years
f1 <- function(x){as.integer(filter(as.logical(x), filter=rep(1,4), sides=1))}
all.data.four <- t(apply(all.data,1,f1))
local.data.four <- all.data.four * sites.raster.dates.run
period.summaries.four <- lapply(2:length(pueblo.subperiod.breaks), function(i){
  all.describe <- psych::describe(as.integer(as.matrix(all.data.four[,(pueblo.subperiod.breaks[i-1]+1):pueblo.subperiod.breaks[i]])), na.rm=T, type=2)
  local.describe <- psych::describe(as.integer(as.matrix(local.data.four[,(pueblo.subperiod.breaks[i-1]+1):pueblo.subperiod.breaks[i]])), na.rm=T, type=2)
  return(list(all=all.describe,local=local.describe))
})
names(period.summaries.four) <- period.labels
write.csv(do.call(rbind,lapply(period.summaries.four,'[[','all')),"../OUTPUT/period.summaries.all.four.csv")
write.csv(do.call(rbind,lapply(period.summaries.four,'[[','local')),"../OUTPUT/period.summaries.local.four.csv")

# Do the same for the entire landscape
if(!file.exists("../OUTPUT/total.summaries.four.Rds")){
  total.summaries.four <- lapply(2:length(pueblo.subperiod.breaks), function(i){
    total.describe <- psych::describe(as.integer(NICHE.SMOOTH[[(pueblo.subperiod.breaks[i-1]+1):pueblo.subperiod.breaks[i]]][]), na.rm=T, type=2)
    return(total.describe)
  })
  names(total.summaries.four) <- period.labels
  saveRDS(total.summaries.four,"../OUTPUT/total.summaries.four.Rds")
}
total.summaries.four <- readRDS("../OUTPUT/total.summaries.four.Rds")

## Create some series for coloring each year based on the mean percent of cutting dates
## in each cell
P_TOTAL.means <- all.cells[Outer_Date_AD %in% 500:1400,mean(P_TOTAL),by=Outer_Date_AD]
P_TOTAL.means <- merge(P_TOTAL.means,data.table(Outer_Date_AD=500:1400), all=T)[["V1"]]
P_TOTAL.means.interp <- approx(P_TOTAL.means,n=length(P_TOTAL.means)*100)
P_TOTAL.colors <- rev(level.colors.ramp)[round(100*P_TOTAL.means)+1]
P_TOTAL.colors.interp <- rev(level.colors.ramp)[round(100*P_TOTAL.means.interp[['y']])+1]
local.data.mean.interp <- approx(local.data.mean[500:1400],n=length(local.data.mean[500:1400])*100)

## Create a png of the mean, smoothed percent of cells in the niche, through time, with
## 95% conficence intervals colored by the mean percent of cutting dates in each cell.
## This facilitates quick plotting later.
quartz(file='../FIGURES/niche_plot.png', width=5.275, height=1.25, bg='transparent', type='png', family="Helvetica Bold", pointsize=8, dpi=600)
par(mai=c(0,0,0,0))
plot(1, type='n', xlab="", ylab="",xaxs='i',yaxs='i', xlim=c(500,1400), ylim=c(0,1), axes=FALSE, main='')
for(year.plot in 500:1400){
  if(!is.na(P_TOTAL.colors[[year.plot-499]])){
    polygon(x=c(year.plot-0.5,year.plot+0.5,year.plot+0.5,year.plot-0.5),y=c(local.data.se.upper[year.plot],local.data.se.upper[year.plot],local.data.se.lower[year.plot],local.data.se.lower[year.plot]), col=P_TOTAL.colors[[year.plot-499]], border=NA)
  }
}
lines(x=local.data.mean.interp[['x']]+499,y=local.data.mean.interp[['y']], col='white', lwd=1)
dev.off()
plot.png <- readPNG('../FIGURES/niche_plot.png')

##### CONCLUDE JOINT NICHE/TREE-RING ANALYSIS #####


##### FIG_1.pdf #####
# This increases the number of vertices of the polygon for appropriate transformation
SWUS.poly.lam <- SpatialPolygons(list(Polygons(list(Polygon(coordinates(suppressWarnings(spsample(as(SWUS.poly,"SpatialLines"),n=1000000,type="regular"))))),ID="1")),proj4string=CRS(projection(SWUS.poly)))
SWUS.poly.lam <- spTransform(SWUS.poly.lam,CRS("+proj=lcc +lat_1=37 +lon_0=-109.045225"))

# (Down)Load the states shapefile form the National Atlas
if(!dir.exists("../DATA/NATIONAL_ATLAS/statesp010g")){
  dir.create("../DATA/NATIONAL_ATLAS/", showWarnings = F, recursive = T)
  download.file("http://dds.cr.usgs.gov/pub/data/nationalatlas/statesp010g.shp_nt00938.tar.gz", destfile="../DATA/NATIONAL_ATLAS/statesp010g.shp_nt00938.tar.gz", mode='wb')
  untar("../DATA/NATIONAL_ATLAS/statesp010g.shp_nt00938.tar.gz", exdir="../DATA/NATIONAL_ATLAS/statesp010g")
}
states <- readOGR("../DATA/NATIONAL_ATLAS/statesp010g", layer='statesp010g')
states <- states[states$NAME %in% c("Colorado","Utah","New Mexico","Arizona"),]
FOUR.extent <- extent(states)
states <- spTransform(states,"+proj=lcc +lat_1=37 +lon_0=-109.045225")

# Define the plot area as a polygon around the states
sim.poly <- polygon_from_extent(rgeos::gBuffer(states, width=20000, quadsegs=1000))
sim.extent <- extent(sim.poly)

# Define the aspect ratio of the plot, and its dimensions
plot.ratio <- (ymax(sim.extent)-ymin(sim.extent))/(xmax(sim.extent)-xmin(sim.extent))
fig.width <- 4.6
legend.height <- 0.0
plot.width <- fig.width
plot.height <- plot.width/plot.ratio
fig.height <- plot.height + legend.height

# define some colors for the niche percentage
colors <- paste0(colorRampPalette(brewer.pal(9,"Greens"))(11),c("00","0D","1A","26","33","40","4D","59","66","73","80")[11])
color.breaks <- seq(from=0,to=1,by=0.1)
# Create a png of the niche percentage
if(!file.exists("../FIGURES/NICHE.PERCENT.png")){
  NICHE.PERCENT <- raster("../OUTPUT/NICHE_PERCENT.tif")
  NICHE.PERCENT <- projectRaster(NICHE.PERCENT,crs=CRS("+proj=lcc +lat_1=37 +lon_0=-109.045225"))
  NICHE.PERCENT <- mask(NICHE.PERCENT,states)
  
  quartz(file='../FIGURES/NICHE.PERCENT.png', width=6*plot.ratio, height=6, antialias=FALSE, bg="transparent", type='png', family="Gulim", pointsize=8, dpi=300)
  par(mai=c(0,0,0,0))
  plot(1, type='n', xlab="", ylab="", xlim=c(xmin(sim.extent),xmax(sim.extent)), ylim=c(ymin(sim.extent),ymax(sim.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
  plot(NICHE.PERCENT, zlim=c(0,1), maxpixels=ncell(NICHE.PERCENT), breaks=color.breaks, col=colors, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
  dev.off()
}
NICHE.PERCENT.PNG <- readPNG('../FIGURES/NICHE.PERCENT.png')

# Create a png for the elevation background
# The GTOPO30 data are available from 
if(!file.exists("../FIGURES/FOUR_Background.png")){
  
  FOUR.GTOPO <- lapply(c("w100n40","w100n90","w140n40","w140n90"),function(f){
    if(!file.exists(paste0("../DATA/GTOPO30/",f,".zip"))){
      download.file(paste0("http://www.webgis.com/GTOPO30/",f,".zip"),destfile=paste0("../DATA/GTOPO30/",f,".zip"))
      unzip(paste0("../DATA/GTOPO30/",f,".zip"), exdir="../DATA/GTOPO30/")
    }
    return(raster(paste0("../DATA/GTOPO30/",toupper(f),".DEM")))
  })
  
  FOUR.GTOPO <- raster::crop(do.call(raster::merge,FOUR.GTOPO),FOUR.extent, snap='out')*2
  slope <- raster::terrain(FOUR.GTOPO, opt='slope')
  aspect <- raster::terrain(FOUR.GTOPO, opt='aspect')
  FOUR.GTOPO.hill <- raster::hillShade(slope, aspect, 40, 230)
  
  FOUR.GTOPO.hill <- projectRaster(FOUR.GTOPO.hill,crs=CRS("+proj=lcc +lat_1=37 +lon_0=-109.045225"))
  FOUR.GTOPO.hill <- mask(FOUR.GTOPO.hill,states)
  
  quartz(file='../FIGURES/FOUR_Background.png', width=6*plot.ratio, height=6, antialias=FALSE, bg="white", type='png', family="Gulim", pointsize=8, dpi=300)
  par(mai=c(0,0,0,0))
  plot(1, type='n', xlab="", ylab="", xlim=c(xmin(sim.extent),xmax(sim.extent)), ylim=c(ymin(sim.extent),ymax(sim.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
  plot(FOUR.GTOPO.hill, maxpixels=ncell(FOUR.GTOPO.hill), col=grey(60:100/100), useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
  dev.off()
  rm(FOUR.GTOPO,FOUR.GTOPO.hill,slope,aspect); gc(); gc()
}
FOUR.background <- readPNG('../FIGURES/FOUR_Background.png')

sites.sp.slim <- sites.sp[sites.sp$Outer_Date_AD %in% 500:1400 & !duplicated(sites.sp@coords),]


pdf(file='../FIGURES/FIG_1.pdf', width=fig.width, height=fig.height, bg="white", pointsize=8, version="1.7")
par(bg='white',fg='black',col.lab='black', col.main='black', col.axis='black', font=2, lend='round',ljoin='round')

par(mai=c(0,0,0,0), lend='round', ljoin='round', xpd=T)
plot(1, type='n', xlab="", ylab="",xaxs='i',yaxs='i', xlim=c(0,1), ylim=c(0,1), axes=FALSE, main='')

par(mai=c(legend.height,0,0,0), xpd=F, new=T)
plot(1, type='n', xlab="", ylab="",xlim=c(xmin(sim.extent),xmax(sim.extent)), ylim=c(ymin(sim.extent),ymax(sim.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
rasterImage(FOUR.background, xleft=xmin(sim.extent), xright=xmax(sim.extent), ybottom=ymin(sim.extent), ytop=ymax(sim.extent), interpolate=F)
rasterImage(NICHE.PERCENT.PNG, xleft=xmin(sim.extent), xright=xmax(sim.extent), ybottom=ymin(sim.extent), ytop=ymax(sim.extent), interpolate=F)

points(spTransform(sites.sp.slim,CRS("+proj=lcc +lat_1=37 +lon_0=-109.045225")), pch=19, cex=0.4)

plot(states, add=T)
plot(SWUS.poly.lam, add=T, lty=3)

label.coords <- coordinates(states)
text(label.coords, labels=states$STATE, cex=1)

inch.x <- (sim.extent@xmax-sim.extent@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
inch.y <- (sim.extent@ymax-sim.extent@ymin)/(fig.height-par('mai')[1]-par('mai')[3])
scalebar.new(d=100000, lonlat=F, cex=1, font=2, side='right',lab.side='right', height=0.05*inch.y, label="100 km", line.offset=c(0.05*inch.x,0.05*inch.y), xy=c(xmin(sim.extent),ymin(sim.extent)), lwd=4, lend=1)

dev.off()
distill('../FIGURES/FIG_1.pdf')
# distill('../FIGURES/FIG_1_EDITED.pdf') ## Manually edited the map that appears in publication.


##### FIG_2.pdf #####
phi <- (1+sqrt(5))/2
mai <- c(0.25,0.25,0,0)
nplots <- 4
fig.width <- 4.6
fig.height <- nplots*fig.width/(phi^2)
between <- 0.2
margins <- 0.5
plot.height <- (fig.height - (margins*1.5) - (between*(nplots-1)))/nplots
legend.cex <- 0.75

quartz(file="../FIGURES/FIG_2.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=8, dpi=600)

this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,125), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
polygon(x=c(500:1400,1400:500),y=c(dates.stack.smooth[2,500:1400],rep(0,length(dates.stack.smooth[2,500:1400]))),col=level.colors[3], border=NA)
polygon(x=c(500:1400,1400:500),y=c(dates.stack.smooth[1,500:1400],rep(0,length(dates.stack.smooth[1,500:1400]))),col=level.colors[1], border=NA)

this.plot <- 2
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,40), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
for(i in nrow(cells.stack.smooth):1){
  polygon(x=c(500:1400,1400:500),y=c(cells.stack.smooth[i,500:1400],rep(0,length(cells.stack.smooth[i,500:1400]))),col=level.colors.ramp[i], border=NA)
}

this.plot <- 3
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=F, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,1), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
rasterImage(plot.png, 500,0,1400,1, interpolate=F)
lines(x=500:1400, y=all.niche.percent.smooth[500:1400])
lines(x=500:1400, y=all.data.mean[500:1400], lty=3)


this.plot <- 4
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=F, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,100), xlim=c(500,1400), axes=FALSE, main='', xpd=T)

rads <- sapply(hats.results,'[[','rad')
sigs <- sapply(hats.results,'[[','sig')

rads.noNA <- na.approx(rads, na.rm=F)

hats.results.smoothed <- filter(rads.noNA,filter=dist)
hats.results.smoothed[!sigs] <- NA

lines(y=rads[500:1400]/1000, x=500:1400, lwd=0.5)
lines(y=hats.results.smoothed[500:1400]/1000, x=500:1400, col=level.colors[1], lwd=2)

par(mai=c(margins,margins*1.25,margins*0.5,margins*0.5), xpd=F, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,1),xlim=c(500,1400), axes=FALSE, main='')
abline(v=pueblo.subperiod.breaks[c(3,5,7,9)], lwd=1.5,lty=2)
abline(v=pueblo.subperiod.breaks[c(2,4,6,8)], lwd=1,lty=3)
text(x=diff(pueblo.subperiod.breaks[c(1,3,5,7,9,10)])/2 + pueblo.subperiod.breaks[c(1,3,5,7,9)],y=par('usr')[4],labels=c("BM III","P I","P II","P III","P IV"),adj=c(0.5,-1), xpd=T)

## LEGENDS and AXES
this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,125), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
legend("topleft", fill=c(level.colors[c(1,3)]), legend=c("Cutting and near-cutting","Non-cutting"), bg='#FFFFFFBF', box.col=NA, border=NA)
axis(1, at=seq(500,1400,100), labels=NA, cex.axis=0.8)
axis(2, at=seq(0,125,25), cex.axis=0.8)
mtext("Number of dates", side=2, line=2.25)
par(xpd=T)
text(x=375,y=125, labels=LETTERS[this.plot], font=2, cex=1.75)

this.plot <- 2
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,40), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
color.breaks <- seq(from=525,to=875,by=(875-525)/length(level.colors.ramp))
rect(xleft=500, ybottom=22.5, xright=900, ytop=40, col='#FFFFFFBF', border=NA)
rect(xleft=color.breaks[1:(length(color.breaks)-1)], ybottom=32.5, xright=color.breaks[2:length(color.breaks)], ytop=36, col=rev(level.colors.ramp), border=NA)
rect(xleft=color.breaks[1], ybottom=32.5, xright=color.breaks[length(color.breaks)], ytop=36, col=NA, border="black", lwd=0.5)
text(x=seq(525,875,length.out=5),y=32.5,labels=paste0(seq(0,100,25),"%"),adj=c(0.5,2),xpd=T, cex=0.8)
text(x=700, y=32.5, labels="Percent cutting and near-cutting", adj=c(0.5,3.5), xpd=T, cex=1)
axis(2, at=seq(0,40,10), cex.axis=0.8)
mtext("Number of cells", side=2, line=2.25)

axis(1, at=seq(500,1400,100), labels=NA, cex.axis=0.8)
par(xpd=T)
text(x=375,y=40, labels=LETTERS[this.plot], font=2, cex=1.75)

this.plot <- 3
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,1), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
par(lend='butt')
legend("bottomleft", col=c('#FFFFFF00','#FFFFFF00',level.colors.ramp[50]), legend=c(expression("All cells"),expression("Cells with dates"),expression("Local cells ("%+-%'95% CI)')),text.col='#FFFFFF00', lwd=10, bg='#FFFFFFBF', box.col=NA, border=NA)
legend("bottomleft", col=c('black','black','white'), legend=c(expression("All cells"),expression("Cells with dates"),expression("Local cells ("%+-%'95% CI)')), lwd=1, lty=c(1,3,1), bg='#FFFFFF00', box.col=NA, border=NA)
par(lend='round')
axis(2, at=seq(0,1,0.2), labels=seq(0,100,20), cex.axis=0.8)
mtext("Percent in niche", side=2, line=2.25)

axis(1, at=seq(500,1400,100), labels=NA, cex.axis=0.8)

par(xpd=T)
text(x=375,y=1, labels=LETTERS[this.plot], font=2, cex=1.75)

this.plot <- 4
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,100), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
axis(2, at=seq(0,100,20), labels=seq(0,100,20), cex.axis=0.8)
mtext("Domain radius (km)", side=2, line=2.25)

axis(1, at=seq(500,1400,100), labels=seq(500,1400,100), cex.axis=0.8)
mtext("Year AD", side=1, line=2.25)
par(xpd=T)
text(x=375,y=100, labels=LETTERS[this.plot], font=2, cex=1.75)

dev.off()
distill("../FIGURES/FIG_2.pdf")



##### FIG_5.pdf #####
## NICHE COMPARISON MAPS
## First, get the National Elevation Dataset for the study area.
SWUS.extent <- extent(-113,-105,32,38)
SWUS.poly <- polygon_from_extent(SWUS.extent,proj4string="+proj=longlat +datum=WGS84")
EXP2.NICHES.RASTS <- lapply(EXP2.NICHES, crop, y=SWUS.extent, snap='out')
if(!file.exists("../OUTPUT/SWUS.NED.small.tif")){
  SWUS.NED <- get_ned(template=SWUS.poly,label="SWUS",res='1',raw.dir='../DATA/NED/', extraction.dir='../DATA/NED/EXTRACTIONS/')
  SWUS.NED.small <- aggregate(SWUS.NED, 10)
  raster::writeRaster(SWUS.NED.small, "../OUTPUT/SWUS.NED.small.tif", datatype = "FLT4S", options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite = T, setStatistics = FALSE)
}
SWUS.NED.small <- crop(raster("../OUTPUT/SWUS.NED.small.tif"),SWUS.poly,snap='out')*2
slope <- terrain(SWUS.NED.small, opt='slope')
aspect <- terrain(SWUS.NED.small, opt='aspect')
hill <- hillShade(slope, aspect, 40, 230)

## Then, create a png of the elevation background for quicker plotting
# SWUS background
plot.ratio <- 1.093363
fig.width <- 4.6
between <- 0.05
nrow <- 5
ncol <- 2
legend.height <- 0.0
year.height <- 0.0
plot.width <- (fig.width - (ncol-1)*between - year.height)/ncol
plot.height <- plot.width/plot.ratio
fig.height <- (plot.height * nrow) + legend.height + year.height + (between * (nrow-1))

if(!file.exists('../FIGURES/SWUS_Background.png')){
  quartz(file='../FIGURES/SWUS_Background.png', width=6*plot.ratio, height=6, antialias=FALSE, bg="white", type='png', family="Gulim", pointsize=8, dpi=300)
  par(mai=c(0,0,0,0))
  plot(1, type='n', xlab="", ylab="", xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
  plot(hill, maxpixels=ncell(hill), col=grey(60:100/100), useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
  dev.off()
}
SWUS.background <- readPNG('../FIGURES/SWUS_Background.png')

## Finally, create some colors
colors <- paste0(colorRampPalette(brewer.pal(9,"Greens"))(11),c("00","0D","1A","26","33","40","4D","59","66","73","80")[11])
color.breaks <- seq(from=0,to=1,by=0.1)
titles <- names(EXP2.NICHES)
years <- paste0("AD ",sapply(2:length(pueblo.subperiod.breaks),function(i){paste0((pueblo.subperiod.breaks[i-1]+1),"-",(pueblo.subperiod.breaks[i]))}))
cols <- colorRampPalette(level.colors)(101)
cols.breaks <- seq(from=0,to=1,by=0.01)

## Seperate out the period points
points <- lapply(2:length(pueblo.subperiod.breaks),function(i){
  out.cutting <- sites[Outer_Date_AD %in% (pueblo.subperiod.breaks[i-1]+1):(pueblo.subperiod.breaks[i]) & C_Level %in% c(2,3) & C_Level != 1,.N,by=CELL]
  out.non_cutting <- sites[Outer_Date_AD %in% (pueblo.subperiod.breaks[i-1]+1):(pueblo.subperiod.breaks[i]) & !(C_Level %in% c(2,3)) & C_Level == 1,.N,by=CELL]
  out.both <- merge(out.cutting,out.non_cutting,all=T,by="CELL")
  out.both$N.x[is.na(out.both$N.x)] <- 0
  out.both$N.y[is.na(out.both$N.y)] <- 0
  out.both$P.y <- round((out.both$N.y/(out.both$N.x + out.both$N.y))*100)+1
  
  out.both.sp <- sites.raster.points[match(out.both$CELL,sites.raster.points$CELL),]
  out.both.sp$COLOR <- cols[out.both$P.y]
  out.both.sp <- out.both.sp[order(out.both$P.y, decreasing=T),]
  
  return(out.both.sp)
})

## And create the figure
pdf(file='../FIGURES/FIG_5.pdf', width=fig.width, height=fig.height, bg="white", pointsize=8, version="1.7")
par(bg='white',fg='black',col.lab='black', col.main='black', col.axis='black', font=2, lend='round',ljoin='round')

locs <- expand.grid(1:nrow,1:ncol)
locs <- locs[order(locs[,1]),]

par(mai=c(0,0,0,0), lend='round', ljoin='round', xpd=T)
plot(1, type='n', xlab="", ylab="",xaxs='i',yaxs='i', xlim=c(0,1), ylim=c(0,1), axes=FALSE, main='')

for(i in 1:(ncol*nrow)){
  loc <- locs[i,]
  
  if(length(EXP2.NICHES.RASTS)<i) break
  
  par(mai=c(legend.height + (nrow-loc[1])*(plot.height) + (nrow-loc[1])*(between),year.height+((ncol-1)-(ncol-loc[2]))*(plot.width+between),year.height+((loc[1]-1)*(plot.height))+((loc[1]-1)*(between)),(ncol-loc[2])*(plot.width+between)), xpd=F, new=T)
  plot(1, type='n', xlab="", ylab="",xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
  rasterImage(SWUS.background, xleft=xmin(SWUS.extent), xright=xmax(SWUS.extent), ybottom=ymin(SWUS.extent), ytop=ymax(SWUS.extent), interpolate=F)
  
  quartz(file='../FIGURES/temp.png', width=6*plot.ratio, height=6, antialias=FALSE, bg="transparent", type='png', family="Gulim", pointsize=8, dpi=300)
  par(mai=c(0,0,0,0))
  plot(1, type='n', xlab="", ylab="", xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
  plot(EXP2.NICHES.RASTS[[i]], zlim=c(0,1), maxpixels=ncell(EXP2.NICHES.RASTS[[i]]), breaks=color.breaks, col=colors, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
  dev.off()
  rasterImage(readPNG('../FIGURES/temp.png'), xleft=xmin(SWUS.extent), xright=xmax(SWUS.extent), ybottom=ymin(SWUS.extent), ytop=ymax(SWUS.extent), interpolate=F)
  
  plot(points[[i]], col=points[[i]]$COLOR, pch=19, cex=0.5, add=T)
  
  inch.x <- (SWUS.extent@xmax-SWUS.extent@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
  inch.y <- (SWUS.extent@ymax-SWUS.extent@ymin)/(fig.height-par('mai')[1]-par('mai')[3])

  text(x=xmin(SWUS.extent)+(0.075*inch.x), y=ymin(SWUS.extent)+(0.225*inch.y), labels=titles[[i]], col="black", font=2, cex=1, adj=c(0,0), xpd=T)  
  text(x=xmin(SWUS.extent)+(0.075*inch.x), y=ymin(SWUS.extent)+(0.075*inch.y), labels=years[[i]], col="black", font=2, cex=1, adj=c(0,0), xpd=T)  
  
  text(x=xmin(SWUS.extent)+(0.075*inch.x), y=ymax(SWUS.extent)-(0.075*inch.y), labels=LETTERS[[i]], col="black", font=2, cex=2, adj=c(0,1), xpd=T)  
  
}

unlink('../FIGURES/temp.png')

par(mai=c(legend.height + (nrow-loc[1])*(plot.height) + (nrow-loc[1])*(between),year.height+((ncol-1)-(ncol-loc[2]))*(plot.width+between),year.height+((loc[1]-1)*(plot.height))+((loc[1]-1)*(between)),(ncol-loc[2])*(plot.width+between)), xpd=F, new=T)
plot(1, type='n', xlab="", ylab="",xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
inch.x <- (SWUS.extent@xmax-SWUS.extent@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
inch.y <- (SWUS.extent@ymax-SWUS.extent@ymin)/(fig.height-par('mai')[1]-par('mai')[3])
scalebar.new(d=100, cex=1, font=2, side='left',lab.side='left', height=0.05*inch.y, label="100 km", line.offset=c(-0.05*inch.x,0.05*inch.y), xy=c(xmax(SWUS.extent),ymin(SWUS.extent)), lwd=3, lend=1)
par(mai=c(legend.height + (nrow-loc[1])*(plot.height) + (nrow-loc[1])*(between),0.25+year.height+((ncol-1)-(ncol-loc[2]))*(plot.width+between),year.height+((loc[1]-1)*(plot.height))+((loc[1]-1)*(between)),0.25+(ncol-loc[2])*(plot.width+between)), xpd=T, new=T)
plot(1, type='n', xlab="", ylab="",xlim=c(0,1), ylim=c(0,plot.height), xaxs="i", yaxs="i", axes=FALSE, main='')
rect(xleft=color.breaks[1], ybottom=1.5, xright=color.breaks[length(color.breaks)], ytop=1.6, col="gray90", border=NA)
rect(xleft=color.breaks[1:(length(color.breaks)-1)], ybottom=1.5, xright=color.breaks[2:length(color.breaks)], ytop=1.6, col=colors, border=NA)
text(x=seq(0,1,0.2),y=1.4,labels=c("0%","20%","40%","60%","80%","100%"),adj=c(0.5,0.5),xpd=T, cex=0.75)
text(x=0.5, y=1.3, labels="Percent of years in niche", adj=c(0.5,1),xpd=T, cex=1)

rect(xleft=cols.breaks[1], ybottom=0.9, xright=cols.breaks[length(cols.breaks)], ytop=1, col="gray90", border=NA)
rect(xleft=cols.breaks[1:(length(cols.breaks)-1)], ybottom=0.9, xright=cols.breaks[2:length(cols.breaks)], ytop=1, col=rev(cols), border=NA)
text(x=seq(0,1,0.2),y=0.8,labels=c("0%","20%","40%","60%","80%","100%"),adj=c(0.5,0.5),xpd=T, cex=0.75)
text(x=0.5, y=0.7, labels="Percent cutting and near-cutting", adj=c(0.5,1),xpd=T, cex=1)

dev.off()
distill('../FIGURES/FIG_5.pdf')



##### FIG_3.pdf #####
## Monte Carlo resampling of the tree-ring dates
cutting.mc <- do.call(cbind,mclapply(unique(sites.montecarlo[['Run']]),function(i){
  sites.montecarlo[C_Level %in% 2:3 & Run==i,my.counter(Outer_Date_AD,1:2000)]
}, mc.cores=8))

non_cutting.mc <- do.call(cbind,mclapply(unique(sites.montecarlo[['Run']]),function(i){
  sites.montecarlo[C_Level %in% 1 & Run==i,my.counter(Outer_Date_AD,1:2000)]
}, mc.cores=8))

cutting.mean <- filter(apply(cutting.mc,1,mean, na.rm=T),filter=dist)
cutting.ci <- filter(apply(cutting.mc,1,se)*1.96,filter=dist)

non_cutting.mean <- filter(apply(non_cutting.mc,1,mean, na.rm=T),filter=dist)
non_cutting.ci <- filter(apply(non_cutting.mc,1,se)*1.96,filter=dist)

remove.na = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)),i:=0,with=FALSE]
}


cells.cutting <-  sites.montecarlo[C_Level %in% c(2,3),.N,by=.(Outer_Date_AD,CELL,Run)]
cells.noncutting <-  sites.montecarlo[C_Level %in% c(1),.N,by=.(Outer_Date_AD,CELL,Run)]
cells <- merge(cells.cutting,cells.noncutting,by=c("Outer_Date_AD","CELL","Run"), all=T)
data.table::setnames(cells,c("Outer_Date_AD","CELL","Run","N.cutting","N.non_cutting"))
remove.na(cells)

cells[,TOTAL_N:=N.cutting + N.non_cutting]
cells[,TOTAL_P:=N.cutting/TOTAL_N]

mean.cells <- cells[,list(.N,mean(TOTAL_P)),by=.(Outer_Date_AD,Run)][,.(sum(N)/999,mean(V2)),by=Outer_Date_AD]
setnames(mean.cells,c("Outer_Date_AD","COUNT","MEAN_P"))

cells.stack <- do.call(rbind,
                       lapply(rev(seq(0,1,by=0.01)), function(prob){
                         temp <- vector('numeric',2000)
                         temp[mean.cells[MEAN_P >= prob & MEAN_P < (prob+0.01)][["Outer_Date_AD"]]] <- mean.cells[MEAN_P >= prob & MEAN_P < (prob+0.01)][["COUNT"]]
                         return(temp)
                       }))

cells.stack.smooth <- t(apply(cells.stack,1,function(x){filter(x=x, filter=dist)}))
cells.stack.smooth[is.na(cells.stack.smooth)] <- 0
cells.stack <- apply(cells.stack,2,cumsum)
cells.stack.smooth <- apply(cells.stack.smooth,2,cumsum)


level.colors.ramp <- colorRampPalette(level.colors)(101)
pueblo.subperiod.breaks <- readRDS("../OUTPUT/pueblo.subperiod.breaks.Rds")

phi <- (1+sqrt(5))/2
mai <- c(0.25,0.25,0,0)
nplots <- 2
fig.width <- 4.6
fig.height <- nplots*fig.width/(phi^2)
between <- 0.2
margins <- 0.5
plot.height <- (fig.height - (margins*1.5) - (between*(nplots-1)))/nplots
legend.cex <- 0.75

quartz(file="../FIGURES/FIG_3.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=8, dpi=600)
par(ljoin=0)
this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=F)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,2.5), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
polygon(x=c(500:1400,1400:500),y=c(non_cutting.mean[500:1400]+non_cutting.ci[500:1400], rev(non_cutting.mean[500:1400]-non_cutting.ci[500:1400])),col=paste0(level.colors[3],"80"), border=NA, ljoin=0)
lines(y=non_cutting.mean[500:1400], x=500:1400, col=level.colors[3], lwd=0.75)
polygon(x=c(500:1400,1400:500),y=c(cutting.mean[500:1400]+cutting.ci[500:1400], rev(cutting.mean[500:1400]-cutting.ci[500:1400])),col=paste0(level.colors[1],"80"), border=NA, ljoin=0)
lines(y=cutting.mean[500:1400], x=500:1400, col=level.colors[1], lwd=0.75)

this.plot <- 2
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,3.5), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
for(i in nrow(cells.stack.smooth):1){
  polygon(x=c(500:1400,1400:500),y=c(cells.stack.smooth[i,500:1400],rep(0,length(cells.stack.smooth[i,500:1400]))),col=level.colors.ramp[i], border=NA)
}

par(mai=c(margins,margins*1.25,margins*0.5,margins*0.5), xpd=F, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,1),xlim=c(500,1400), axes=FALSE, main='')
abline(v=pueblo.subperiod.breaks[c(3,5,7,9)], lwd=1.5,lty=2)
abline(v=pueblo.subperiod.breaks[c(2,4,6,8)], lwd=1,lty=3)
text(x=diff(pueblo.subperiod.breaks[c(1,3,5,7,9,10)])/2 + pueblo.subperiod.breaks[c(1,3,5,7,9)],y=par('usr')[4],labels=c("BM III","P I","P II","P III","P IV"),adj=c(0.5,-1), xpd=T)

## LEGENDS and AXES
this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,2.5), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
par(lend='butt')
legend("topleft", col=paste0(level.colors[c(1,3)],80), lwd=10, legend=c("Cutting and near-cutting","Non-cutting"), text.col='#FFFFFF00', bg='#FFFFFFBF', box.col=NA, border=NA)
legend("topleft", col=level.colors[c(1,3)], lwd=1, legend=c("Cutting and near-cutting","Non-cutting"), bg='#FFFFFF00', box.col=NA, border=NA)
par(lend='round')
axis(1, at=seq(500,1400,100), labels=NA, cex.axis=0.8)
axis(2, at=seq(0,2.5,.5), cex.axis=0.8)
mtext("Number of dates", side=2, line=2.25)
text(x=375,y=2.5, labels=LETTERS[this.plot], font=2, cex=1.75)

this.plot <- 2
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins*1.25, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,3.5), xlim=c(500,1400), axes=FALSE, main='', xpd=T)

color.breaks <- seq(from=525,to=875,by=(875-525)/length(level.colors.ramp))
rect(xleft=500, ybottom=1.75, xright=900, ytop=3.5, col='#FFFFFFBF', border=NA)
rect(xleft=color.breaks[1:(length(color.breaks)-1)], ybottom=2.85, xright=color.breaks[2:length(color.breaks)], ytop=3.25, col=rev(level.colors.ramp), border=NA)
rect(xleft=color.breaks[1], ybottom=2.85, xright=color.breaks[length(color.breaks)], ytop=3.25, col=NA, border="black", lwd=0.5)
text(x=seq(525,875,length.out=5),y=2.85,labels=paste0(seq(0,100,25),"%"),adj=c(0.5,2),xpd=T, cex=0.8)
text(x=700, y=2.85, labels="Percent cutting and near-cutting", adj=c(0.5,3.5), xpd=T, cex=1)

axis(1, at=seq(500,1400,100), labels=seq(500,1400,100), cex.axis=0.8)
axis(2, at=seq(0,3.5,0.5), cex.axis=0.8)
mtext("Number of 800-m cells", side=2, line=2.25)
mtext("Year AD", side=1, line=2.25)
text(x=375,y=3.5, labels=LETTERS[this.plot], font=2, cex=1.75)

dev.off()
distill("../FIGURES/FIG_3.pdf")


##### FIG_4.pdf #####
## Prediction error of PaleoCAR reconstruction
SWUS.extent <- extent(-113,-105,32,38)
SWUS.poly <- polygon_from_extent(SWUS.extent,proj4string="+proj=longlat +datum=WGS84")
if(!file.exists("../OUTPUT/SWUS.NED.small.tif")){
  SWUS.NED <- get_ned(template=SWUS.poly,label="SWUS",res='1',raw.dir='../DATA/NED/', extraction.dir='../DATA/NED/EXTRACTIONS/')
  SWUS.NED.small <- aggregate(SWUS.NED, 10)
  raster::writeRaster(SWUS.NED.small, "../OUTPUT/SWUS.NED.small.tif", datatype = "FLT4S", options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite = T, setStatistics = FALSE)
}
SWUS.NED.small <- crop(raster("../OUTPUT/SWUS.NED.small.tif"),SWUS.poly,snap='out')*2
slope <- terrain(SWUS.NED.small, opt='slope')
aspect <- terrain(SWUS.NED.small, opt='aspect')
hill <- hillShade(slope, aspect, 40, 230)

phi <- (1+sqrt(5))/2
plot.ratio <- 1.093363

## Export PNGs of the spatial errors for quick plotting
quartz(file='../FIGURES/PPT.ERROR.SPACE.png', width=6*plot.ratio, height=6, antialias=FALSE, bg="transparent", type='png', family="Gulim", pointsize=8, dpi=300)
par(mai=c(0,0,0,0))
plot(1, type='n', xlab="", ylab="", xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
colors <- paste0(colorRampPalette(brewer.pal(9,"Greens"))(40),c("00","0D","1A","26","33","40","4D","59","66","73","80")[11])
color.breaks <- seq(from=0,to=40,by=1)
plot(hill, maxpixels=ncell(hill), col=grey(60:100/100), useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(PPT.ERROR.SPACE, zlim=c(0,40), maxpixels=ncell(PPT.ERROR.SPACE), breaks=color.breaks, col=colors, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
dev.off()

quartz(file='../FIGURES/GDD.ERROR.SPACE.png', width=6*plot.ratio, height=6, antialias=FALSE, bg="transparent", type='png', family="Gulim", pointsize=8, dpi=300)
par(mai=c(0,0,0,0))
plot(1, type='n', xlab="", ylab="", xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
colors <- paste0(colorRampPalette(brewer.pal(9,"Greens"))(20),c("00","0D","1A","26","33","40","4D","59","66","73","80")[11])
color.breaks <- seq(from=5,to=25,by=1)
plot(hill, maxpixels=ncell(hill), col=grey(60:100/100), useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(GDD.ERROR.SPACE, zlim=c(5,25), maxpixels=ncell(GDD.ERROR.SPACE), breaks=color.breaks, col=colors, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
dev.off()

fig.width <- 6.5
between <- 0.3
margins <- 0.5

plot.width <- (fig.width - (between*2) - margins*0.25)/2
map.height <- plot.width/plot.ratio
plot.height <- (plot.width/phi)*0.75
fig.height <- plot.height + map.height + between*1.75

legend.cex <- 0.75



quartz(file="../FIGURES/FIG_4.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=8, dpi=600)
par(bg='white',fg='black',col.lab='black', col.main='black', col.axis='black', font=2, lend='round',ljoin='round')

par(mai=c(0,0,0,0), lend='round', ljoin='round', xpd=T)
plot(1, type='n', xlab="", ylab="",xaxs='i',yaxs='i', xlim=c(0,fig.width), ylim=c(0,fig.height), axes=FALSE, main='')
text(x=c(between/2,between/2,between+plot.width+between/2,between+plot.width+between/2),y=c(fig.height-between/2, fig.height - map.height - between*1.75, fig.height-between/2, fig.height - map.height - between*1.75), labels=LETTERS[1:4], font=2, cex=1.75)

par(mai=c(plot.height + between*1.75, between, 0, plot.width+between+margins*0.25), xpd=F, new=T)
plot(1, type='n', xlab="", ylab="",xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
rasterImage(readPNG('../FIGURES/PPT.ERROR.SPACE.png'), xleft=xmin(SWUS.extent), xright=xmax(SWUS.extent), ybottom=ymin(SWUS.extent), ytop=ymax(SWUS.extent), interpolate=F)
par(xpd=T)
colors <- paste0(colorRampPalette(brewer.pal(9,"Greens"))(40),c("00","0D","1A","26","33","40","4D","59","66","73","80")[11])
color.breaks <- seq(from=xmin(SWUS.extent),to=xmax(SWUS.extent),length.out=41)
rect(xleft=xmin(SWUS.extent), ybottom=31.7, xright=xmax(SWUS.extent), ytop=31.9, col='#FFFFFFBF', border=NA)
rect(xleft=color.breaks[1:(length(color.breaks)-1)], ybottom=31.7, xright=color.breaks[2:length(color.breaks)], ytop=31.9, col=colors, border=NA)
rect(xleft=color.breaks[1], ybottom=31.7, xright=color.breaks[length(color.breaks)], ytop=31.9, col=NA, border="gray80", lwd=0.5)
text(x=seq(xmin(SWUS.extent),xmax(SWUS.extent),length.out=9),y=31.7,labels=seq(0,40,5),adj=c(0.5,2),xpd=T, cex=0.8)
text(x=mean(c(xmin(SWUS.extent),xmax(SWUS.extent))), y=31.7, labels="Mean Error (mm)", adj=c(0.5,3.0), xpd=T, cex=1)

par(mai=c(margins, between + margins*0.75, map.height + between*1.75, plot.width + between + margins*0.25), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(8.4,9.4), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
lines(y=PPT.ERROR.TIME, x=500:1400)

axis(1, at=seq(500,1400,100), labels=seq(500,1400,100), cex.axis=0.8)
axis(2, at=seq(8.4,9.4,0.2), labels=c(8.4,"","","","",9.4), cex.axis=0.8)
mtext("Mean Error (mm)", side=2, line=2.25)
mtext("Year AD", side=1, line=2.25)


par(mai=c(plot.height + between*1.75, plot.width + (between*2),0,margins*0.25), xpd=F, new=T)
plot(1, type='n', xlab="", ylab="",xlim=c(xmin(SWUS.extent),xmax(SWUS.extent)), ylim=c(ymin(SWUS.extent),ymax(SWUS.extent)), xaxs="i", yaxs="i", axes=FALSE, main='')
rasterImage(readPNG('../FIGURES/GDD.ERROR.SPACE.png'), xleft=xmin(SWUS.extent), xright=xmax(SWUS.extent), ybottom=ymin(SWUS.extent), ytop=ymax(SWUS.extent), interpolate=F)
par(xpd=T)
colors <- paste0(colorRampPalette(brewer.pal(9,"Greens"))(20),c("00","0D","1A","26","33","40","4D","59","66","73","80")[11])
color.breaks <- seq(from=xmin(SWUS.extent),to=xmax(SWUS.extent),length.out=20)
rect(xleft=xmin(SWUS.extent), ybottom=31.7, xright=xmax(SWUS.extent), ytop=31.9, col='#FFFFFFBF', border=NA)
rect(xleft=color.breaks[1:(length(color.breaks)-1)], ybottom=31.7, xright=color.breaks[2:length(color.breaks)], ytop=31.9, col=colors, border=NA)
rect(xleft=color.breaks[1], ybottom=31.7, xright=color.breaks[length(color.breaks)], ytop=31.9, col=NA, border="gray80", lwd=0.5)
text(x=seq(xmin(SWUS.extent),xmax(SWUS.extent),length.out=5),y=31.7,labels=seq(5,25,5),adj=c(0.5,2),xpd=T, cex=0.8)
text(x=mean(c(xmin(SWUS.extent),xmax(SWUS.extent))), y=31.7, labels="Mean Error (GDD)", adj=c(0.5,3.0), xpd=T, cex=1)

par(mai=c(margins, plot.width + between*2 + margins*0.75, map.height + between*1.75, margins*0.25), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(12,13), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
lines(y=GDD.ERROR.TIME, x=500:1400)

axis(1, at=seq(500,1400,100), labels=seq(500,1400,100), cex.axis=0.8)
axis(2, at=seq(12,13,0.2), labels=c(12,"","","","",13), cex.axis=0.8)
mtext("Mean Error (GDD)", side=2, line=2.25)
mtext("Year AD", side=1, line=2.25)

dev.off()
distill("../FIGURES/FIG_4.pdf")


##### FIG_6.pdf #####

phi <- (1+sqrt(5))/2
mai <- c(0.25,0.25,0,0)
nplots <- 1
fig.width <- 4.6
fig.height <- nplots*fig.width/phi
between <- 0.2
margins <- 0.5
plot.height <- (fig.height - (margins*1.5) - (between*(nplots-1)))/nplots
legend.cex <- 0.75

quartz(file="../FIGURES/FIG_6.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=8, dpi=600)

this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=F)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,4), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
mids <- pueblo.subperiod.breaks[-1] - diff(pueblo.subperiod.breaks)/2
spread <- 15
i <- 1

for(summaries in total.summaries.four){
  segments(x0=c(mids[[i]]),
           x1=c(mids[[i]]),
           y0=c(summaries[["mean"]]-summaries[["sd"]]),
           y1=c(summaries[["mean"]]+summaries[["sd"]]),
           col='black',
           lwd=2)
  points(y=c(summaries[["mean"]]),
         x=c(mids[[i]]),
         col='black',
         cex=2,
         pch=19)
  i <- i+1
}

i <- 1
for(summaries in period.summaries.four){
  segments(x0=c(mids[[i]]-spread,mids[[i]]+spread),
           x1=c(mids[[i]]-spread,mids[[i]]+spread),
           y0=c(summaries$all[["mean"]]-summaries$all[["sd"]],summaries$local[["mean"]]-summaries$local[["sd"]]),
           y1=c(summaries$all[["mean"]]+summaries$all[["sd"]],summaries$local[["mean"]]+summaries$local[["sd"]]),
           col=level.colors[c(1,3)],
           lwd=2)
  points(y=c(summaries$all[["mean"]],summaries$local[["mean"]]),
         x=c(mids[[i]]-spread,mids[[i]]+spread),
         col=level.colors[c(1,3)],
         cex=2,
         pch=19)
  i <- i+1
}



par(mai=c(margins,margins,margins*0.5,margins*0.5), xpd=F, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,4),xlim=c(500,1400), axes=FALSE, main='')
abline(v=pueblo.subperiod.breaks[c(3,5,7,9)], lwd=1.5,lty=2)
abline(v=pueblo.subperiod.breaks[c(2,4,6,8)], lwd=1,lty=3)
text(x=diff(pueblo.subperiod.breaks[c(1,3,5,7,9,10)])/2 + pueblo.subperiod.breaks[c(1,3,5,7,9)],y=par('usr')[4],labels=c("BM III","P I","P II","P III","P IV"),adj=c(0.5,-1), xpd=T)

## LEGENDS and AXES
this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=T, new=T)
plot(1, type='n', xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0,4), xlim=c(500,1400), axes=FALSE, main='', xpd=T)
legend("bottomleft", fill=c(level.colors[c(1,3)],'black'), legend=c("Cells with dates","Local cells","All cells"), bg='#FFFFFFBF', box.col=NA, border=NA)
axis(1, at=seq(500,1400,100), labels=seq(500,1400,100), cex.axis=0.8)
axis(2, at=seq(0,4,1), cex.axis=0.8)
mtext(bquote(.("Count of past four years in niche (") ~ bar(x) %+-% σ ~ .(")")), side=2, line=2.09)
mtext("Year AD", side=1, line=2.25)

dev.off()
distill("../FIGURES/FIG_6.pdf")


##### FIG_7.pdf #####
## Histogram of site-level counts of tree-ring dates

phi <- (1+sqrt(5))/2
mai <- c(0.25,0.25,0,0)
nplots <- 1
fig.width <- 4.6
fig.height <- nplots*fig.width/phi
between <- 0.2
margins <- 0.5
plot.height <- (fig.height - (margins*1.5) - (between*(nplots-1)))/nplots
legend.cex <- 0.75

quartz(file="../FIGURES/FIG_7.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=8, dpi=600)

this.plot <- 1
par(mai=c(margins + (plot.height * (nplots-this.plot)) + (between * (nplots-this.plot)), margins, (margins * 0.5) + (plot.height * (this.plot-1)) + (between * (this.plot-1)), margins*0.5), xpd=F)
this.hist <- hist(table(sites[Outer_Date_AD %in% 500:1400][['Seq_Number']]), breaks=100, plot=F)
this.hist$count[this.hist$count==0] <- NA
plot(y=this.hist$count,x=this.hist$mids, type='h', log='y', lwd=5, lend=1, xaxs="i",yaxs="i", xlab="", ylab="", ylim=c(0.9,600), xlim=c(0,1000), axes=FALSE, main='', xpd=T)

axis(1, at=seq(0,1000,200), labels=seq(0,1000,200), cex.axis=0.8)
axis(2, at=c(1,10,100,200,600),labels=c(1,10,100,200,600), cex.axis=0.8)
mtext("Number of sites", side=2, line=2.09)
mtext("Number of dates", side=1, line=2.25)

dev.off()
distill("../FIGURES/FIG_7.pdf")


 ##### MOVIE_S1.pdf #####
## This section contains code that creates a pseudo-3D animation
## of the joint-gensis of the maize dry-farming niche and the tree-ring dates.
# First, create a directory for the animation frames.
dir.create("../FIGURES/MOVIE_S1/",showWarnings = F, recursive = T)

## Create the legend
open3d()
par3d(
  "windowRect"= c(0,0,345,50),
  "zoom"= 0.15,
  "FOV"= 0,
  "userMatrix"=matrix(c(1,0,0,0,0,0.35,0.95,0,0,-0.95,0.35,0,0,0,0,1), ncol=4, byrow=T),
  "viewport"=c(0,0,345,50)
)
spheres3d(x=seq(-113,-105,length.out=5),y=c(2,2,2,2,2),z=0,radius=seq(0.02,0.1,by=0.02)*(1600/345),col=colorRampPalette(rev(level.colors))(5))
rgl.snapshot('../FIGURES/legend_color_sizes.png')
rgl.close()
legend_color_sizes.png <- readPNG('../FIGURES/legend_color_sizes.png')


dist.2 <- dnorm(seq(-10,10,1), sd=2)
year.range <- 501:1400
i <- 0
for(year in year.range){
  i <- i+1
  if(file.exists(paste('../FIGURES/MOVIE_S1/image',formatC(i, width = 4, format = "d", flag = "0"),'.png',sep=''))) next
  smooth.rast <- focal(NICHE.SMOOTH[[year]],w=crossprod(t(dist.2),dist.2), fun=sum, pad=T, na.rm=T)

  plot3D(smooth.rast,maxpixels=ncell(smooth.rast), col=colorRampPalette(brewer.pal(9,"Greys")[1:7])(1000), zlim=c(0,4), adjust=F, zfac=0.05, lit=T, specular="black", useLegend=F)
  par3d("ignoreExtent"=T)
  this.spheres <- all.cells[Outer_Date_AD==year]
  if(nrow(this.spheres)!=0){
    setorder(this.spheres,-SIZE,-COLOR)
    this.spheres[,x:=xFromCell(FOUR_CORNERS.rast,this.spheres[["CELL"]])]
    this.spheres[,y:=yFromCell(FOUR_CORNERS.rast,this.spheres[["CELL"]])]
    this.spheres[,z:=(extract(smooth.rast,cbind(x,y))*0.05)+SIZE]
    spheres3d(x=this.spheres[["x"]],y=this.spheres[["y"]],z=this.spheres[["z"]],radius=this.spheres[["SIZE"]],col=this.spheres[["COLOR"]])
  }
  
  par3d(
    "windowRect"= c(0,0,1600,1200),
    "zoom"= 0.65,
    "FOV"= 30,
    "userMatrix"=matrix(c(1,0,0,0,0,0.8,0.65,0,0,-0.65,0.8,0,0,0,0,1), ncol=4, byrow=T),
    "viewport"=c(0,0,1600,1200)
  )
  
  rgl.snapshot('../FIGURES/MOVIE_S1/temp.png')
  rgl.close()
  
  temp.png <- readPNG('../FIGURES/MOVIE_S1/temp.png')
  
  quartz(file=paste('../FIGURES/MOVIE_S1/image',formatC(i, width = 4, format = "d", flag = "0"),'.png',sep=''), width=8, height=6, bg="white", type='png', family="Helvetica Bold", pointsize=8, dpi=200)
  
  par(mai=c(0,0,0,0))
  plot(1, type='n', xlab="", ylab="", xlim=c(0,1600),ylim=c(0,1200), xaxs="i", yaxs="i", axes=FALSE, main='')
  rasterImage(temp.png[201:1000,1:1600,], 0,400,1600,1200, interpolate=F)
  
  text(25,25, labels=paste0("AD ",year), adj=c(0,0), font=2, cex=4)
  
  
  # PLOT
  par(mai=c(0.5,2.5,4.25,0.225), new=T, xpd=F)
  plot(1, type='n', xlab="", ylab="",xaxs='i',yaxs='i', xlim=c(500,1400), ylim=c(0,1), axes=FALSE, main='')
  abline(v=year, col='gray50')
  rasterImage(plot.png, 500,0,1400,1, interpolate=F)
  axis(1, at=seq(500,1400,100), labels=seq(500,1400,100), cex.axis=1, font=2)
  axis(2, at=seq(0,1,0.2), labels=seq(0,100,20), las=1,cex.axis=1,font=2)

  mtext("Percent in niche", side=2, line=2.5, cex=1.25, font=2)
  mtext("Year AD", side=1, line=2.5, cex=1.25, font=2)
  
  abline(v=pueblo.subperiod.breaks[c(3,5,7,9)], lwd=1.5,lty=2)
  abline(v=pueblo.subperiod.breaks[c(2,4,6,8)], lwd=1,lty=3)
  text(x=diff(pueblo.subperiod.breaks[c(1,3,5,7,9,10)])/2 + pueblo.subperiod.breaks[c(1,3,5,7,9)],y=par('usr')[4],labels=c("BM III","P I","P II","P III","P IV"),adj=c(0.5,-1), xpd=T, font=2)
  
  
  # LEGENDS
  par(mai=c(0.6,0.125,4,6.15), new=T, xpd=T)
  plot(1, type='n', xlab="", ylab="",xaxs='i',yaxs='i', xlim=c(0,345), ylim=c(0,280), axes=FALSE, main='')
  rasterImage(legend_color_sizes.png, 0,175,345,225, interpolate=F)
  
  text(x=172.5, y=225, labels="Size: Number of dates (logged)", adj=c(0.5,-1.5), font=2)
  text(x=c(11.5, 88, 165.5, 242.5, 320),y=225,labels=c(1,4,12,41,141), adj=c(0.5,0),cex=0.75, font=2)
  text(x=c(11.5, 88, 165.5, 242.5, 320),y=175,labels=c("0%","25%","50%","75%","100%"), adj=c(0.5,1),cex=0.75, font=2)
  text(x=172.5, y=175, labels="Color: Percent cutting dates", adj=c(0.5,2.5), font=2)
  
  colors <- colorRampPalette(brewer.pal(9,"Greys")[2:8])(100)
  color.breaks <- seq(from=10,to=335,by=345/length(colors))
  text(x=172.5,y=65,labels="Elevation & grayscale:",adj=c(0.5,-1),xpd=T, font=2)
  rect(xleft=color.breaks[1:(length(color.breaks)-1)], ybottom=40, xright=color.breaks[2:length(color.breaks)], ytop=65, col=colors, border=NA)
  rect(xleft=color.breaks[1], ybottom=40, xright=color.breaks[length(color.breaks)], ytop=65, col=NA, border="black", lwd=0.5)
  text(x=seq(10,335,length.out=5),y=40,labels=seq(0,4,1),adj=c(0.5,2),xpd=T, cex=0.75, font=2)
  text(x=172.5, y=40, labels="Count of past four years in niche", adj=c(0.5,3), xpd=T,font=2)
  
  dev.off()
}
unlink('../FIGURES/MOVIE_S1/temp.png')

## Compile the images into video.
## Create videos of different sizes and compression levels.
time <- 120 # Length in seconds of video 
fps <- length(year.range)/time
system(paste0("ffmpeg -r ",fps," -i ../FIGURES/MOVIE_S1/image%04d.png -s:v ",1600,"x",1200," -c:v libx264 -profile:v high444 -crf 30 -pix_fmt yuv420p ../FIGURES/MOVIE_S1.mov -y"))
system(paste0("ffmpeg -r ",fps," -i ../FIGURES/MOVIE_S1/image%04d.png -s:v ",800,"x",600," -c:v libx264 -profile:v high444 -crf 30 -pix_fmt yuv420p ../FIGURES/MOVIE_S1_SMALL.mov -y"))
system(paste0("ffmpeg -r ",fps," -i ../FIGURES/MOVIE_S1/image%04d.png -s:v ",640,"x",480," -c:v libx264 -profile:v high444 -crf 30 -pix_fmt yuv420p ../FIGURES/MOVIE_S1_TINY.mov -y"))
