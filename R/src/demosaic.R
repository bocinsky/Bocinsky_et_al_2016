demosaic <- function(raster_brick, corners, out_dir, corner_loc="SW", datatype="INT2S"){
  dir.create(out_dir, showWarnings=F)
  
  corners <- lapply(1:nrow(corners),function(i){as.numeric(corners[i,])})
  
  cl <- parallel::makeCluster(8)
  parallel::clusterEvalQ(cl, {library(raster)})
  junk <- parallel::parLapply(cl, corners, function(raster_brick, out_dir, datatype, corner){
    if(file.exists(paste0(out_dir,0-corner[1],"W",corner[2],"N",".tif"))) return()
    out <- raster::crop(raster_brick,raster::extent(corner[1],corner[1]+1,corner[2],corner[2]+1), filename=paste0(out_dir,0-corner[1],"W",corner[2],"N",".tif"), datatype=datatype, options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
    return()
  }, raster_brick=raster_brick, out_dir=out_dir, datatype=datatype)
  parallel::stopCluster(cl)
  
  if(class(raster_brick)=="RasterLayer"){
    return(lapply(list.files(out_dir, full.names=T), raster::raster))
  }else{
    return(lapply(list.files(out_dir, full.names=T), raster::stack, quick=T))
  }
}