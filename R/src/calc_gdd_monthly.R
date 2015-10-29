calc_gdd_monthly <- function(tmin_brick, tmax_brick, t.base, t.cap=NULL, multiplier=1, to_fahrenheit=T, output.dir='./gdd/'){
  if(nlayers(tmin_brick)!=nlayers(tmax_brick)){
    stop("tmin and tmax bricks must have same number of layers!")
  }
  
  dir.create(output.dir, recursive=T, showWarnings=F)
  
  t.base <- t.base*multiplier
  if(!is.null(t.cap)){
    t.cap <- t.cap*multiplier
  }
  
  files <- names(tmin_brick)
  
  cl <- parallel::makeCluster(8)
  
  junk <- parallel::parLapply(cl, files, function(tmin_brick, tmax_brick, output.dir, t.base, t.cap, multiplier,to_fahrenheit, file){
    if(file.exists(paste0(output.dir,file,".tif"))) return()
    
    tmin <- raster::subset(tmin_brick,subset=file)
    tmax <- raster::subset(tmax_brick,subset=file)
    
    # Floor tmax and tmin at Tbase
    tmin <- raster::calc(tmin,function(x) { x[x<t.base] <- t.base; return(x) })
    tmax <- raster::calc(tmax,function(x) { x[x<t.base] <- t.base; return(x) })
    
    # Cap tmax and tmin at Tut
    if(!is.null(t.cap)){
      tmin <- raster::calc(tmin,function(x) { x[x>t.cap] <- t.cap; return(x) })
      tmax <- raster::calc(tmax,function(x) { x[x>t.cap] <- t.cap; return(x) })
    }
    
    GDD <- ((tmin+tmax)/2)-t.base
    
    # Multiply by days per month, and convert to Fahrenheit GDD
    GDD <- GDD * Hmisc::monthDays(as.Date(paste0(file,"D01"),format="Y%YM%mD%d")) / multiplier
    
    if(to_fahrenheit){
      GDD <- GDD * 1.8
    }
    
    GDD <- round(GDD)
    
    raster::writeRaster(GDD,paste0(output.dir,file,".tif"), datatype="INT2U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
    return()
  }, tmin_brick=tmin_brick, tmax_brick=tmax_brick,output.dir=output.dir, t.base=t.base, t.cap=t.cap, multiplier=multiplier, to_fahrenheit=to_fahrenheit)
  
  parallel::stopCluster(cl)
  
  return(stack(list.files(output.dir, full.names=T), quick=T))
}