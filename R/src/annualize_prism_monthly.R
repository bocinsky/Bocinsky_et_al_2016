annualize_prism_monthly <- function(prism.brick, months=c(1:12), fun, out_dir){
  dir.create(out_dir, showWarnings=F)
  brick.names <- names(prism.brick)
  brick.years <- as.numeric(lapply(strsplit(brick.names,"[A-Z]"),'[[',2))
  brick.months <- as.numeric(lapply(strsplit(brick.names,"[A-Z]"),'[[',3))
  
  # If more months than a year, break
  if(length(months)>12){
    stop("ERROR! Too many months.")
  }
  
  # Process the months to get months from current, previous, and future years
  previous.year <- 12+months[months<1]
  current.year <- months[months>=1 & months<=12]
  next.year <- months[months>12]-12
  no.year <- (1:12)[!(1:12 %in% c(previous.year,current.year,next.year))]
  
  brick.years[brick.months %in% previous.year] <- brick.years[brick.months %in% previous.year]+1
  brick.years[brick.months %in% next.year] <- brick.years[brick.months %in% next.year]-1
  brick.years[brick.months %in% no.year] <- 0
  
  prism.brick <- prism.brick[[which(brick.years %in% as.numeric(names(table(brick.years))[table(brick.years)==length(months)]))]]
  brick.years <- brick.years[which(brick.years %in% as.numeric(names(table(brick.years))[table(brick.years)==length(months)]))]
  
  fun.type <- fun
  
  cl <- parallel::makeCluster(8)
  junk <- parallel::parLapply(cl, unique(brick.years), function(out_dir, prism.brick, brick.years, fun.type, year){
    #     cat("\nCalculating year",year)
    if(file.exists(paste0(out_dir,paste0("Y",year),'.tif'))) return()
    
    if(fun.type=="sum"){
      out <- raster::calc(raster::subset(prism.brick, subset=which(brick.years %in% year)), fun=sum, filename=paste0(out_dir,paste0("Y",year),'.tif'), datatype="INT2U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
    }else if(fun.type=="mean"){
      out <- raster::calc(raster::subset(prism.brick, subset=which(brick.years %in% year)), fun=mean, filename=paste0(out_dir,paste0("Y",year),'.tif'), datatype="INT2S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
    }
    return()
  }, out_dir=out_dir, prism.brick=prism.brick, brick.years=brick.years, fun.type=fun.type)
  
  parallel::stopCluster(cl)
  
  return(raster::stack(list.files(out_dir, full.names=T), quick=T))
}
