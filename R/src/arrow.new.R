arrow.new <- function(x,y,height,width,adj){
  arrows(x0=xmin(plot.extent)+0.2*inch.x,y0=ymin(plot.extent)+0.3*inch.y,x1=xmin(plot.extent)+0.2*inch.x,y1=ymin(plot.extent)+0.6*inch.y, length=0.1, lwd=1.5, lend=1)
  text(labels="N",x=xmin(plot.extent)+0.2*inch.x,y=ymin(plot.extent)+0.3*inch.y,adj=c(0.5,0),cex=1.5, font=2)
  
}