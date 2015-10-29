distill <- function(file){
  system(paste0("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH -dCompatibilityLevel=1.4 -sColorConversionStrategy=/CMYK -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=false -dAutoRotatePages=/None -sOutputFile=./temp.pdf ",file))
  system(paste0("rm ",file))
  system(paste0("mv ./temp.pdf ",file))
  system(paste0("convert -density 600 ",file," ",gsub(".pdf",".png",file)))
}