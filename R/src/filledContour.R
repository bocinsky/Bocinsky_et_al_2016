filledContourAdd <- function (x, y = 1, maxpixels = 1e+05, levels, col) 
{
  if (nlayers(x) > 1) {
    y <- min(max(1, y), nlayers(x))
    x <- raster(x, y)
  }
  x <- sampleRegular(x, maxpixels, asRaster = TRUE, useGDAL = TRUE)
  X <- xFromCol(x, 1:ncol(x))
  Y <- yFromRow(x, nrow(x):1)
  Z <- t(matrix(getValues(x), ncol = x@ncols, byrow = TRUE)[nrow(x):1, 
                                                            ])
  .filled.contour(x = X, y = Y, z = Z, levels=levels, col=col)
}