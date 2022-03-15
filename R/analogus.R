# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Dec. 2021
# Version 2.1
# Licence GPL v3
#--------



.analogusClimate <- function(k1,k2) {
  r <- raster(k1)
  u <- unique(k1[])
  u <- u[!is.na(u)]
  for (v in u) {
    w1 <- length(which(k1[] == v))
    w2 <- length(which(k2[] == v)) 
    w <- which(k1[] == v)
    r[w] <- ((w2-w1)/w1) * 100
  }
  r
}
#--------------
.disAnalogus <- function(k1,k2,.lon=FALSE) {
  
  if (.getProj(k1) == 'longlat') .lon <- TRUE
  
  w <- which(!is.na(k1[]))
  r <- raster(k1)
  .getProj(k1)
  
  
  for (.c in w) {
    .xy <- xyFromCell(k1,.c)
    .xy2 <- xyFromCell(k1,which(k1[] == k1[.c]))
    .d1 <- pointDistance(.xy,.xy2,lonlat = .lon)
    .xy2 <- xyFromCell(k2,which(k2[] == k1[.c]))
    .d2 <- pointDistance(.xy,.xy2,lonlat = .lon)
    q <- quantile(c(.d1,.d2),0.1)
    v <- median(.d2[.d2 < q],na.rm=TRUE) - median(.d1[.d1 < q],na.rm=TRUE)
    r[.c] <- v
  }
  
  r
}