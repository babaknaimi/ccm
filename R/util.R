# Authors: Shirin Taheri, taheri.shi@gmail.com; Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Sep. 2021
# Version 1.2
# Licence GPL v3
#--------



.getProj <- function(x) {
  if (!is.na(projection(x))) strsplit(strsplit(projection(x),'\\+proj=')[[1]][2],' ')[[1]][1]
  else {
    if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
    else 'projected'
  }
}
#=-===============
.is.projected <- function(x) {
  e <- as.vector(extent(x))
  all(e > -180 & e < 180)
}