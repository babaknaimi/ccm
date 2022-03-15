# Authors: Shirin Taheri, taheri.shi@gmail.com; Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Dec. 2021
# Version 1.6
# Licence GPL v3
#--------



.getProj <- function(x) {
  if (inherits(x,'Raster')) {
    if (!is.na(projection(x))) strsplit(strsplit(projection(x),'\\+proj=')[[1]][2],' ')[[1]][1]
    else {
      if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
      else 'projected'
    }
  } else {
    if (!is.na(crs(x))) strsplit(strsplit(crs(x,proj=TRUE),'\\+proj=')[[1]][2],' ')[[1]][1]
    else {
      if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
      else 'projected'
    }
  }
  
}
#=-===============
.is.projected <- function(x) {
  if (inherits(x,'SpatRaster')) {
    e <- as.vector(ext(x))
  } else e <- as.vector(extent(x))
  
  !all(e >= -180 & e <= 180)
}
#-----

.rad <- function (degree) {
  (degree * pi) / 180
}
#---
.deg <-  function (radian) {
  (radian * 180) / pi
}



# following two functions are copied from the package vocc in "https://github.com/cbrown5/vocc"
.mnwm <- function(d1, d2, d3, d4, d5, d6){
  X <- sum(c(d1, d2*2, d3, d4, d5*2, d6), na.rm = T)
  w <- sum(c(1,2,1,1,2,1) * is.finite(c(d1, d2, d3, d4, d5, d6)))
  return(X/w)
}
#-----
.ang <- function(dx, dy){
  ifelse(dy < 0, 180 + .deg(atan(dx/dy)),
         ifelse(dx < 0, 360 + .deg(atan(dx /dy )), .deg(atan(dx/dy))))
}
#--------
######################



# this function stretch each variable between 0 and 1 (in T1, also in T2 but based on parameters of T1), 
# and then the multiple variable are averaged to have a single variable for each time.
.getScaledMultiVariateIntoOne <- function(x,...,t1,t2) {
  lst <- list(x,...)
  xt1 <- xt2 <- list()
  
  if (length(t1) > 1 | length(t2) > 1) .multi <- TRUE
  else .multi <- FALSE
  
  if (inherits(x,'Raster')) .raster <- TRUE
  else .raster <- FALSE
  
  
  
  for (i in 1:length(lst)) {
    if (.multi) {
      if (.raster) {
        xt1[[i]] <- mean(rast(lst[[i]])[[t1]])
        xt2[[i]] <- mean(rast(lst[[i]])[[t2]])
      } else {
        xt1[[i]] <- mean(lst[[i]][[t1]])
        xt2[[i]] <- mean(lst[[i]][[t2]])
      }
    } else {
      if (.raster) {
        xt1[[i]] <- rast(lst[[i]])[[t1]]
        xt2[[i]] <- rast(lst[[i]])[[t2]]
      } else {
        xt1[[i]] <- lst[[i]][[t1]]
        xt2[[i]] <- lst[[i]][[t2]]
      }
    }
  }
  #----
  
  .scaleParam <- vector(mode='list',length=length(xt1))
  names(.scaleParam) <- sapply(xt1,names)
  
  
  for (i in 1:length(xt1)) {                                                                                                                                                                                                                                                                                      
    .scaleParam[[i]] <- minmax(xt1[[i]])[,1]
    xt1[[i]] <- (xt1[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1] )
    xt2[[i]] <- (xt2[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1])
  }
  #-----
  x1 <- xt1[[1]]
  x2 <- xt2[[1]]
  
  if (length(xt1) > 1) {
    for (i in 2:length(xt1)) {
      x1 <- x1 + xt1[[i]]
      x2 <- x2 + xt2[[i]]
    }
    x1 <- x1 / length(xt1)
    x2 <- x2 / length(xt1)
  }
  list(t1=x1,t2=x2)
}

#----------
.getLongLat <- function(x,crs) {
  # x is data.frame
  if (ncol(x) > 2) {
    if (all(c('x','y') %in% colnames(x))) x <- x[,c('x','y')]
    else {
      x <- x[,1:2]
      colnames(x) <- c('x','y')
      warning('".getLongLat function:" the first two columns in the input data.frame is assumed as the coordinate columns!')
    }
  } else if (ncol(x) < 2) stop('".getLongLat function:" input x should be a data.frame with 2 columns, i.e., x and y coordinates!')
  
  
  geom(project(vect(x,geom=colnames(x),crs=crs),"epsg:4326"))[,c('x','y')]
}
#-----
