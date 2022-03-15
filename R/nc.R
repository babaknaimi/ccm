# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Dec. 2021
# Version 2.1
# Licence GPL v3
#--------
.getRasterList <- function(x1,...,t1,t2) {
  xx <- list(x1,...)
  
  xt1 <- list(xx[[1]][[t1]])
  xt2 <- list(xx[[1]][[t2]])
  
  if (length(xx) > 1) {
    for (i in 2:length(xx)) {
      xt1 <- c(xt1,xx[[i]][[t1]])
      xt2 <- c(xt2,xx[[i]][[t2]])
    }
  }
  list(T1=xt1,T2=xt2)
}

#--------------
.nc <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1@raster,mean,na.rm=TRUE)
    mn2 <- app(x2@raster,mean,na.rm=TRUE)
    v <- app(x1@raster,var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .s <- sqrt(.sed)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
#-----------------
.ncR <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1,mean,na.rm=TRUE)
    mn2 <- app(x2,mean,na.rm=TRUE)
    v <- app(x1@raster,var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .s <- sqrt(.sed)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}

# .nc <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
#   xx <- list()
#   if (!is.null(p)) xx <- c(p,xx)
#   if (!is.null(tmin)) xx <- c(tmin,xx)
#   if (!is.null(tmax)) xx <- c(tmax,xx)
#   if (!is.null(tmean)) xx <- c(tmean,xx)
#   
#   for (i in 1:length(xx)) {
#     x <- xx[[i]]
#     x1 <- x[[t1]]
#     x2 <- x[[t2]]
#     mn1 <- app(x1@raster,mean,na.rm=TRUE)
#     mn2 <- app(x2@raster,mean,na.rm=TRUE)
#     v <- app(x1@raster,var,na.rm=TRUE)
#     
#     w <- which(!is.na(mn2[]))
#     ww <- data.frame(mn1[w],v[w])
#     
#     wm <- sapply(mn2[w][,1],function(x,...) {
#       min(((x - ww[,1]) ^ 2) / ww[,2])
#     },na.rm=TRUE)
#     wm2 <- rast(mn2)
#     wm2[w] <- wm
#     # wm <- app(mn2,function(x) {
#     #   min(((x - ww[,1]) ^ 2) / ww[,2])
#     # })
#     
#     if (i == 1) .sed <- wm2
#     else .sed <- .sed + wm2
#   }
#   
#   sqrt(.sed)
# }
#-----------------

#--------
# .ncR <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
#   xx <- list()
#   if (!is.null(p)) xx <- c(p,xx)
#   if (!is.null(tmin)) xx <- c(tmin,xx)
#   if (!is.null(tmax)) xx <- c(tmax,xx)
#   if (!is.null(tmean)) xx <- c(tmean,xx)
#   
#   for (i in 1:length(xx)) {
#     x <- xx[[i]]
#     x1 <- x[[t1]]
#     x2 <- x[[t2]]
#     mn1 <- app(rast(x1),mean,na.rm=TRUE)
#     mn2 <- app(rast(x2),mean,na.rm=TRUE)
#     v <- app(rast(x1),var,na.rm=TRUE)
#     
#     w <- which(!is.na(mn2[]))
#     ww <- data.frame(mn1[w],v[w])
#     
#     wm <- sapply(mn2[w][,1],function(x,...) {
#       min(((x - ww[,1]) ^ 2) / ww[,2])
#     },na.rm=TRUE)
#     wm2 <- rast(mn2)
#     wm2[w] <- wm
#     
#     # wm <- app(mn2,function(x) {
#     #   min(((x - ww[,1]) ^ 2) / ww[,2])
#     # })
#     # 
#     if (i == 1) .sed <- wm2
#     else .sed <- .sed + wm2
#   }
#   
#   sqrt(.sed)
# }
#-----------------