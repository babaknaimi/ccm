# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Dec. 2021
# Version 2.1
# Licence GPL v3
#--------


# standardized local anomalies


.sed <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1@raster,mean,na.rm=TRUE)
    mn2 <- app(x2@raster,mean,na.rm=TRUE)
    v <- app(x1@raster,var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
# -----------
.sedR <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1,mean,na.rm=TRUE)
    mn2 <- app(x2,mean,na.rm=TRUE)
    v <- app(x1,var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
#------
# .sed <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
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
#     if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
#     else .s <- .s + ((mn1 - mn2) ^ 2 / v)
#   }
#   
#   .s <- sqrt(.s)
#   .s <- ifel(is.infinite(.s),0,.s)
#   .s
# }

#--------


# standardized local anomalies
.sedR <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
  xx <- list()
  if (!is.null(p)) xx <- c(p,xx)
  if (!is.null(tmin)) xx <- c(tmin,xx)
  if (!is.null(tmax)) xx <- c(tmax,xx)
  if (!is.null(tmean)) xx <- c(tmean,xx)
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1,mean,na.rm=TRUE)
    mn2 <- app(x2,mean,na.rm=TRUE)
    v <- app(x1,var,na.rm=TRUE)
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
