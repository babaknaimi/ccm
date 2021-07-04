# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  July 2021
# Version 1.4
# Licence GPL v3
#--------


# standardized local anomalies

.sed <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
  xx <- list()
  if (!is.null(p)) xx <- c(p,xx)
  if (!is.null(tmin)) xx <- c(tmin,xx)
  if (!is.null(tmax)) xx <- c(tmax,xx)
  if (!is.null(tmean)) xx <- c(tmean,xx)
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    # mn1 <- calc(x1@raster,mean,na.rm=TRUE)
    # mn2 <- calc(x2@raster,mean,na.rm=TRUE)
    # v <- calc(x1@raster,var,na.rm=TRUE)
    mn1 <- app(rast(x1@raster),mean,na.rm=TRUE)
    mn2 <- app(rast(x2@raster),mean,na.rm=TRUE)
    v <- app(rast(x1@raster),var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  sqrt(.s)
}
#----------------
.nc <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
  xx <- list()
  if (!is.null(p)) xx <- c(p,xx)
  if (!is.null(tmin)) xx <- c(tmin,xx)
  if (!is.null(tmax)) xx <- c(tmax,xx)
  if (!is.null(tmean)) xx <- c(tmean,xx)
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(rast(x1@raster),mean,na.rm=TRUE)
    mn2 <- app(rast(x2@raster),mean,na.rm=TRUE)
    v <- app(rast(x1@raster),var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    # wm <- app(mn2,function(x) {
    #   min(((x - ww[,1]) ^ 2) / ww[,2])
    # })
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  sqrt(.sed)
}
#-----------------
#-----------------
.eeChange <- function(tmp, pre, t1, t2, extreme = 0.95) {
  # local extreme event change
  # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
  tm1 <- tmp[[t1]]
  tm2 <- tmp[[t2]]
  p1 <- pre[[t1]]
  p2 <- pre[[t2]]
  #------
  #.ext1 <- calc(tm1@raster,function(x) quantile(x,extreme,na.rm=TRUE))
  .ext1 <- app(rast(tm1@raster),function(x) quantile(x,extreme,na.rm=TRUE))
  
  .pt <- app(rast(stack(tm2@raster,.ext1)),function(x,...) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      1 - (length(which(x >= .ex)) / length(x))
    } else NA
  })
  #--------
  .ext2 <- app(rast(p1@raster),function(x) quantile(x,1-extreme,na.rm=TRUE))
  
  .pp <- app(rast(stack(p2@raster,.ext2)),function(x,...) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      length(which(x <= .ex)) / length(x)
    } else NA
  })
  
  extreme - ((.pt + .pp) - (.pt * .pp)) 
}
#------------------


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
    # mn1 <- calc(x1,mean,na.rm=TRUE)
    # mn2 <- calc(x2,mean,na.rm=TRUE)
    # v <- calc(x1,var,na.rm=TRUE)
    mn1 <- app(rast(x1),mean,na.rm=TRUE)
    mn2 <- app(rast(x2),mean,na.rm=TRUE)
    v <- app(rast(x1),var,na.rm=TRUE)
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  sqrt(.s)
}
#----------------
#----------------
.ncR <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
  xx <- list()
  if (!is.null(p)) xx <- c(p,xx)
  if (!is.null(tmin)) xx <- c(tmin,xx)
  if (!is.null(tmax)) xx <- c(tmax,xx)
  if (!is.null(tmean)) xx <- c(tmean,xx)
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(rast(x1),mean,na.rm=TRUE)
    mn2 <- app(rast(x2),mean,na.rm=TRUE)
    v <- app(rast(x1),var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    # wm <- app(mn2,function(x) {
    #   min(((x - ww[,1]) ^ 2) / ww[,2])
    # })
    # 
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  sqrt(.sed)
}
#-----------------
#-----------------
.eeChangeR <- function(tmp, pre, t1, t2, extreme = 0.95) {
  # local extreme event change
  # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
  tm1 <- tmp[[t1]]
  tm2 <- tmp[[t2]]
  p1 <- pre[[t1]]
  p2 <- pre[[t2]]
  #------
  .ext1 <- app(rast(tm1),function(x) quantile(x,extreme,na.rm=TRUE))
  .pt <- app(rast(stack(tm2,.ext1)),function(x) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      1 - (length(which(x >= .ex)) / length(x))
    } else NA
  })
  #--------
  .ext2 <- app(rast(p1),function(x) quantile(x,1-extreme,na.rm=TRUE))
  
  .pp <- app(rast(stack(p2,.ext2)),function(x) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      length(which(x <= .ex)) / length(x)
    } else NA
  })
  
  extreme - ((.pt + .pp) - (.pt * .pp)) 
}
#------------------

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
#-------

if (!isGeneric("ccm")) {
  setGeneric("ccm", function(p,tmin,tmax,tmean,stat,t1,t2,extreme,...)
    standardGeneric("ccm"))
}


setMethod('ccm', signature(p='RasterStackBrickTS'),
          function(p,tmin,tmax,tmean,stat,t1,t2,extreme=0.95,...) {
            if (missing(p)) p <- NULL
            if (missing(tmin)) tmin <- NULL
            if (missing(tmax)) tmax <- NULL
            if (missing(tmean)) tmean <- NULL
            
            if (stat == 'sed') {
              .sed(p,tmin,tmax,tmean,t1=t1,t2=t2)
            } else if (stat %in% c('leech','eech','exch','localExtreme')) {
              
              if (is.null(p)) stop('precipitation (p) is not provided...!')
              
              if (missing(extreme)) extreme <- 0.95
              
              if (is.null(tmean)) {
                if (!is.null(tmax) & !is.null(tmin)) {
                  tmean <- (tmin@raster + tmax@raster) / 2
                  tmean <- rts(tmean,index(tmin))
                } else if (!is.null(tmax)) {
                  tmean <- tmax
                  warning('tmean is not provided, tmax is used instead...!')
                }  else if (!is.null(tmin)) {
                  tmean <- tmin
                  warning('tmean is not provided, tmin is used instead...!')
                } else {
                  stop('temperature is not provided...!')
                }
              }
              
              .eeChange(tmean,p,t1=t1,t2=t2,extreme=extreme)
            } else if (stat == 'nc') {
              .nc(p,tmin,tmax,tmean,t1=t1,t2=t2)
            } else if (stat == 'dac') {
              if (is.null(tmean)) {
                if (!is.null(tmax) & !is.null(tmin)) {
                  tmean <- (tmin@raster + tmax@raster) / 2
                  tmean <- rts(tmean,index(tmin))
                } else stop('Both tmin and tmax are needed...!')
              }
              #-----
              tmin1 <- tmin[[t1]]
              tmin2 <- tmin[[t2]]
              tmax1 <- tmax[[t1]]
              tmax2 <- tmax[[t2]]
              tmean1 <- tmean[[t1]]
              tmean2 <- tmean[[t2]]
              prt1 <- p[[t1]]
              prt2 <- p[[t2]]
              
              k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
              k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
              .disAnalogus(k1,k2)
            } else if (stat == 'aac') {
              if (is.null(tmean)) {
                if (!is.null(tmax) & !is.null(tmin)) {
                  tmean <- (tmin@raster + tmax@raster) / 2
                  tmean <- rts(tmean,index(tmin))
                } else stop('Both tmin and tmax are needed...!')
              }
              #-----
              tmin1 <- tmin[[t1]]
              tmin2 <- tmin[[t2]]
              tmax1 <- tmax[[t1]]
              tmax2 <- tmax[[t2]]
              tmean1 <- tmean[[t1]]
              tmean2 <- tmean[[t2]]
              prt1 <- p[[t1]]
              prt2 <- p[[t2]]
              
              k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
              k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
              .analogusClimate(k1,k2)
              
            } else if (stat == 've') {
              
            } else stop('stat is unknown...!')
          }
)



#----------------
setMethod('ccm', signature(p='RasterStackBrick'),
          function(p,tmin,tmax,tmean,stat,t1,t2,extreme=0.95,dates,...) {
            if (missing(p)) p <- NULL
            if (missing(tmin)) tmin <- NULL
            if (missing(tmax)) tmax <- NULL
            if (missing(tmean)) tmean <- NULL
            
            if (stat == 'sed') {
              .sedR(p,tmin,tmax,tmean,t1=t1,t2=t2)
            } else if (stat %in% c('leech','eech','exch','localExtreme')) {
              
              if (is.null(p)) stop('precipitation (p) is not provided...!')
              
              if (missing(extreme)) extreme <- 0.95
              
              if (is.null(tmean)) {
                if (!is.null(tmax) & !is.null(tmin)) {
                  tmean <- (tmin + tmax) / 2
                } else if (!is.null(tmax)) {
                  tmean <- tmax
                  warning('tmean is not provided, tmax is used instead...!')
                }  else if (!is.null(tmin)) {
                  tmean <- tmin
                  warning('tmean is not provided, tmin is used instead...!')
                } else {
                  stop('temperature is not provided...!')
                }
              }
              
              .eeChangeR(tmean,p,t1=t1,t2=t2,extreme=extreme)
            } else if (stat == 'nc') {
              .ncR(p,tmin,tmax,tmean,t1=t1,t2=t2)
            } else if (stat == 'dac') {
              if (missing(dates)) stop('corresponding times/dates for raster layers are needed, should be provided in the "dates" argument...!')
              
              if (length(dates) != nlayers(p)) stop('The number of items in dates is not equal to the number of raster layers...!')
              
              if (is.null(tmean)) {
                if (!is.null(tmax) & !is.null(tmin)) {
                  tmean <- (tmin + tmax) / 2
                } else stop('Both tmin and tmax are needed...!')
              }
              #-----
              tmin1 <- tmin[[t1]]
              tmin2 <- tmin[[t2]]
              tmax1 <- tmax[[t1]]
              tmax2 <- tmax[[t2]]
              tmean1 <- tmean[[t1]]
              tmean2 <- tmean[[t2]]
              prt1 <- p[[t1]]
              prt2 <- p[[t2]]
              
              k1 <- kgc(apply.months(prt1,dates=dates[t1]),apply.months(tmin1,dates=dates[t1]),apply.months(tmax1,dates=dates[t1]),apply.months(tmean1,dates=dates[t1]))
              k2 <- kgc(apply.months(prt2,dates=dates[t2]),apply.months(tmin2,dates=dates[t2]),apply.months(tmax2,dates=dates[t2]),apply.months(tmean2,dates=dates[t2]))
              .disAnalogus(k1,k2)
            } else if (stat == 'aac') {
              if (missing(dates)) stop('corresponding times/dates for raster layers are needed, should be provided in the "dates" argument...!')
              
              if (length(dates) != nlayers(p)) stop('The number of items in dates is not equal to the number of raster layers...!')
              
              if (is.null(tmean)) {
                if (!is.null(tmax) & !is.null(tmin)) {
                  tmean <- (tmin + tmax) / 2
                } else stop('Both tmin and tmax are needed...!')
              }
              #-----
              tmin1 <- tmin[[t1]]
              tmin2 <- tmin[[t2]]
              tmax1 <- tmax[[t1]]
              tmax2 <- tmax[[t2]]
              tmean1 <- tmean[[t1]]
              tmean2 <- tmean[[t2]]
              prt1 <- p[[t1]]
              prt2 <- p[[t2]]
              
              k1 <- kgc(apply.months(prt1,dates=dates[t1]),apply.months(tmin1,dates=dates[t1]),apply.months(tmax1,dates=dates[t1]),apply.months(tmean1,dates=dates[t1]))
              k2 <- kgc(apply.months(prt2,dates=dates[t2]),apply.months(tmin2,dates=dates[t2]),apply.months(tmax2,dates=dates[t2]),apply.months(tmean2,dates=dates[t2]))
              .analogusClimate(k1,k2)
            } else if (stat == 've') {
              
            } else stop('stat is unknown...!')
            
          }
)
