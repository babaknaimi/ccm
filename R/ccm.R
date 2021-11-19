# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Oct. 2021
# Version 2.0
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

.velocM <- function(p1,p2,f1,f2,log=TRUE,...) {
  # based on the script provided in Hamnan et al. (2015)
  w <- which(!is.na(p1[]))
  
  present1 <- as.data.frame(p1)[w,]
  present2 <- as.data.frame(p2)[w,]
  future1  <- as.data.frame(f1)[w,]
  future2  <- as.data.frame(f2)[w,]
  
  idxy <- data.frame(cbind(id=1:length(present1),xyFromCell(p1,w)) )  # data frame of IDs and XY coords
  b <- (max(present1)-min(present1))/120  # bin size for 120 PC1 bins
  
  pr1 <- round(present1/b)              # convert PC1 to 120 bins via rounding
  pr2 <- round(present2/b)              # convert PC2 to <120 bins via rounding
  fu1 <- round(future1/b)               # same for future PC1
  fu2 <- round(future2/b)               # same for future PC2
  p  <- paste(pr1,pr2)                         # PC1/PC2 combinations in present climate
  f  <- paste(fu1,fu2)                         # PC1/PC2 combinations in future climate
  u  <- unique(p)[order(unique(p))]          # list of unique PC1/PC2 combinations
  
  sid <- c()                                 # empty vector for source IDs
  tid <- c()                                 # empty vector for target IDs
  d   <- c()                                 # empty vector for distances
  
  for(i in u){                          # loop for each unique PC1/PC2 combination
    pxy <- idxy[which(p==i),]           # coordinates of i-th combination in present
    fxy <- idxy[which(f==i),]           # coordinates of i-th combination in future
    sid <- c(sid, pxy$id)               # append i-th PC1/PC2 combination to previous 
    
    if(nrow(fxy)>0){                    # kNN search unless no-analogue climate
      knn <- data.frame(yaImpute::ann(as.matrix(fxy[,-1]), as.matrix(pxy[,-1]), k=1,verbose=FALSE)$knnIndexDist)      
      tid <- c(tid, fxy[knn[,1],"id"]) # the IDs of the closest matches  
      d <- c(d, sqrt(knn[,2]))         # their corresponding geographic distances
    } else {                              # else statement for no-analogue climates
      tid <- c(tid, rep(NA,nrow(pxy))) # flag destinations as missing for no analogues
      d <- c(d, rep(Inf,nrow(pxy)))    # flag distances as infinity for no analogues
    }
  }
  
  sxy <- merge(sid, idxy, by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)[2:3]  # source coordinates
  txy <- merge(tid, idxy, by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)[2:3]  # target coordinates
  names(txy)=c("target_x","target_y")
  names(sxy)=c("x","y")
  outtab <- cbind(id=sid, sxy, txy, distance=d)   
  # writes out log10 velocities and distances multiplied by 100 in ESRI ASCII format
  # conversion: -200=0.01km, -100=0.1km, 0=1km, 100=10km, 200=100km etc.
  out=merge(idxy, outtab[,c(2,3,6)], by=c("y","x"), sort=F)
  out$distance[out$distance==Inf] <- 10000  # sets no analogue to 10,000km
  out$distance[out$distance==0] <- 0.5  # sets zero distance to 0.5km (1/2 cell size)
  #out$logDist=round(log10(out$distance)*100)
  #out$logSpeed=round(log10(out$distance/50)*100)
  r <- raster(p1)
  if (log) r[cellFromXY(r,out[,c(2,1)])] <- round(log10(out$distance)*100)
  else r[cellFromXY(r,out[,c(2,1)])] <- out$distance
  r
}

#---------
# .vel <- function(basr,futr,nyears,longlat=TRUE,...) {
#   # The function is based on the script that is kindly provided by Raquel Garcia following the paper Garcia et al. (2014)
#   nc = ncol(basr) 
#   nr = nrow(basr) 
#   resolution = xres(basr) 
#   basm = matrix(getValues(basr), nrow = nr, ncol = nc, byrow = TRUE)  
#   #-------
#   lats = yFromRow(basr,1:nrow(basr)) 
#   longDist = NULL 
#   for(i in 1:length(lats)) {
#     longDist[i] = (pointDistance(c(0,lats[i]),c(resolution,lats[i]),longlat=longlat)) / 1000
#   } 
#   #----
#   spatialg = matrix(NA, nrow=nr, ncol=nc) 
#   if (longlat) {
#     for(i in 2:(nr-1)) {
#       for(j in 2:(nc-1)) {
#         if(!is.na(basm[i,j])) {
#           xDist = longDist[i]
#           
#           NS = (basm[i-1,j] - basm[i+1,j]) / (2*111.3195*resolution)
#           EW = (basm[i,j-1] - basm[i,j+1]) / (2*xDist)
#           
#           spatialg[i,j] = sqrt(NS^2 + EW^2)
#         }
#       }
#     } 
#   } else {
#     for(i in 2:(nr-1)) {
#       for(j in 2:(nc-1)) {
#         if(!is.na(basm[i,j])) {
#           xDist = longDist[i]
#           
#           NS = (basm[i-1,j] - basm[i+1,j]) / (2*resolution)
#           EW = (basm[i,j-1] - basm[i,j+1]) / (2*xDist)
#           spatialg[i,j] = sqrt(NS^2 + EW^2)
#         }
#       }
#     } 
#   }
#   #-------
#   spatialgr <- raster(basr)
#   temporalgr <- futr - basr #overall temporal gradient between the two time periods
#   tyear <- temporalgr / nyears #temporal gradient per year; nyears is the number of years between the two time periods
#   
#   
#   spatialgr <- setValues(spatialgr, spatialg)
#   
#   
#   
#   ### computing climate change velocity
#   
#   ccvel <- abs(tyear) / spatialgr
#   
#   # for grid cells with spatial gradient of zero the velocity will be Inf; so here I truncate the spatial gradient so as not to have zero values. The point at which to truncate (0.00005 in the example below) will depend on your data as you're trying to get a good balance between not having many Inf but also not changing a lot of values.
#   
#   zero <- Which(spatialgr < 0.00005, cells=TRUE)
#   spatialgr[zero] <- 0.00005
#   
#   ccvelt <- abs(tyear) / spatialgr
#   ccvelt
# }
# 
# 





#----------
if (!isGeneric("ccm")) {
  setGeneric("ccm", function(p,tmin,tmax,tmean,stat,t1,t2,extreme,longlat,ny, ...)
    standardGeneric("ccm"))
}

setMethod('ccm', signature(p='RasterStackBrickTS'),
          function(p,tmin,tmax,tmean,stat,t1,t2,extreme=0.95,longlat,ny,...) {
            if (missing(p)) p <- NULL
            if (missing(tmin)) tmin <- NULL
            if (missing(tmax)) tmax <- NULL
            if (missing(tmean)) tmean <- NULL
            
            
            w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
            l1 <- c('p','tmin','tmax','tmean')[w[1]]
            
            if (missing(ny)) {
              ny <- nyears(eval(get(l1))@time)
            }
            
            if (missing(longlat)) {
              longlat <- is.projected(crs(eval(get(l1))@raster))
              if (is.na(longlat)) {
                longlat <- .is.projected(eval(get(l1))@raster)
              }
            }
            
            
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
              nl <- length(which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean))))
              if (nl > 2) {
                stop('This method of velocity is implemented based on the approach in Hamnan et al. (2014) that works based on two climate variable; For multiple variables, you may either use a PCA transformation and take the first two components, or use dVe (dVelocity) stat!')
              }
              #----
              
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              l1 <- c('p','tmin','tmax','tmean')[w[1]]
              l2 <- c('p','tmin','tmax','tmean')[w[2]]
              p1 <- calc(eval(get(l1))[[t1]]@raster, mean)
              f1 <- calc(eval(get(l1))[[t2]]@raster, mean)
              p2 <- calc(eval(get(l2))[[t1]]@raster, mean)
              f2 <- calc(eval(get(l2))[[t2]]@raster, mean)
              .velocM(p1,p2,f1,f2,...)
              
            } else if (stat %in% c('dve','dVe','dv','dVE','dVelocity')) {
              
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              
              if (length(w) == 1) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),t1=t1,t2=t2,ny=ny)
              else if (length(w) == 2) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),t1=t1,t2=t2,ny=ny)
              else if (length(w) == 3) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),get(c('p','tmin','tmax','tmean')[w[3]]),t1=t1,t2=t2,ny=ny)
              
            } else if (stat %in% c('gve','gVe','gv','gVE','gVelocity')) {
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              
              if (length(w) != 1) stop('The gVelocity implemented in this package works based on a single climate variable, so only provide one of the climate variables you wish to get gradiant-based velocity for!')
              
              gVelocity(get(c('p','tmin','tmax','tmean')[w[1]]))
              
              
            } else stop('stat is unknown...!')
          }
)



#----------------
setMethod('ccm', signature(p='RasterStackBrick'),
          function(p,tmin,tmax,tmean,stat,t1,t2,extreme,longlat,ny,...) {
            if (missing(p)) p <- NULL
            if (missing(tmin)) tmin <- NULL
            if (missing(tmax)) tmax <- NULL
            if (missing(tmean)) tmean <- NULL
            
            if (missing(longlat)) {
              longlat <- is.projected(crs(p))
              if (is.na(longlat)) {
                longlat <- .is.projected(p)
              }
            }
            
            
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
              nl <- length(which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean))))
              if (nl > 2) {
                stop('This method of velocity is implemented based on the approach in Hamnan et al. (2014) that works based on two climate variable; For multiple variables, you may either use a PCA transformation and take the first two components, or use dVe (dVelocity) stat!')
              }
              
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              l1 <- c('p','tmin','tmax','tmean')[w[1]]
              l2 <- c('p','tmin','tmax','tmean')[w[2]]
              p1 <- calc(get(l1)[[t1]], mean)
              f1 <- calc(get(l1)[[t2]], mean)
              p2 <- calc(get(l2)[[t1]], mean)
              f2 <- calc(get(l2)[[t2]], mean)
              .velocM(p1,p2,f1,f2,...)
              
            } else if (stat %in% c('dve','dVe','dv','dVE','dVelocity')) {
              
              if (missing(ny)) {
                if (!missing(dates)) ny <- nyears(dates)
                else stop('ny (number of years) should be specified...!')
              }
              
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              
              if (length(w) == 1) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),t1=t1,t2=t2,ny=ny)
              else if (length(w) == 2) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),t1=t1,t2=t2,ny=ny)
              else if (length(w) == 3) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),get(c('p','tmin','tmax','tmean')[w[3]]),t1=t1,t2=t2,ny=ny)
              
            } else if (stat %in% c('gve','gVe','gv','gVE','gVelocity')) {
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              
              if (length(w) != 1) stop('The gVelocity implemented in this package works based on a single climate variable, so only provide one of the climate variables you wish to get gradiant-based velocity for!')
              
              gVelocity(get(c('p','tmin','tmax','tmean')[w[1]]))
              
              
            } else stop('stat is unknown...!')
            
          }
)
