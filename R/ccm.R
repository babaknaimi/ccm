# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Dec. 2021
# Version 2.1
# Licence GPL v3
#--------



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
  setGeneric("ccm", function(x,...,stat,t1,t2,extreme,longlat,ny,dates,names)
    standardGeneric("ccm"))
}





setMethod('ccm', signature(x='SpatRasterTS'),
          function(x,...,stat,t1,t2,extreme=0.95,longlat,ny,dates,names) {
            xx <- list(x,...)
            
            if (missing(ny)) {
              ny <- nyears(xx[[1]]@time)
            }
            
            if (missing(longlat)) {
              longlat <- .is.projected(xx[[1]]@raster)
            }
            
            if (stat == 'sed') {
              .sed(xx,t1=t1,t2=t2)
            } else if (stat %in% c('leech','eech','exch','localExtreme')) {
              
              if (missing(extreme)) stop('extreme is needed...!')
              else if (length(extreme) == 1) {
                if (extreme > 1 || extreme < 0) stop('extreme should be within the range of 0 and 1')
                else {
                  extreme <- c(extreme,1-extreme)
                  warning('Only one value is provided for extreme that is assumed for the first variable, and for the second variable, 1-extreme is considered...!')
                }
              } else if (length(extreme) > 2) {
                extreme <- extreme[1:2]
                warning('More than two values are provided for extreme; Only the first two values are considered!')
              }
              
              if (length(xx) < 2) stop('two variables should be provided to ')
              else if (length(xx) > 2) {
                xx <- xx[[1:2]]
                warning('For eeChange metric, two variables are required (e.g., temperature and precipitation); Only the first two variables are considered!')
              }
              
              .eeChange(xx,t1=t1,t2=t2,extreme=extreme)
              
            } else if (stat == 'nc') {
              .n <- .nc(xx,t1=t1,t2=t2)
              .q <- global(a,quantile,probs=0.995,na.rm=TRUE)[1,1]
              .n <- ifel(.n >= .q,.q,.n)
              .n
            } else if (stat %in% c('dac','aac')) {
              
              if (missing(names)) stop('The argument of "names" is missing! names of the input variables should be provided...!')
              
              if (length(names) != length(xx)) stop('The length of the provided "names" is not the same of the number of the input variables!')
              
              
              names <- lower(names)
              
              
              w <- which(names %in% c('precipitation','prec','p','pr','precip'))
              
              if (length(w) == 0) stop('Precipitation is not provided (or its name is not identified)!')
              else if (length(w) > 1) stop('It seems two variables are provided that are related to precipitation; one precipitation variable is needed!')
              
              names[w] <- 'prec'
              #---
              w <- which(names %in% c('temp.min','tmin','tmn','tempmin','mintemp','min.temp'))
              if (length(w) == 0) stop('Minimum Temperature (tmin) is not provided (or its name is not identified)!')
              else if (length(w) > 1) stop('It seems two variables are provided that are related to Minimum Temperature (tmin); one variable is needed!')
              
              names[w] <- 'tmin'
              
              w <- which(names %in% c('temp.max','tmax','tmx','tempmax','maxtemp','max.temp'))
              if (length(w) == 0) stop('Maximum Temperature (tmax) is not provided (or its name is not identified)!')
              else if (length(w) > 1) stop('It seems two variables are provided that are related to Maximum Temperature (tmax); one variable is needed!')
              
              names[w] <- 'tmax'
              #----------
              w <- which(names %in% c('temp.mean','tmean','taverage','tempmean','meantemp','mean.temp','temp.average'))
              
              if (length(w) > 1) stop('It seems two variables are provided that are related to Mean Temperature (tmean); one variable is needed!')
              else if (length(w) == 1) names[w] <- 'tmean'
              else if (length(w) == 0) {
                names(xx) <- names
                
                .tmean <- (xx$tmin@raster + xx$tmax@raster) / 2
                .tmean <- rts(.tmean,index(xx$tmin))
                xx$tmean <- .tmean
                names <- c(names,'tmean')
              }
              #----------
              
              names(xx) <- names
              
              tmin1 <- xx$tmin[[t1]]
              tmin2 <- xx$tmin[[t2]]
              tmax1 <- xx$tmax[[t1]]
              tmax2 <- xx$tmax[[t2]]
              tmean1 <- xx$tmean[[t1]]
              tmean2 <- xx$tmean[[t2]]
              prt1 <- xx$prec[[t1]]
              prt2 <- xx$prec[[t2]]
              
              k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
              k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
              
              if (stat == 'dac') .disAnalogus(k1,k2)
              else .analogusClimate(k1,k2)
              
            } else if (stat == 've') {
              
              if (length(xx) == 1) stop('This method of velocity is implemented based on the approach in Hamnan et al. (2014) that works based on two climate variable; For single variable, you may use dVe (dVelocity) stat!')
              
              if (length(xx) > 2) stop('This method of velocity is implemented based on the approach in Hamnan et al. (2014) that works based on two climate variable; For multiple variables, you may either use a PCA transformation and take the first two components, or use dVe (dVelocity) stat!')
              #----
              p1 <- app(xx[[1]][[t1]]@raster, mean,na.rm=TRUE)
              f1 <- app(xx[[1]][[t2]]@raster, mean,na.rm=TRUE)
              
              p2 <- app(xx[[2]][[t1]]@raster, mean,na.rm=TRUE)
              f2 <- app(xx[[2]][[t2]]@raster, mean,na.rm=TRUE)
              
              .velocMTerra(p1,p2,f1,f2,...)
              
            } else if (stat %in% c('dve','dVe','dv','dVE','dVelocity')) {
              
              if (length(xx) == 1) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),t1=t1,t2=t2,ny=ny)
              else if (length(w) == 2) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),t1=t1,t2=t2,ny=ny)
              else if (length(w) == 3) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),get(c('p','tmin','tmax','tmean')[w[3]]),t1=t1,t2=t2,ny=ny)
              
            } else if (stat %in% c('gve','gVe','gv','gVE','gVelocity')) {
              w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
              
              if (length(w) != 1) stop('The gVelocity implemented in this package works based on a single climate variable, so only provide one of the climate variables you wish to get gradiant-based velocity for!')
              
              gVelocity(get(c('p','tmin','tmax','tmean')[w[1]]))
              
              
            } else stop('stat is unknown...!')
          }
)



# 
# 
# 
# setMethod('ccm', signature(p='SpatRasterTS'),
#           function(p,tmin,tmax,tmean,stat,t1,t2,extreme=0.95,longlat,ny,...) {
#             if (missing(p)) p <- NULL
#             if (missing(tmin)) tmin <- NULL
#             if (missing(tmax)) tmax <- NULL
#             if (missing(tmean)) tmean <- NULL
#             
#             
#             w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
#             l1 <- c('p','tmin','tmax','tmean')[w[1]]
#             
#             if (missing(ny)) {
#               ny <- nyears(eval(get(l1))@time)
#             }
#             
#             if (missing(longlat)) {
#               longlat <- is.projected(crs(eval(get(l1))@raster))
#               if (is.na(longlat)) {
#                 longlat <- .is.projected(eval(get(l1))@raster)
#               }
#             }
#             
#             
#             if (stat == 'sed') {
#               .sed(p,tmin,tmax,tmean,t1=t1,t2=t2)
#             } else if (stat %in% c('leech','eech','exch','localExtreme')) {
#               
#               if (is.null(p)) stop('precipitation (p) is not provided...!')
#               
#               if (missing(extreme)) extreme <- 0.95
#               
#               if (is.null(tmean)) {
#                 if (!is.null(tmax) & !is.null(tmin)) {
#                   tmean <- (tmin@raster + tmax@raster) / 2
#                   tmean <- rts(tmean,index(tmin))
#                 } else if (!is.null(tmax)) {
#                   tmean <- tmax
#                   warning('tmean is not provided, tmax is used instead...!')
#                 }  else if (!is.null(tmin)) {
#                   tmean <- tmin
#                   warning('tmean is not provided, tmin is used instead...!')
#                 } else {
#                   stop('temperature is not provided...!')
#                 }
#               }
#               
#               .eeChange(tmean,p,t1=t1,t2=t2,extreme=extreme)
#             } else if (stat == 'nc') {
#               .nc(p,tmin,tmax,tmean,t1=t1,t2=t2)
#             } else if (stat == 'dac') {
#               if (is.null(tmean)) {
#                 if (!is.null(tmax) & !is.null(tmin)) {
#                   tmean <- (tmin@raster + tmax@raster) / 2
#                   tmean <- rts(tmean,index(tmin))
#                 } else stop('Both tmin and tmax are needed...!')
#               }
#               #-----
#               tmin1 <- tmin[[t1]]
#               tmin2 <- tmin[[t2]]
#               tmax1 <- tmax[[t1]]
#               tmax2 <- tmax[[t2]]
#               tmean1 <- tmean[[t1]]
#               tmean2 <- tmean[[t2]]
#               prt1 <- p[[t1]]
#               prt2 <- p[[t2]]
#               
#               k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
#               k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
#               .disAnalogus(k1,k2)
#             } else if (stat == 'aac') {
#               if (is.null(tmean)) {
#                 if (!is.null(tmax) & !is.null(tmin)) {
#                   tmean <- (tmin@raster + tmax@raster) / 2
#                   tmean <- rts(tmean,index(tmin))
#                 } else stop('Both tmin and tmax are needed...!')
#               }
#               #-----
#               tmin1 <- tmin[[t1]]
#               tmin2 <- tmin[[t2]]
#               tmax1 <- tmax[[t1]]
#               tmax2 <- tmax[[t2]]
#               tmean1 <- tmean[[t1]]
#               tmean2 <- tmean[[t2]]
#               prt1 <- p[[t1]]
#               prt2 <- p[[t2]]
#               
#               k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
#               k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
#               .analogusClimate(k1,k2)
#               
#             } else if (stat == 've') {
#               nl <- length(which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean))))
#               if (nl > 2) {
#                 stop('This method of velocity is implemented based on the approach in Hamnan et al. (2014) that works based on two climate variable; For multiple variables, you may either use a PCA transformation and take the first two components, or use dVe (dVelocity) stat!')
#               }
#               #----
#               
#               w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
#               l1 <- c('p','tmin','tmax','tmean')[w[1]]
#               l2 <- c('p','tmin','tmax','tmean')[w[2]]
#               p1 <- calc(eval(get(l1))[[t1]]@raster, mean)
#               f1 <- calc(eval(get(l1))[[t2]]@raster, mean)
#               p2 <- calc(eval(get(l2))[[t1]]@raster, mean)
#               f2 <- calc(eval(get(l2))[[t2]]@raster, mean)
#               .velocM(p1,p2,f1,f2,...)
#               
#             } else if (stat %in% c('dve','dVe','dv','dVE','dVelocity')) {
#               
#               w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
#               
#               if (length(w) == 1) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),t1=t1,t2=t2,ny=ny)
#               else if (length(w) == 2) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),t1=t1,t2=t2,ny=ny)
#               else if (length(w) == 3) dVelocity(get(c('p','tmin','tmax','tmean')[w[1]]),get(c('p','tmin','tmax','tmean')[w[2]]),get(c('p','tmin','tmax','tmean')[w[3]]),t1=t1,t2=t2,ny=ny)
#               
#             } else if (stat %in% c('gve','gVe','gv','gVE','gVelocity')) {
#               w <- which(c(!is.null(p),!is.null(tmin),!is.null(tmax),!is.null(tmean)))
#               
#               if (length(w) != 1) stop('The gVelocity implemented in this package works based on a single climate variable, so only provide one of the climate variables you wish to get gradiant-based velocity for!')
#               
#               gVelocity(get(c('p','tmin','tmax','tmean')[w[1]]))
#               
#               
#             } else stop('stat is unknown...!')
#           }
# )



#----------------


setMethod('ccm', signature(x='RasterStackBrickTS'),
          function(x,...,stat,t1,t2,extreme=0.95,longlat,ny,dates,names) {
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
setMethod('ccm', signature(x='RasterStackBrick'),
          function(x,...,stat,t1,t2,extreme,longlat,ny,dates,names) {
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
