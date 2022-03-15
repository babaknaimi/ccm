# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  Dec. 2021
# Version 2.1
# Licence GPL v3
#--------


.eeChange <- function(xx, t1, t2, extreme) {
  # local extreme event change
  # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
  x1.t1 <- xx[[1]][[t1]]
  x1.t2 <- xx[[1]][[t2]]
  x2.t1 <- xx[[2]][[t1]]
  x2.t2 <- xx[[2]][[t2]]
  #------
  .ext1 <- app(x1.t1@raster,function(x) quantile(x,extreme[1],na.rm=TRUE))
  
  .pt <- app(c(x1.t2@raster,.ext1),function(x,...) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      1 - (length(which(x >= .ex)) / length(x))
    } else NA
  })
  #--------
  .ext2 <- app(x2.t1@raster,function(x) quantile(x,extreme[2],na.rm=TRUE))
  
  .pp <- app(c(x2.t2@raster,.ext2),function(x,...) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      length(which(x <= .ex)) / length(x)
    } else NA
  })
  
  extreme[1] - ((.pt + .pp) - (.pt * .pp))
}
#-------
.eeChangeR <- function(xx, t1, t2, extreme) {
  # local extreme event change
  # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
  x1.t1 <- xx[[1]][[t1]]
  x1.t2 <- xx[[1]][[t2]]
  x2.t1 <- xx[[2]][[t1]]
  x2.t2 <- xx[[2]][[t2]]
  #------
  .ext1 <- app(x1.t1@raster,function(x) quantile(x,extreme[1],na.rm=TRUE))
  
  .pt <- app(c(x1.t2@raster,.ext1),function(x,...) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      1 - (length(which(x >= .ex)) / length(x))
    } else NA
  })
  #--------
  .ext2 <- app(x2.t1,function(x) quantile(x,extreme[2],na.rm=TRUE))
  
  .pp <- app(c(x2.t2,.ext2),function(x,...) {
    x <- x[!is.na(x)]
    l <- length(x)
    if (l > 5) {
      .ex <- x[l]
      x <- x[1:c(l-1)]
      length(which(x <= .ex)) / length(x)
    } else NA
  })
  
  extreme[1] - ((.pt + .pp) - (.pt * .pp))
}
# 
# .eeChange <- function(tmp, pre, t1, t2, extreme = 0.95) {
#   # local extreme event change
#   # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
#   tm1 <- tmp[[t1]]
#   tm2 <- tmp[[t2]]
#   p1 <- pre[[t1]]
#   p2 <- pre[[t2]]
#   #------
#   #.ext1 <- calc(tm1@raster,function(x) quantile(x,extreme,na.rm=TRUE))
#   .ext1 <- app(rast(tm1@raster),function(x) quantile(x,extreme,na.rm=TRUE))
#   
#   .pt <- app(rast(stack(tm2@raster,.ext1)),function(x,...) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       1 - (length(which(x >= .ex)) / length(x))
#     } else NA
#   })
#   #--------
#   .ext2 <- app(rast(p1@raster),function(x) quantile(x,1-extreme,na.rm=TRUE))
#   
#   .pp <- app(rast(stack(p2@raster,.ext2)),function(x,...) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       length(which(x <= .ex)) / length(x)
#     } else NA
#   })
#   
#   extreme - ((.pt + .pp) - (.pt * .pp)) 
# }
#------------------


#-----------------
# .eeChangeR <- function(tmp, pre, t1, t2, extreme = 0.95) {
#   # local extreme event change
#   # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
#   tm1 <- tmp[[t1]]
#   tm2 <- tmp[[t2]]
#   p1 <- pre[[t1]]
#   p2 <- pre[[t2]]
#   #------
#   .ext1 <- app(rast(tm1),function(x) quantile(x,extreme,na.rm=TRUE))
#   .pt <- app(rast(stack(tm2,.ext1)),function(x) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       1 - (length(which(x >= .ex)) / length(x))
#     } else NA
#   })
#   #--------
#   .ext2 <- app(rast(p1),function(x) quantile(x,1-extreme,na.rm=TRUE))
#   
#   .pp <- app(rast(stack(p2,.ext2)),function(x) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       length(which(x <= .ex)) / length(x)
#     } else NA
#   })
#   
#   extreme - ((.pt + .pp) - (.pt * .pp)) 
# }
#------------------