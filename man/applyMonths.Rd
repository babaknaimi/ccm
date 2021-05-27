\name{apply.months}
\alias{apply.months}
\alias{apply.months,RasterStackBrickTS-method}
\alias{apply.months,RasterStackBrick-method}
\alias{apply.months,SpatRaster-method}

\title{Apply a function on rasters of each month over years}

\usage{
  apply.months(x,FUN,...)
}

\arguments{
  \item{x}{A Raster or Raster Time Series object}
  \item{FUN}{A function to apply on rasters of each month over different years; The defauls if \code{mean}}
  \item{...}{additional argument: if x is a Raster or SpatRaster object, corresponsing dates should be specified in the \code{dates} argument}
  
}

\description{
the output contains 12 raster layers corresponding to the 12 months (Jan-Dec.). For each of the available 12 months (from January to December), the corresponding rasters are selected and the \code{FUN} is applied. For example if the \code{FUN} is \code{mean}, the first layer in the output would be the mean of the rasters in January over different years, and the last later (12th layer) would be the mean of the rasters in December over different years.

}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
