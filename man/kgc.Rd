\name{kgc}
\alias{kgc}
\alias{kgc,RasterStackBrick,RasterStackBrick,RasterStackBrick-method}
\alias{kgc,SpatRaster,SpatRaster,SpatRaster-method}

\title{Koggen Geiger climate classification}

\usage{
  kgc(p,tmin, tmax, tmean)
}

\arguments{
  \item{p}{A Raster object contains 12 months precipitation}
  \item{tmin}{A Raster object contains 12 months minimum temperature}
  \item{tmax}{A Raster object contains 12 months maximum temperature}
  \item{tmean}{A Raster object contains 12 months mean temperature}
}

\description{
If the input object is not a 12 month Raster object but a time series, the function first calls the \code{apply.months} to get the 12 month raster object.


}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
