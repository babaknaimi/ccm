\name{sed}
\alias{sed}
\alias{sed,SpatRasterTS-method}
\alias{sed,SpatRaster-method}
\alias{sed,RasterStackBrickTS-method}
\alias{sed,RasterStackBrick-method}

\title{Standardized Local Anomalies}

\usage{
  sed(x,...,t1,t2)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of a climate variable}
  \item{...}{additional climate variables can be included as a Raster object or a Raster Time Series}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
}

\description{
Standardized Local Anomalies is one of the climate change metrics to quantify the change between two time periods.

}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
