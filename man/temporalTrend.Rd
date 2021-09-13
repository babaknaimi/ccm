\name{temporalTrend}
\alias{temporalTrend}
\alias{temporalTrend,RasterStackBrickTS-method}
\alias{temporalTrend,RasterStackBrick-method}

\title{Temporal Trend (Slope) of a climate variable}

\usage{
  temporalTrend(x,...)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of climate variable}
  
  \item{...}{additional arguments; not implemented}
}

\description{
The method calculates the trend of a climate variale changes over time. The function use a Raster time series as input and returns a Raster object that represent the trend at the pixel level.
}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
