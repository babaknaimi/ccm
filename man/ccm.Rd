\name{ccm}
\alias{ccm}
\alias{ccm,RasterStackBrickTS-method}
\alias{ccm,RasterStackBrick-method}

\title{Climate Change Metrics}

\usage{
  ccm(p,tmin,tmax,tmean,stat,t1,t2,extreme,...)
}

\arguments{
  \item{p}{a Raster object or a Raster Time Series of precipitation}
  \item{tmin}{a Raster object or a Raster Time Series of minimum temperature }
  \item{tmax}{a Raster object or a Raster Time Series of maximum temperature }
  \item{tmean}{a Raster object or a Raster Time Series of mean temperature }
  \item{stat}{a character specifying the metric used to generate climate change metrics, can be either of 'sed', 'leech'}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{extreme}{default: 0.95; the percentile used for stat='leech', local climate extreme event change}
  \item{...}{additional arguments; if stat is 'dec' and the input climate data are not raster time series, the corresponding dates/times should be provided in the \code{dates} argument}
}

\description{
multiple climate change metrics are available to quantify multiple dimensions of climate change between two times.

The metrics include:

- sed: Standardized Local Anomalies

- leech: Changes in probability of local extremes

- aac : Changes in areas of analogous climates

- nc : Novel climates

- dac: Changes in distance to analogous climates

- ve: Climate change velocity


}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
