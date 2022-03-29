\name{novelClimate}
\alias{novelClimate}
\alias{novelClimate,SpatRasterTS-method}
\alias{novelClimate,SpatRaster-method}
\alias{novelClimate,RasterStackBrickTS-method}
\alias{novelClimate,RasterStackBrick-method}
\alias{novelClimate,missing-method}

\title{Novel Climates}

\usage{
  novelClimate(x,...,t1,t2)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of a climate variable}
  \item{...}{additional climate variables can be included as a Raster object or a Raster Time Series}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1; if \code{x} is missing, \code{t1} and \code{t2} can be list of climate variables in time 1 and time 2, respectively.}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2;  if \code{x} is missing, \code{t1} and \code{t2} can be list of climate variables in time 1 and time 2, respectively.}
}

\description{
Novel climate is one of the climate change metrics to quantify the dissimilarities between climate parameters between two time periods. This metric uses standard euclidean distance between each cell in Time 2 and all cells in Time 1 and  retains the minimum of those distances. The inter-annual standard deviation for each parameter is used for the standardization. The larger the score, the more dissimilar the climate in Time 2 is in relation to the global pool of potential climatic analogues.

}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
