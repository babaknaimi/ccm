\name{localExtreme}
\alias{localExtreme}
\alias{localExtreme,SpatRasterTS-method}
\alias{localExtreme,SpatRaster-method}
\alias{localExtreme,RasterStackBrickTS-method}
\alias{localExtreme,RasterStackBrick-method}
\alias{localExtreme,list,list-method}
\alias{localExtreme,missing,missing,list,list-method}

\title{Changes in probability of local extremes}

\usage{
  localExtreme(x1,x2,t1,t2,extreme)
}

\arguments{
  \item{x1}{Time Series of the first climate variable as a Raster or Raster Time Series Object}
  \item{x2}{Time Series of the second climate variable as a Raster or Raster Time Series Object; \code{x1} and \code{x2} can be lists of time series of climate variables in time 1 and time 2, respectively}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
}

\description{
Probability of changes in local extremes is calculated at each cell based on either one or two climate variables (e.g., precipitation and temperature). For each climate variable, the extreme value (percentile of the extreme in the variable distribution) should be specified by user in the extreme parameter. For instance, if temperature and precipitation are specified as the first and second inputs (x1 and x2), extreme=c(0.95,0.05) may be used that consideres 95th and 5th percentiles of the distribution of the temperature and precipitation time series, respectively over the first time period, and then it calculates the probability that the time1 percentile will be exceeded in time2. 
If two climate variables are used (e.g., temperature and precipitation), for each cell, the two probabilities are summed from which the product of the two probabilities are subtracted to obtain a measure of the probability of occurrence of either of the two extreme events (to avoid counting probabilities twice). The  probability of time_1 extreme climates in time_2, is then subtracted from the probability of time_1 extreme climates to obtain the changes in the probability of extreme climate variables. Positive values indicated increased probability in time 2, whereas negative values indicated a decrease. 

}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
