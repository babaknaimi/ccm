\name{ccm}
\alias{ccm}
\alias{ccm,RasterStackBrickTS-method}
\alias{ccm,RasterStackBrick-method}

\title{Climate Change Metrics}

\usage{
  ccm(p,tmin,tmax,tmean,stat,t1,t2,extreme,longlat,ny,...)
}

\arguments{
  \item{p}{a Raster object or a Raster Time Series of precipitation}
  \item{tmin}{a Raster object or a Raster Time Series of minimum temperature }
  \item{tmax}{a Raster object or a Raster Time Series of maximum temperature }
  \item{tmean}{a Raster object or a Raster Time Series of mean temperature }
  \item{stat}{a character specifying the metric used to generate climate change metrics, can be either of 'sed', 'leech', 'aac', 'dac', 've', 'gVe', 'dVe'}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{extreme}{default: 0.95; the percentile used for stat='leech', local climate extreme event change}
  \item{longlat}{TRUE or FALSE, specifies whether the input data is geographic or projected}
  \item{ny}{integer; specifies the number of years between time 1 and time2, if the input is Raster time series, it will be identified by the function}
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

- ve: Climate change velocity (based on the method represented in the paper Hamnan et al. (2015))

- dVe: Distance-based climate change velocity (based on the method used in Garcia et al. (2014))

- gVe: Gradiant-based climate change velocity (based on the method represented in Burrow et al. (2011))


References:

- Burrows MT, Schoeman DS, Buckley LB, Moore P, Poloczanska ES, Brander KM, Brown C, Bruno JF, Duarte CM, Halpern BS, Holding J. (2011) The pace of shifting climate in marine and terrestrial ecosystems. Science. Nov 4;334(6056):652-5.

- Garcia, R.A., Cabeza, M., Rahbek, C. and Araújo, M.B. (2014). Multiple dimensions of climate change and their implications for biodiversity. Science, 344(6183).

- Hamann, A., Roberts, D.R., Barber, Q.E., Carroll, C. and Nielsen, S.E. (2015). Velocity of climate change algorithms for guiding conservation and management. Global Change Biology, 21(2), pp.997-1004.



}
\examples{
#-------
path <- system.file("external/", package="ccm") # path to the dataset folder



}
