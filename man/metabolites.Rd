\name{metabolites}
\alias{metabolites}
\docType{data}
\title{
The raw peak intensity data from Kalesh Sasidharan, Tomoyoshi Soga, Masaru Tomita, Douglas B. Murray(2012)[1]
}
\description{
The dataset consists of 20 time-series samples obtained during the respiratory oscillation, that were extracted[1] and measured using capillary electrophoresis time-of-flight mass spectrometry (CE-ToFMS).
}
\usage{data(metabolites)}
\format{
  The data is a matrix of identified and aligned peaks using the \link{xcms} package:
  Row names have a unique peak identifier "\code{pol_mz_rt}", where \code{pol} describes the polarity (n or p), \code{mz} describes the mass-charge ratio and \code{rt} describes the retention time.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
[1] A yeast metabolite extraction protocol optimised for time-series analyses. Sasidharan K, Soga T, Tomita M, Murray DB. PLoS One. 2012;7(8):e44283. doi: 10.1371/journal.pone.0044283. Epub 2012 Aug 29.
}
\references{
}
\examples{
data(metabolites)

## plot the distribution
hist(metabolites, breaks=256)

## apply a filter
dset <-metabolites[apply(metabolites, 1, min) >= 500, ] #remove rows with low (noisy) values, for CE-MS this was around 500
hist(log2(dset), breaks=64)

## scale and log the data
dset <-log(as.numeric(unlist(dset/rowMeans(dset))), 2)

## plot the distribution
hist(dset, breaks = 64)

## generate a normal distribution
rset <- rnorm(length(dset), sd=sd(dset)) 

## check the distribution
qqplot(dset, rset)
abline(0,1)
}
\keyword{datasets}
