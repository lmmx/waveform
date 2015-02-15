\name{cellcycle}
\alias{cellcycle}
\docType{data}
\title{
The G1 aligned DNA concentration data from Kalesh Sasidharan, Cornelia Amariei, Masaru Tomita and Douglas B. Murray[1].
}
\description{
The dataset consists of flow cytometric measurements of peak height from 64 samples taken during the yeast respiratory oscillation. Samples were taken 2.5 min apart.
}
\usage{data(cellcycle)}
\format{
  A matrix of 257 DNA concentration bins by 64 time series samples. Values represent the number of samples in the bin.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
[1] Rapid DNA, RNA and protein extraction protocols optimized for slow continuously growing yeast cultures. Sasidharan K, Amariei C, Tomita M, Murray DB. Yeast. 2012 Aug;29(8):311-22. doi: 10.1002/yea.2911. Epub 2012 Jul 5.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(cellcycle)

## plot the distribution
hist(cellcycle, breaks=256)

## apply a filter
dset <-cellcycle[apply(cellcycle, 1, min) >=5, ]
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
