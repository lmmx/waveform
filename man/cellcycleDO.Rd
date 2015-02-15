\name{cellcycleDO}
\alias{cellcycleDO}
\docType{data}
\title{
The residual dissolved oxygen data from data from Kalesh Sasidharan, Cornelia Amariei, Masaru Tomita and Douglas B. Murray[1]
}
\description{
The dataset consists of the residual dissolved oxygen data measured during time-series sampling during the respiratory oscillation. Flow cytometric measurements of peak height from 64 samples taken during the yeast respiratory oscillation. Measurements were taken during this experiment.
}
\usage{data(cellcycleDO)}
\format{
  A data frame with 457 observations on the following 4 variables.
  \describe{
    \item{\code{time.sec}}{time in seconds}
    \item{\code{DO}}{residual dissolved oxygen concentration (\%)}
    \item{\code{angle}}{the phase angle of the cycle wrt the minimum firts derivative of \code{DO} (°)}
    \item{\code{sample}}{a numeric vector indicating the sample ID's}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
[1] Rapid DNA, RNA and protein extraction protocols optimized for slow continuously growing yeast cultures. Sasidharan K, Amariei C, Tomita M, Murray DB. Yeast. 2012 Aug;29(8):311-22. doi: 10.1002/yea.2911. Epub 2012 Jul 5.
}
\references{

}
\examples{
data(cellcycleDO)
## plot an overview
plot(cellcycleDO)
}
\keyword{datasets}
