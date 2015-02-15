\name{metabolitesDO}
\alias{metabolitesDO}
\docType{data}
\title{
The residual dissolved oxygen data from  from Kalesh Sasidharan, Tomoyoshi Soga, Masaru Tomita, Douglas B. Murray(2012)[1]
}
\description{
The dataset consists of the residual dissolved oxygen data measured during time-series sampling during the respiratory oscillation. Samples were extracted[1] and measured using capillary electrophoresis time-of-flight mass spectrometry (CE-ToFMS).
}
\usage{data(metabolitesDO)}
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
[1] A yeast metabolite extraction protocol optimised for time-series analyses. Sasidharan K, Soga T, Tomita M, Murray DB. PLoS One. 2012;7(8):e44283. doi: 10.1371/journal.pone.0044283. Epub 2012 Aug 29.
}
\references{

}
\examples{
data(metabolitesDO)
## plot an overview
plot(metabolitesDO)
}
\keyword{datasets}
