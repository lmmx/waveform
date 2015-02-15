\name{transcriptsDO}
\alias{transcriptsDO}
\docType{data}
\title{
The residual dissolved oxygen data from Li & Klevecz(2005)[1]
}
\description{
The respiratory oscillation in yeast that involves the majority of the coding transcriptome, was perturbed with phenelzine (an antidepressant drug used in the treatment of affective disorders in humans) after 48 min[1]. This produced a lengthening in the period during the experiment. The experiment consists of online dissolved oxygen measurements taken during the time-series sampling
}
\usage{data(transcriptsDO)}
\format{
  A data frame with 1130 observations on the following 4 variables.
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
[1] A rapid genome-scale response of the transcriptional oscillator to perturbation reveals a period-doubling path to phenotypic change. Li CM, Klevecz RR. Proc Natl Acad Sci U S A. 2006 Oct 31;103(44):16254-9. Epub 2006 Oct 16.
}
\references{
}
\examples{
data(transcriptsDO)
## plot an overview
plot(transcriptsDO)
}
\keyword{datasets}

