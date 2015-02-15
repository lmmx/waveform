\name{transcripts}
\alias{transcripts}
\docType{data}
\title{
The raw transcript abundance data from Li & Klevecz(2005)[1]
}
\description{The respiratory oscillation in yeast[2] that involves the majority of the coding transcriptome, was perturbed with phenelzine (an antidepressant drug used in the treatment of affective disorders in humans) after 48 min[3]. This produced a lengthening in the period during the experiment. The experiment consists of 48 samples spaced at 4 min intervals
}
\usage{data(transcripts)}
\format{
  A datatable frame with 5570 Genes on the following 48 timepoints.
  \describe{
    \item{\code{Gene ID}}{row names}
    \item{\code{Sample ID}}{column names}
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
data(transcripts)

## plot transcript used in the manuscript
plot(unlist(transcripts["YAL067C",]))

## plot the distribution
hist(transcripts, breaks=256)

## apply a filter
dset <- transcripts #place any data filter here
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
