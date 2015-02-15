\name{sampledSq}
\alias{sampledSq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A sampled square wave generator
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampledSq(n, frequency, period, high = 0.5, DC = 0, amplitude = 1, stdev = 0, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
 number of samples
}
  \item{frequency}{
 sampling frequency in seconds
}
  \item{period}{
  period of the waveform in seconds
}
  \item{high}{
 the ratio in the on state (0-1)
}
  \item{DC}{
 amplitude offset
}
  \item{amplitude}{
 peak-to-peak amplitude of the wave
}
  \item{stdev}{
noise injection (0 no noise 1, 1 standard deviation)
}
  \item{plot}{
plot the output
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix of sample time (seconds) and generated values.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sampledSin}}
}
\examples{
## Generate a square waveform 40 min sampled every 4 min
sampledSq(32,4*60,40*60)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{waveform}
\keyword{generator}% __ONLY ONE__ keyword per line

