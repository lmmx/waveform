\name{movingAverage}
\alias{movingAverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A simple moving average function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
movingAverage(x, n = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{x} A vector of values
}
  \item{n}{
\code{n} number of points to smooth
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector of smoothed values
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, n = 5) 
{
    filter(x, rep(1/n, n), sides = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
