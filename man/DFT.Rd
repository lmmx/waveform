\name{DFT}
\alias{DFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DFT is a wrapper around R's default fft command
}
\description{
\code{DFT} provides a wrapper to combine R's fft and mvfft functions and to output complex(default) peak-to-peak amplitude, phase angle, mean (DC component) or ordered powers of a vector or matrix. Edge effects can be minimised by adding data windows. Noise can be filtered and data can be detrended.
}
\usage{
DFT(dset, inverse = F, output = "complex", window = "none", detrend = FALSE, denoise.level = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dset}{A numerical vector or matrix with time-series datasets}
  \item{inverse}{Perform the inverse Fourier transform (complex)}
  \item{output}{Specifies the output type of the \code{DFT}:
    \itemize{
    \item{'complex' is the default fft (complex)}
    \item{'angle' is the phase angle (atan2(-Im(x),Re(x))),  ( values between -pi and pi)}
    \item{'amplitude' the peak-to-peak amplitude (abs(x)*4/ncol(dset))}
    \item{'DC' the mean or DC component}
    \item{'power.order' ranked order of the amplitudes}
    }
  }
  \item{window}{
    A data window to minimise spectral leaking:
      \itemize{
        \item{'hanning'}
        \item{'hamming'} 
        \item{'rectangle'}
        }
  }
  \item{detrend}{set the first component, non-DC components of the fft to 0(boolean)}
  \item{denoise.level}{set n higher frequency components of the fft, to 0 (numeric)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric vector or matrix depending on input and output selection}
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
\code{\link{fft}}, \code{\link{oscilGet}}, \code{\link{waveform}}, \code{\link{hamming.window}}, \code{\link{hanning.window}}, \code{\link{rectangle.window}}
}
\examples{
## perform classical mvfft on a transcriptome experiment [2]
data(transcripts)
complex <- DFT(transcripts)
## get angles, amplitudes or mean
angles <- DFT(transcripts, output = 'angle')
amplitude <- DFT(transcripts, output = 'amplitude')
mean <- Re(DFT(transcripts, output = 'DC'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Fast}
\keyword{Fourier}
\keyword{Transform}
\keyword{Wrapper}
\keyword{Discrete}
\keyword{Angle}
\keyword{Amplitude}
\keyword{DC}% __ONLY ONE__ keyword per line
