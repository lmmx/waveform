\name{oscilGet}
\alias{oscilGet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
oscilGet is a function for the statistical analyses of oscillatory time-series
}
\description{
This function performs the main statistical calculations and for permutation based analyses can be run in parallel by setting the nslaves.
}
\usage{
oscilGet(dset, statistic = c("OS", "SN", "ACF", "Fisher", "Robust", "Box"), permutations = "model", sd = NULL, denoise.level = 0, detrend = F, verbose = F, nSlaves = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dset}{
\code{dset} a vector or matrix of test values
}
  \item{statistic}{
\code{statistic} The statistics to return in the oscilGet list 
\itemize{
\item \code{'SN'}  signal-noise-ratio
\item \code{'OS'}  oscillation strength
\item \code{'ACF'} Autocorrelation
\item \code{'Fisher'} Fisher-g-test
\item \code{'Robust'}  Robust Fisher-g-test (not implemented correctly)
\item \code{'Box'}  Ljung Box test
}
}
  \item{permutations}{
can be an integer value to specify the number of variables or 'model' for a log normal distribution (needs to be verified for each dataset). Linear regression is used to estimate lower significance values in non-model mode (integer or 'model')
}
  \item{sd}{
\code{sd} standard deviation for the model (above), if not specified the sd for the input data is calculated (double or NULL)
}
  \item{denoise.level}{
 set n higher frequency components of the fft, to 0 (numeric)
}
  \item{detrend}{
 set n higher frequency components of the fft, to 0 (numeric)
}
  \item{verbose}{
 speak to the user during the process
}
  \item{nSlaves}{
 number of slave processes to start for the permutation analysis
}
  \item{\dots}{
 parameters to pass to DFT
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  List of the following structure:
  \item{$raw}{raw data} 
  \item{$filtered.data}{data filtered according to denoising and detrending} 
  \item{$DC}{average or DC component(s)} 
  \item{$residual}{raw data minus the DC component}
  \item{$rows.sd}{standard deviation of each row}
  \item{$window}{Fourier window used}
  \item{$denoise.level}{denoise level}
  \item{$detrend}{was the data detrended?}
  \item{$permutations}{permutation number or model based}
  \item{$sd}{standard deviation of the whole dataset}
  \item{$amplitude}{peak-to-peak amplitudes of the fourier components}
  \item{$angle}{phase-angles of the Fourier components }
  \item{$complex}{the fft output}
  \item{$power.order}{amplitude rank}
  \item{$SN$statistic}{the signal-to-noise ratio} 
  \item{$SN$p.value}{the significance of SN}
  \item{$OS$statistic}{the oscillation strength}
  \item{$OS$p.value}{the significance of OS}
  \item{$ACF$statistic}{the autocorrelation statistic for different lag values (see ACF in the 'stats' package)}
  \item{$ACF$p.value}{the significance of the ACF for different lag values}
  \item{$Box$statistic}{the Ljung-Box statistic for different lag values (see Box.test in the 'stats' package)}
  \item{$Box$p.value}{the significance of the Ljung-Box for different lag values}
  \item{$Fisher$dominant.freqs}{the major power detected (see fisher.g.test in the 'GeneCycle' package)}
  \item{$Fisher$p.value}{The significance of the dominant power}
  }
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
Phase-angles are, by construct, in reference to the first timepoint (timepoint 1 = phase-angle 0).

Please note that in the Fourier-based analyses (SN, OS) the column numbers indicate the frequency (column 2 = 2 cycles per dataset), while in the Autocorrelation methods (ACF, Box) the column numbers indicate the cycle period ( column 2 = a cycle every 2 timepoints).

The P values for SN, OS and ACF are calculated by this package. For Box and Fisher, the p-values are calculated by the original functions and the "permutations" argument has no influence. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{waveform}},\code{\link{DFT}},\code{\link{acf}}, \code{\link{Box.test}}, \code{\link{fisher.g.test}}
}
\examples{



## Perform a full oscilGet on the transcript dataset[1] using a log normal distribution model
data(transcripts)
mset <- oscilGet(transcripts)

## Find the dominant frequencies in the dataset
signif.counts = colSums(mset$SN$p.value < 0.01) # counts the number of genes with a significant oscillation at each frequency
names(signif.counts)= 1:ncol(mset$SN$p.value)
signif.counts 

## Compare significant 4-cycle oscillators identified by different methods
table(SN = mset$SN$p.value[,4] <0.01, OS = mset$OS$p.value[,4] < 0.01)
table(SN = mset$SN$p.value[,4] <0.01, ACF = mset$ACF$p.value[,12] < 0.01) #4 cycles per dataset of 48 timepoints  = 12 timepoints per cycle (frequency 4 = lag 12)

## Phase-angle histogram for genes with significant 4-cycle component 
oscil4 = mset$SN$p.value[,4] < 0.01
angles4 = mset$angle[oscil4, 4]
hist(angles4, xlim=c(-pi,pi), breaks=24)

## Adjust phase-angles with respect to the DO angle reference

ref.angle = transcriptsDO$angle[transcriptsDO$sample == 1] #the phase-angle of the first timepoint with respect to the respiratory cycle

library(circular) #for minusPiPlusPi function which returns angles in [-pi,pi)
angles4 = as.numeric(minusPiPlusPi(angles4+ref.angle)) 

## Re-plot the phase-angle histogram with adjusted angles and DO trace
par(mar=c(5,4,4,4))
plot(transcriptsDO$angle, transcriptsDO$DO, col="grey50",pch=20, xlim=c(-pi,pi) ,axes=F,ylab="",xlab="")
axis(side=4)
mtext("DO", side=4, line=2)
par(new=T)
hist(angles4, xlim=c(-pi,pi), breaks=24, xlab = "Phase-angle", ylab="Gene count", main ="Phase-angles of 4-cycle oscillators")


#Calculate the false discovery rate on these data
library(fdrtool)
fdr <- fdrtool(x=mset$SN$p.value[,4],statistic ='pvalue', plot=T)
summary(fdr$lfdr)

# With Bonferroni correction for p-values
pvalue <- p.adjust(mset$SN$p.value[,4],method ='bonferroni')
fdr = fdrtool(x=pvalue,statistic ='pvalue', plot=T)
summary(fdr$lfdr)

## Perform a signal-to-noise ratio oscilGet on the transcript dataset[1] using a permutation based analysis. Adjust the nSlaves to the number of cores on your machine. Takes a wee while.
pset <- oscilGet(transcripts, statistic='SN', permutations=10000, nSlaves=2)

#Compare plots the result for 4 cycles
plot(1-pset$SN$p.value[,4], 1-mset$SN$p.value[,4])
plot(log(pset$SN$p.value[,4],10), log(mset$SN$p.value[,4],10), ylim=c(-10,0), xlim=c(-10,0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ oscillation }
\keyword{ statistic }
\keyword{ DFT }
\keyword{ acf }
\keyword{ autocorrelation }% __ONLY ONE__ keyword per line
