\name{waveform}
\alias{waveform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a spectrally filtered model and correlates this with the real data.
}
\description{
Waveform uses oscilGet in either a targetted or untargetted approach. In the untargetted approach the user supplies a numerical vector or matrix which then undergoes frequency decomposition. The significant powers (either specified by a signal-to-noise level or the significance for this) are then used to reconstruct a model based waveform. This model is then linearly correlated with the real data to produce a goodness of fit (r^2) value and a Pearson correlation. From these values cut-offs. In a targetted approach, the user can specify targets and omissions and do a harmonic search so that the goodness of fit for waveform combinations can be quantified in oscillatory datasets.
}
\usage{
waveform(dset, cutoff, cutoff.method = "p.val", omit = NULL, target = NULL, search.harmonics = T, statistic = "SN", verbose = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dset}{
A numerical vector or matrix with time-series datasets
}
  \item{cutoff}{
hp-value or SN cutoff to use to create the model (numeric)
}
  \item{cutoff.method}{
 whether p-value or SN cutoffs are to be used ('p.value' or 'SN')
}
  \item{omit}{
omit certain frequencies
}
  \item{target}{
target certain frequencies
}
  \item{search.harmonics}{
searches the higher and neighbouring frequencies for harmonics
}
  \item{statistic}{
 extra statistics to return (see oscilGet, always returns SN)
}
  \item{verbose}{
speak to the user during the process
}
  \item{\dots}{
arguments to pass onto oscilGet
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns all oscilGet output plus a list of the following structure:
\item{$waveform$filtered.data}{the model time-series}
\item{$waveform$spectra}{the model spectra}
\item{$waveform$component.used}{logical values indicating whether a component was preserved or not}
\item{$waveform$pearson$correlation}{the model pearson correlation}
\item{$waveform$pearson$p.value}{the model person correlation significance}
\item{$waveform$r.squared$fit$$r.squared}{coefficient of determination (R2) of the model}
\item{$waveform$r.squared$fit$p.value}{coefficient of determination significance}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{oscilGet}},\code{\link{DFT}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
## The function is currently defined as


##Untargeted waveform algorithm. 
data(transcripts)
dset = waveform(transcripts)

##check goodness of fit
hist(dset$waveform$fit$r.squared)
poor.fit = dset$waveform$fit$r.squared < 0.5


##Visualize model 
library(RColorBrewer) #for spectral colors
library(scales) #for rescale

angle.order = order(dset$angle[!poor.fit,4]) #get order of genes by phase-angle, helps with visualisation
raw.data = dset$raw[!poor.fit,]
model = dset$waveform$filtered.data[!poor.fit,]

x=1:ncol(raw.data)
y=1:nrow(raw.data)
#rescale each line for visualisation
z.raw.data=apply(raw.data, 1, function(x) rescale(x, to=c(0,1)))
z.model=apply(model, 1, function(x) rescale(x, to=c(0,1)))

image(x,y,z.raw.data[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Raw data") 

image(x,y,z.model[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Model")

##identify major powers kept and group genes by their major components
comp.used = dset$waveform$component.used[!poor.fit,]
gene.count = colSums(comp.used)
maj.comp = which(gene.count >  100) #use components that are found in at least 100 genes

groups = apply(comp.used[,maj.comp], 1, function(x) paste(maj.comp[x], collapse=".") ) 
sort(table(groups), decreasing=T) #number of genes containing each combination of major powers

target.group = "1.3.4"
image(x,1:sum(groups == target.group ),z.raw.data[,groups == target.group ], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Raw data")
image(x,1:sum(groups == target.group ),z.model[,groups == target.group ], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Model")


##cluster by major components
to.cluster = cbind(Re(dset$waveform$spectra[!poor.fit, maj.comp+1]), Im(dset$waveform$spectra[!poor.fit,maj.comp+1]))
hc = hclust(dist(to.cluster),  method="ward.D")

image(x,y,z.raw.data[,hc$order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Raw data")
image(x,y,z.model[,hc$order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Model")


## Targeted approach


## Search for peaks with a 2-cycle component
data(metabolites)
dset = waveform(metabolites, target=2)

##remove very poor fits
hist(dset$waveform$fit$r.squared)
poor.fit = dset$waveform$fit$p.value > 0.01

##Plot identified oscillatory peaks
angle.order = order(dset$angle[!poor.fit,2]) #get order of genes by phase-angle, helps with visualisation
raw.data = dset$raw[!poor.fit,]
model = dset$waveform$filtered.data[!poor.fit,]

x=1:ncol(raw.data)
y=1:nrow(raw.data)
#rescale each line for visualisation
z.raw.data=apply(raw.data, 1, function(x) rescale(x, to=c(0,1)))
z.model=apply(model, 1, function(x) rescale(x, to=c(0,1)))

image(x,y,z.raw.data[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Raw data") 

image(x,y,z.model[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Gene ID", main="Model")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Oscillation }
\keyword{ Spectra }
\keyword{ Regression }
\keyword{ Linear }
\keyword{ Model }
\keyword{ pvalue }% __ONLY ONE__ keyword per line
