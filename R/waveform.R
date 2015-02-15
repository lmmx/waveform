###Functions to generate sampled waveforms###
## R-SCRIPT START ##
#Filter a spectra to a target S:N ratio. If target is NA then the largest power is used. If targetted then the target is always kept
#1 is the DV component as is ignored
#should return a flat signal if nothing is above the user defined noise or if targetted the signal comprising of just the target.
#target is the number of cycles expected in the signal.
#user can specify multiple targets
#omit powers less than a specific frequency
#MERGED with oscilGet. R 2014-06-04 #source('~/Documents/R-scripts/oscilGet.R')
library('GeneCycle')
library('matrixStats')
library('foreach')
harmonics <- function(dset, spectra, target=NULL, search.harmonics= T, cutoff=2)
{
  if(!is.null(target) ) 
  {
    harmonics <- NULL
    for(i in 1:length(target)) # add the target if it's > the cutoff
    {
      if(dset[target[i]] > cutoff)
      {
        
        harmonics <- c(harmonics, target[i])
      }
    }
    if(search.harmonics)
    {
      #if(length(target) == 1) #if only one target specified try to find the harmonics and drift
      #catch adjacent closest power above the snr
      for(i in 1:length(target))
        if(target[i]-1 > 0)
        {
          if(dset[target[i]-1] > cutoff & dset[target[i]] > dset[target[i]-1]) 
          {
            harmonics <- c(target[i]-1, harmonics)
          }
          #now add back in all the downstream powers above the snr
          filter <- vector(mode='logical', length(dset))
          # this only targets those harmonics > the specified signal to noise ratio or p.val 
          # and adds the condition that the harmonic must have less power than the target
          ranges <- (target[i]+1):length(dset)
          filter[ranges] <- dset[ranges] >= cutoff & dset[ranges] < dset[target[i]]
          harmonics <- c(harmonics,  which(filter))
        }
      harmonics <- unique(harmonics)
    }
  }
  else
  {
    if(search.harmonics)
    {
      harmonics <- which(dset>cutoff)
    }
    else
    {
      if(max(dset) > cutoff)
      {
        harmonics <- which.max(dset)
      }
    }
  }
  stash <- vector(mode='numeric', length=length(spectra)) #makes a storage vector
  stash[1 + harmonics] <- spectra[1 + harmonics]  
  stash[1 + length(spectra) - harmonics] <- spectra[1 + length(spectra) - harmonics]
  return(stash)
}
## A waveform model generator based on signal to noise or SN p-value cutoffs
## Inputs dset is the dataset
## Statistic are the other statistics to pass onto the (see oscilGet). SN is required here
## method whether it is a signal-to-noise ratio (SN) or p.value cutoff (SN is the default)
## the required cutoff for each Fourier power. 0.05 is the p.val default. 2 is the p.val default
## omit are the powers to leave out of the calulations
## target are the powers to target. NULL uses an untargetted approach
## ... arguments passed onto the DFT (see DFT and oscilGet)
waveform <- function(dset, cutoff, cutoff.method='p.val', omit=NULL, target=NULL, search.harmonics=T, statistic='SN', verbose=F, ...)
{
  if(is.null(nrow(dset))) #transform a vector into a matrix
  {
    dset <- t(as.matrix(dset))
  }
  #calculates the waveform statistics
  if(!'SN' %in% statistic)
  {
    statistic <- c('SN', statistic)
  }
  dset <- oscilGet(dset=dset, statistic = statistic,  ...) # gets the SN
  #checks the method that needs to be used
  if(cutoff.method=='SN')
  {
    sset <- dset$SN$statistic
    if(missing(cutoff))
    {
      cutoff <- 2
    }
  }
  else if(cutoff.method=='p.val')
  {
    sset <- 1/dset$SN$p.value
    if(missing(cutoff))
    {
      cutoff <- 0.05
    }
    cutoff <- 1/cutoff
  }
  else
  {
    stop("invalid cutoff method. 'SN' or 'p.val'")
  }
  if(length(omit) >= 1){sset[,omit] <- 0} #omits the user defined periods
  #creates a list of harmonics in the data
  stash <-  t(sapply(1:nrow(dset$raw), function(i) harmonics(sset[i,,drop=F], dset$complex[i,,drop=F], target=target, cutoff=cutoff, search.harmonics = search.harmonics)))
  stash[,1] <- dset$complex[,1,drop=F] #adds back the DC component
  dset$waveform$filtered.data <- Re(DFT(stash,inverse=T))/ncol(dset$raw)
  colnames(dset$waveform$filtered.data) = colnames(dset$raw)
  dset$waveform$spectra <- stash
  colnames(dset$waveform$spectra) = paste("Power", 0:(ncol(dset$waveform$spectra)-1), sep=".")
  dset$waveform$component.used <- abs(dset$waveform$spectra[,2:floor(ncol(dset$waveform$spectra-1)/2)])>0
  #now run correlation tests
  dset$waveform$pearson$correlation <- sapply(1:nrow(dset$raw), function(i) suppressWarnings(cor(as.numeric(dset$raw[i,]), as.numeric(dset$waveform$filtered.data[i,]))))
  dset$waveform$pearson$p.value<- sapply(1:nrow(dset$raw), function(i) suppressWarnings(cor.test(as.numeric(dset$raw[i,]), as.numeric(dset$waveform$filtered.data[i,]))$"p.value"))
  dset$waveform$fit$p.value <- sapply(1:nrow(dset$raw), function(i) (fit.model(as.numeric(dset$raw[i,]), as.numeric(dset$waveform$filtered.data[i,]))))
  dset$waveform$fit$r.squared <- sapply(1:nrow(dset$raw), function(i) (summary(lm(as.numeric(dset$raw[i,]) ~ as.numeric(dset$waveform$filtered.data[i,]))))$r.squared)
  names(dset$waveform$fit$r.squared) = names(dset$waveform$fit$p.value) = names(dset$waveform$pearson$p.value) = names(dset$waveform$pearson$correlation) =  rownames(dset$waveform$filtered.data) = rownames(dset$waveform$spectra) = rownames(dset$waveform$component.used) = rownames(dset$raw)
  return(dset)
}
#returns the p.value on an r^2 linear model
fit.model <- function(x,y)
{
  statistic <- lm(formula = x ~ y)
  if(nrow(summary.lm(statistic)$coefficients) == 2)
  {
    return(summary.lm(statistic)$coefficients[2,4])
  }
  else
  {
    return(1)
  }
}

stat.chooser <- function(dset, statistic) #chooses the statistic and outputs it given a dataset
{
  if(statistic=='OS')
  {
    Areas <- sapply(1:ncol(dset$amplitude),function(x) dset$amplitude[,x,drop=F]*x*4)
    rset <- matrix(Areas/IntAbsRes(dset$filtered.data), ncol=ncol(dset$amplitude))
  }
  else if(statistic=='SN')
  {
    rset <-dset$amplitude/noise(dset$amplitude) ### ueda formula
    #rset <-(dset$amplitude/noise(dset$amplitude))^2 ### wiki formula
  }
  else if(statistic=='ACF')
  {
    rset <- (t(foreach(i=1:nrow(dset$filtered.data), .combine=cbind) %do% acf(unlist(dset$filtered.data[i,]), lag.max=ncol(dset$filtered.data), plot=FALSE)$acf[2:(ncol(dset$filtered.data))]))
colnames(rset) = paste("lag",1:ncol(rset),sep=".")
  }
else
{stop}
rownames(rset) <- rownames(dset$raw)
return(rset)
}

stat.perm <- function(dset, ...) #uses a list based dset to calculate a permutation matrix ... are fourier parameters
{
  stat.rows <- c(which(names(dset) == 'OS'), which(names(dset) == 'SN'), which(names(dset) == 'ACF'))
  r.dset <- NULL
  r.dset$filtered.data <- t(sapply(1:nrow(dset$filtered.data), function(i) sample(dset$filtered.data[i,]))) #rowperm(dset$filtered.data)
  r.dset$amplitude <- DFT(r.dset$filtered.data, output='amplitude', ...)
  model <- array(0,c(nrow(dset$amplitude), ncol(dset$amplitude), length(stat.rows)))
  for(i in 1:length(stat.rows))
  {
    index <- stat.rows[i]
    if(names(dset)[index] == 'ACF') #there was a bug here in the acf now its fixed pls check
    {
      model[,,i] <- abs(dset[[index]]$statistic) <= abs(stat.chooser(r.dset, names(dset)[index]))
    }
    else
    {
      model[,,i] <- dset[[index]]$statistic <= stat.chooser(r.dset, names(dset)[index]) #stat.chooser(r.dset, names(dset)[index])
    }
  }
  return(model)
}

compute.pvals <- function(dset, nSlaves=0, verbose=FALSE, sd, ...) #chooses wheter to do a normal based stat or a permutation stat
{
  stat.rows <- c(which(names(dset) == 'OS'), which(names(dset) == 'SN'), which(names(dset) == 'ACF'))
  if (is.integer(dset$permutations) | is.double(dset$permutations))
  {
    p.val <- function(x)
    {			
      if(nSlaves < 1) {nSlaves <- 1}
      cat("Calculating p values for",  
          names(x)[stat.rows], 
          "by", 
          x$method, 
          "using", 
          nSlaves, 
          "cores and", 
          x$permutations,
          'permutations. \n')
      #start.time = proc.time()
      #library("doSNOW", quietly=TRUE)
      cl <- makeCluster(nSlaves, type = "SOCK")
      clusterExport(cl, c("stat.perm", "DFT", "clean.data", "fourier.window", "stat.chooser", "IntAbsRes", "noise"))
      registerDoSNOW(cl)			
      pb <- txtProgressBar(max=x$permutations, style=3)
      p.vals <- stat.perm(x)
      permIndex <- 2
      while(permIndex <= x$permutations)
      {
        p.vals <- p.vals + foreach(1:nSlaves, .combine="+") %dopar% stat.perm(x)
        permIndex <- permIndex + nSlaves
        setTxtProgressBar(pb, permIndex)
      }
      close(pb)
      stopCluster(cl)
      #detach("package:doSNOW")
      p.vals <- p.vals/permIndex
      cat(dim(p.vals), '\n')
      for(i in 1:length(stat.rows))
      {
        index <- stat.rows[i]
        x[[index]]$p.value <- p.vals[,,i]
        colnames(x[[index]]$p.value) = colnames(x[[index]]$statistic)
        rownames(x[[index]]$p.value) = rownames(x[[index]]$statistic)
      }
      return(x)
    }
  }
  else if (dset$permutations == "model")
  {
    nrow <- nrow(dset$filtered.data)
    nsamples<-10000
    rdset <- 2^rnorm(n=nsamples*ncol(dset$filtered.data), sd=sd)
    rdset <- matrix(rdset, nrow=nsamples)
    rdset <- list(filtered.data=rdset, 
                  DC=rowMeans(rdset),
                  residual= data.frame(rdset - rowMeans(rdset)), 
                  rows.sd=apply(rdset, 1, sd), 
                  sd=sd,
                  amplitude=DFT(rdset, output='amplitude', ...),
                  angle=DFT(rdset, output='angle', ...),
                  complex=DFT(rdset, output='complex', ...),
                  power.order=DFT(rdset, output='power.order', ...))
    p.val <- function(x)
    {
      for(i in 1:length(stat.rows))
      {
        index <- stat.rows[i]
        #Make a quick get a  model distribution
        cat("Calculating p values for",  
            names(dset)[index], '\n')
        model <- stat.chooser(rdset, names(dset)[index])
        data <- x[[index]]$statistic
        temp <- matrix(NA, ncol=ncol(data), nrow = nrow(data))
        foreach(j=1:ncol(data)) %do% {
          sd <- sd(as.numeric(model[,j]))
          mean <- mean(as.numeric(model[,j]))
          #x[[index]]$p.value[,j] <- pnorm(data[,j], sd=sd, mean=mean, lower.tail=F)
          temp[,j] <-  pnorm(data[,j], sd=sd, mean=mean, lower.tail=F)
          if(verbose){pm(data[,j], mean=mean, sd=sd)}
        }
        x[[index]]$p.value = temp
      }
      colnames(x[[index]]$p.value) = colnames(x[[index]]$statistic)
      rownames(x[[index]]$p.value) = rownames(x[[index]]$statistic)
      return(x)
    }
  }
  else
  {
    stop("invalid method")
  }
  return(p.val(dset))
}

pm <- function(dset, mean, sd) #plots the model data
{
  x=dset
  hist((x), freq=FALSE, breaks=32)
  curve(pnorm(x, sd=sd, mean=mean, lower.tail=FALSE), add=TRUE, col='red')
  curve(dnorm(x, sd=sd, mean=mean), add=TRUE, col='blue')
}
false.positive <- function(statistic, statistic.id, p.val, verbose=FALSE) #generates false positives for data
{
  false.p <- fdrtool(p.val, statistic='pvalue', plot=verbose, verbose=verbose)
  lfdr <- false.p$lfdr
  q.val <- false.p$qval
  results <- cbind(statistic, p.val, q.val, lfdr)
  colnames(results) <- paste(statistic.id, colnames(results), sep = ".")
  return(results)
}

noise <- function(dset) #gets the noise from a given spectra clean up in an apply
{
  #output <- NULL
  
  #for(target in 1:ncol(dset))
  #	{
  #	output <- cbind(output, rowMeans(dset[,-target,drop=F]))
  #	}
  return(sapply(1:ncol(dset), function(i) rowMeans(dset[,-i, drop=F])))
}

IntAbsRes <- function(dset) # Returns the integral of the absolute residual values
{
  #if(is.matrix(dset) | is.array(dset) | is.list(dset) & !is.null(dim(dset)))
  #{
    return(rowSums(abs(dset - rowMeans(dset))))
  #}
  #else
  #{
  #  return(sum(abs(dset-mean(dset))))
  #}
}

DFT <- function(dset, inverse=F, output = 'complex', window='none', detrend=FALSE, denoise.level=0) #chooses the correct fourier to use for the dataset outputs can me 'complex', 'angle' or 'amplitude' or 'DC'
{
  start = 1
  dset <- clean.data(dset * fourier.window(ncol(dset), window), detrend, denoise.level)
  if(is.matrix(dset) | is.array(dset) | is.list(dset) & !is.null(dim(dset)))
  {
    divider = ncol(dset)
  }
  else
  {
    divider = length(dset)
  }
  if(output == 'complex')
  {
    FUN1 <- function(x){x}
    ranges <- function(x){1:x}
    FUN2 <- function(x){x}
    start=0
  }
  else if(output == 'angle')
  {
    FUN1 <- function(x){(atan2(-Im(x),Re(x)))}
    ranges <- function(x){2:(floor(x/2))}
    inverse <- FALSE
    FUN2 <- function(x){x}
  }
  else if(output == 'amplitude')
  {
    FUN1 <- function(x){abs(x)*4/divider}
    ranges <- function(x){2:(floor(x/2))}
    inverse=FALSE
    FUN2 <- function(x){x}
  }
  else if(output == 'DC')
  {
    FUN1 <- function(x){x}
    inverse=FALSE
    ranges <- function(x){0*x +1}
    FUN2 <- function(x){x}
  }
  else if(output == 'power.order')
  {
    FUN1 <- function(x){abs(x)*4/divider}
    inverse=FALSE
    ranges <- function(x){2:(floor(x/2))}
    FUN2 <- function(x){t(apply(x, 1, order, decreasing=T))} 
  }
  else if(output == 'inverse')
  {
    FUN1 <- function(x){x}
    ranges <- function(x){1:x}
    FUN2 <- function(x){x}
  }
  else
  {
    stop("Invalid output selected. Try DC, angle, amplitude or complex")
  }
  tmp <- FUN2(t(FUN1(mvfft(t(dset), inverse)))[,ranges(ncol(dset)), drop=F])
  colnames(tmp) <- paste(output, start:(ncol(tmp)-(1-start)), sep=".")
  rownames(tmp) <- rownames(dset)
  return (tmp)
}

#cyclohedron.order <- function(dset, verbose = FALSE) #Calculates the cyclohedron order for a given array
#	{
#	topoGraph.dir <- '/home/dougie/R-scripts/functions/topoGraph.R'
#	if(file.exists(topoGraph.dir))
#		{
#		source(topoGraph.dir)
#		cat("Calculating Cyclohedron order.", '\n')
#		return(cycleCounts(dset/rowMeans(dset))) #remove amplitude biases
#		}
#	else
#		{
#		cat('Please download the cyclohedron source code and point to the right directory', '\n')
#		return(NULL)
#		}
#	}

fourier.window <- function(ncol, window = 'none') # generates windows for data to remove edge effects of fourier transforms
{
  if(window == 'hanning')
  {
    return(matrix(e1071::hanning.window(ncol), nrow=1, ncol = ncol))
  }
  else if(window == 'hamming')
  {
    return(matrix(e1071::hamming.window(ncol), nrow=1, ncol = ncol))
  }
  else if(window == 'rectangle')
  {
    return(matrix(e1071::rectangle.window(ncol), nrow=1, ncol = ncol))
  }
  else
  {
    return(1)
  }
}

clean.data <- function(dset, detrend=	TRUE, denoise.level=1, verbose=FALSE ) # denoise and/or detrend an array
{
  if(!detrend & denoise.level==0)
  {
  }
  else
  {
    fft <- mvfft(t(dset)) #perform fourier
    if(detrend) # detrend
    {
      fft[2,,drop=F]  <- 0
      fft[nrow(fft),,drop=F]  <- 0
    }
    if(denoise.level > 0)  # Remove requested levels
    {
      lower.filter <- (floor(nrow(fft)/2)-denoise.level+1)
      upper.filter <- (ceiling(nrow(fft)/2) + denoise.level-1)
      fft[lower.filter:upper.filter, ] <- 0
    }
    dset <- t(Re(mvfft(fft, inverse=TRUE)))/nrow(fft) # reconstruct denoised/detrended data
  }
  return(dset)
}

predict.pvals <- function(value, pvalue, coverage = 0.2, verbose = TRUE) #uses linear regression fit of the value vs its log statistic to estimate p-vals =0
{
  x <- value[pvalue !=0] #removes zero values
  y <- log10(pvalue[pvalue !=0]) #removes zero values
  n.y <- ceiling(length(y)*coverage) #get the number of points to use for the fit
  if(n.y == 1)
  {
    n.y = length(y) #refactor y if too many 0's to use the entire length of the data
  }
  target <- NULL
  for(i in 1:n.y) #find the targets that encompass the coverage
  {
    target <- c(target, which.min(y))
    y[which.min(y)] = max(y) #temp move min to max
  }
  y <- log10(pvalue[pvalue !=0]) #resets y
  y <- y[target]
  x <- x[target]
  new <- data.frame(x = value)
  try({
    pred.w.plim <- predict(lm(y ~ x), new, interval="prediction")
    pred.w.clim <- predict(lm(y ~ x), new, interval="confidence")
    pvalue <- ifelse(pvalue !=0 , pvalue, 10^pred.w.plim)
    if(verbose)
    {
      matplot(new$x ,cbind(pred.w.clim, pred.w.plim[,-1]),lty=c(1,2,2,3,3), type="l", ylab="predicted y")
      points(value, log10(pvalue), col="red")
      points(x , y, col="grey", pch = 20)
    }
  })
  return(pvalue)
}

# Takes an double matrix with the time-series as columns
# and experimental values as rows, as well as the number
# of cycles suspected (it finds the maximum number of cycles if 0 is specified.
# and calculates phase angle, amplitude, oscillation strength,
# signal-to-noise ratio	dset.names  <- rownames(dset)os, maximal power, DC component, and optionally
# diverse statistics. The diverse statistics are supplied in a vector of strings.
# pvals for the oscillation strength (OS) and signal-to-noise ratios (SN) can
# either be calculated/passed based on the random.data function according to a normal
# distribution or be calculated from a permuted matrix of dset (random = "PERM").
# The random.size option is either the size of the random data passed (nrow(random)),
# the size of the generated array or the number of iterations of the permutation loop.
# WARNING:: the permutation calculation is computationally heavy.

oscilGet <- function(dset, 
                     statistic = c('OS', 'SN', 'ACF', 'Fisher',  'Box'),  #'Robust' not implemented correctly
                     permutations = 'model',
                     sd = NULL,
                     denoise.level =0,
                     detrend = F,
                     verbose = F, 
                     nSlaves = 0, 
                     ...)
{
  if(is.null(nrow(dset))) #transform a vector into a matrix
  {
    dset <- t(as.matrix(dset))
  }
  else
  {
    dset <- as.matrix(dset)
  }
{# DATA-CHECKING handle NAs Removes rows that have NAs
  if( any( is.na(rowSums(dset))) | any(is.nan(rowSums(dset)))) 
  {
    if(verbose){cat("WARNING: Check data integrity, data contains", sum(is.na(rowSums(dset))), "rows of NAs!", " ")
                cat("Omitting the following probes from calculation (the data will be reintegrated at the end of the analyses):", '\n')
                cat( rownames(dset[is.na(rowSums(dset)), , drop=F]), '\n' )}
    dset.na <- dset[ is.na(rowSums(dset)), ,drop=F]
    dset.order <- rownames(dset)
    dset <- dset[!is.na(rowSums(dset)), , drop=F]
    reconstruct <- TRUE
  }
  else
  {
    reconstruct <- FALSE
  }
}
#get to the nitty gritty. In this section dset becomes a list of results
cat("Calculating Fourier transform, noise-cut:", ifelse(denoise.level==0, "NO", denoise.level),", detrend:", detrend, '\n')
if(length(sd) == 0)
{
  sd <- sd(log2(as.numeric(dset/rowMeans(dset))[as.numeric(dset/rowMeans(dset)) >0]),  na.rm=T)
}
dset <- list(raw=dset,
             filtered.data=abs(DFT(DFT(dset, output='complex', ...), inverse=T))/ncol(dset), 
             DC=rowMeans(dset),
             residual= data.frame(dset - rowMeans(dset)), 
             rows.sd=apply(dset, 1, sd), 
             window=window, 
             denoise.level=denoise.level, 
             detrend=detrend,  
             permutations=permutations, 
             sd=sd,
             amplitude=DFT(dset, output='amplitude', ...),
             angle=DFT(dset, output='angle', ...),
             complex=DFT(dset, output='complex', ...),
             power.order=DFT(dset, output='power.order', ...))
#Now do the statistics
if('OS' %in% statistic)
{
  dset$OS$statistic <- stat.chooser(dset, 'OS')
}
if('SN' %in% statistic)
{
  dset$SN$statistic <- stat.chooser(dset, 'SN')
}
if('ACF' %in% statistic)
{
  dset$ACF$statistic <- stat.chooser(dset, 'ACF')
}
dset <- compute.pvals(dset, nSlaves=nSlaves, sd=sd, verbose=verbose, ...)
#	if(sum(statistic == 'Cyclo')==1)  #Oscillation Cyclohedron
#		{
#		dset$cyclo <- cyclohedron.order(dset$filtered.data, verbose=FALSE)
#		}
if('Box' %in% statistic)
{
  cat("Calculating p values for Ljung-Box\n")
  box.res<- lapply(1:nrow(dset$filtered.data), function(i)
    lapply(1:(ncol(dset$filtered.data)-1), function(x) Box.test(as.numeric(dset$filtered.data[i,]), type='Ljung', lag=x)[c("statistic","p.value")]))
  dset$Box$statistic = t(sapply(box.res, function(x) sapply(x, function(y) y$statistic)))
  dset$Box$p.value = t(sapply(box.res, function(x) sapply(x, function(y) y$p.value)))
  colnames(dset$Box$statistic) = colnames(dset$Box$p.value)= paste("lag",1:(ncol(dset$filtered.data)-1), sep=".")
  rownames(dset$Box$statistic) = rownames(dset$Box$p.value) = rownames(dset$filtered.data)
}
if('Fisher' %in% statistic)
{
  cat("Calculating p values for Fisher-g-test\n")
  dset$Fisher$dominant.freqs <- as.numeric(dominant.freqs(t(dset$filtered.data)))*ncol(dset$raw)
  dset$Fisher$p.value <- fisher.g.test(t(dset$filtered.data))
  
}
if('Robust' %in% statistic)
{
  cat("Calculating p values for Robust-g-test\n")
  dset$robust$statistic <- robust.spectrum(t(dset$filtered.data))
  dset$robust$p.value <- foreach(i=1:nrow(dset$robust$statistic), .combine="cbind") %do% robust.g.test(dset$robust$statistic, index=i)
  dset$robust$statistic <- t(dset$robust$statistic)
  rownames(dset$robust$statistic) = rownames(dset$robust$p.value) = rownames(dset$filtered.data)
  colnames(dset$robust$p.value) = colnames(dset$robust$statistic)
  unlink(paste("g_pop_length_",ncol(dset$filtered.data),".txt",		sep=""))
  unlink(paste("g_pop_length_",ncol(dset$filtered.data),"indexed.txt", sep=""))
}
if(reconstruct) 
{
  if('OS' %in% statistic)
  {
    dset$OS$statistic <- reconstructor(dset$OS$statistic, dset.na, dset.order, T)
    dset$OS$p.value <- reconstructor(dset$OS$p.value, dset.na, dset.order, T) 
  }
  if('SN' %in% statistic)
  {
    dset$SN$statistic  <- reconstructor(dset$SN$statistic, dset.na, dset.order, T) 
    dset$SN$p.value  <- reconstructor(dset$SN$p.value, dset.na, dset.order, T)
  }
  if('ACF' %in% statistic)
  {
    dset$ACF$statistic  <- reconstructor(dset$ACF$statistic, dset.na, dset.order, T) 
    dset$ACF$p.value  <- reconstructor(dset$ACF$p.value, dset.na, dset.order, T)
  }
  #		if(sum(statistic == 'Cyclo')==1)  #Oscillation Cyclohedron
  #			{
  #			dset$cyclo  <- reconstructor(dset$cyclo , dset.na, dset.order, T)
  #			}
  if ('Box' %in% statistic)
  {
    dset$Box$statistic  <- reconstructor(dset$Box$statistic, dset.na, dset.order, T) 
    dset$Box$p.value  <- reconstructor(dset$Box$p.value, dset.na, dset.order, T)
  }
  if('Fisher' %in% statistic)
  {
    dset$Fisher$statistic  <- reconstructor(dset$Fisher$statistic, dset.na, dset.order, T) 
    dset$Fisher$p.value  <- reconstructor(dset$Fisher$p.value, dset.na, dset.order, T)    
  }
  if('Robust' %in% statistic)
  {
    dset$robust$statistic  <- reconstructor(dset$robust$statistic, dset.na, dset.order, T) 
    dset$robust$p.value  <- reconstructor(dset$robust$p.value, dset.na, dset.order, T)
  }
  dset$filtered.data <- reconstructor(dset$filtered.data, dset.na, dset.order, F)
  dset$angle <- reconstructor(dset$angle, dset.na, dset.order, T)
  dset$complex <- reconstructor(dset$complex, dset.na, dset.order, T)
  dset$rows.sd <- reconstructor(dset$rows.sd, dset.na, dset.order, T)
  dset$power.order <- reconstructor(dset$power.order, dset.na, dset.order, T)
  dset$amplitude <- reconstructor(dset$amplitude, dset.na, dset.order, T)
  dset$DC <- reconstructor(dset$DC, dset.na, dset.order, T)
  dset$residual <- reconstructor(dset$residual, dset.na, dset.order, T)
}
return(dset)
}
reconstructor <- function(dset, dset.na, dset.order, fill.na)
{
  if(fill.na)
  {
    naMatrix <- matrix(NA, ncol = ncol(as.matrix(dset)), nrow=nrow(as.matrix(dset.na)))
    rownames(naMatrix) <- rownames(dset.na)
    dset.na <- naMatrix
  }
  dset <- rbind(as.matrix(dset), dset.na)
  return (dset[match(dset.order, rownames(dset)),,drop=F]) #put back in original order
}
movingAverage <- function(x, n=5){
  filter(x, rep(1/n,n), sides=2)} #moving average function 
sampledSin <- function(n,frequency, period, DC=0, amplitude=1, stdev=0, plot=T) #n is the number of samples, frequency is tha sampling frequency in seconds, period of the waveform in seconds
{
  duration <- n*frequency
  cycles <- duration/period
  tseq <- seq(from=0, to=(duration-frequency/100), by=frequency/100)
  #generate a radian sequence
  radseq <- seq(from=0, to=(2*pi*cycles-2*pi*cycles/n/100), by=2*pi*cycles/n/100)
  #sample sequence
  samseq <- seq(from=1, to=length(tseq), by=100)
  sim <- sin(radseq)
  #cat(rnorm(sim, sd))
  sim <- (sim+rnorm(length(sim), sd=stdev))*amplitude +DC #adds amplitude, noise and DC
  #now inject noise
  if(plot==T) 
  {
    plot(tseq, sim, pch=".", col='red', xlab="Time (s)", ylab="Amplitude")
    points(tseq[samseq], sim[samseq], col='blue')
  }
  return(cbind(tseq[samseq], sim[samseq]))
}

sampledSq <- function(n,frequency, period, high=0.5, DC=0, amplitude=1, stdev=0, plot=T) # high is ratio in the on state
{
  duration <- n*frequency
  cycles <- duration/period
  tseq <- seq(from=0, to=duration, by=frequency/100)
  #sample sequence
  samseq <- seq(from=1, to=length(tseq), by=100)
  samseq <- samseq[-length(samseq)]
  sim  <- NULL
  #counter <- 0
  for(i in 1:length(tseq))
  {
    #check to see which cycle
    if(i==1) 
    {
      cycle <- 1
      sim <- c(sim,1)
    }
    else if(cycle != ceiling(tseq[i]/period))
    {
      cycle <- ceiling(tseq[i]/period)
      sim <- c(sim,-1)
    }
    else
    {
      if(tseq[i] < period*(1-high) + period * (cycle-1))
      {
        sim <- c(sim,-1)
      }
      else
      {
        sim <- c(sim,1)
      }
    }
  }
  #make it more realistic by running a moving average
  sim <- movingAverage(sim, 100)
  sim <- (sim+rnorm(length(sim), sd=stdev))*amplitude +DC #adds amplitude, noise and DC
  if(plot==T) 
  {
    plot(tseq, sim, pch=".", col='red', xlab="Time (s)", ylab="Amplitude")
    points(tseq[samseq], sim[samseq], col='blue')
  }
  return(cbind(tseq[samseq], sim[samseq]))
}


## NOT RUN
#				window='none' 
#				denoise.level=0 
#				detrend=F 
#				permutations=10000 

## R-SCRIPT-END ##
## NOT RUN ##
#library(mclust)
#  Li <- read.csv('/home/dougie/Documents/Data-for paper/Li2006/annotated_data_vals.txt', sep=' ')
#  Li.data <- as.matrix(Li[,3:ncol(Li)])
#  #Li.stats <- oscilGet(Li.data)
#  #Li.fft <- t(mvfft(t(Li.data)))
#  mydata <- waveform(scale(Li.data), target=4, sn=0.075)
#  myclust <- Mclust(mydata$waveform, G=3:25)
#  pdf('clusters.pdf')
#	for(i in 1:myclust$G)
#		{
#		matplot(t(mydata$dset[myclust$classification==i,]), type='l')
#		}
#  dev.off()
#  
#  
#	data3c= get(load("./SavedObjects/data3c-3-35.rda"))
#	ds= data3c$red/data3c$green
#	ds = ds/(rowMedians(data3c$mnase)/rowMedians(data3c$green))
#	ds=avgStrands(ds)
#	ds=ds[str_detect(rownames(ds), "00P"),]
#  #Li.stats <- oscilGet(Li.data)
#  #Li.fft <- t(mvfft(t(Li.data)))
#  ld = li.data[complete.cases(li.data),]
#  ld=sweep(ld, 1,apply(ld,1,mean),"/")
#  #dout = filter.spectrum(as.matrix(ld))
#  mydata <- waveform(as.matrix(ld))
#   #load( file = paste("./SavedObjects/mydata.raw.rda"))
#  #myclust <- Mclust(mydata$waveform, G=3:7)
#  rownames(mydata$waveform) = rownames(mydata$dset)
#  #test = apply(mydata$waveform.spectra, 1, function(y) any(y[3:(ncol(mydata$waveform.spectra)/2)] != 0))
#  #hc = hclust(dist(mydata$waveform), method="ward.D")
#for(k in 4:10){
#	cluster=cutree(hc,k)
#		for(i in 1:k)
#			{
#			matplot(t(mydata$dset[test,][cluster==i,]), type='l')
#			}
#		matplot(t(mydata$dset[!test,]), type='l')
#	  dev.off()
#	  
#	  
#  }



## END NOT RUN


