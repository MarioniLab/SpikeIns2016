spikeParam <- function(spike.counts) 
# This computes the parameters of the spike-in counts, including the
# fitted values and the dispersions used for rescaling later.
#
# written by Aaron Lun
# created 14 October 2015
# last modified 29 January 2016    
{
	spike.counts <- as.matrix(spike.counts)
	totals <- colSums(spike.counts)
	dispersion <- estimateDisp(DGEList(spike.counts, lib.size=totals), 
		design=cbind(rep(1, length(totals))), tagwise=FALSE)$trended.dispersion # Fits are usually good.
	central.mu <- exp(mglmOneGroup(spike.counts, dispersion=dispersion, offset=log(totals)))
	cur.mus <- outer(central.mu, totals)
	return(list(counts=spike.counts, lib.size=totals, mean=central.mu, fitted=cur.mus, dispersion=dispersion))
}

resampleSpikes <- function(param, var.log)
# Removes the putative spike-in variability from the fitted values of the
# original data, and then randomly adds it back in (assuming variance refers to
# the log-volume with base '2'). Counts are rescaled such that quantiles are
# preserved between the original and new fitted values, given the fitted
# parameters to the spike-ins.
#
# written by Aaron Lun
# created 14 October 2015
# last modified 4 May 2016 
{
	# Removing the effect of spike-in variability from the original library sizes.
	log.lib <- log2(param$lib.size)
	cur.var <- var(log.lib)
	if (cur.var < var.log) { warning("variance in original library sizes is lower than specified") }
	new.log.lib <- log.lib * sqrt(max(0, 1 - var.log/cur.var))
	new.log.lib <- new.log.lib - mean(new.log.lib) + mean(log.lib)
	new.lib.size <- 2^new.log.lib

	# Adding spike-in variability back in, randomly.
	fc <- 2^rnorm(length(new.lib.size), 0, sd=sqrt(var.log))
	new.fitted <- outer(param$mean, fc*new.lib.size)
	new.counts <- edgeR::q2qnbinom(param$counts, input=param$fitted, output=new.fitted, dispersion=param$dispersion)

    new.counts <- round(pmax(new.counts, 0)) # Making it reasonably count-like.
    dim(new.counts) <- dim(param$counts)
    dimnames(new.counts) <- dimnames(param$counts)
	return(new.counts)
}

