resampleSpikes <- function(param, var.log)
# Removes the putative spike-in variability from the fitted values of the
# original data, and then randomly adds it back in (assuming variance refers to
# the log-volume with base '2'). Counts are rescaled such that quantiles are
# preserved between the original and new fitted values, given the fitted
# parameters to the spike-ins.
#
# written by Aaron Lun
# created 14 October 2015
# last modified 23 June 2017
{
	# Removing the effect of spike-in variability from the original library sizes.
	log.lib <- log2(param$totals)
	cur.var <- var(log.lib)
	if (cur.var < var.log) { warning("variance in original library sizes is lower than specified") }
    new.log.lib <- (log.lib - mean(log.lib)) * sqrt(max(0, 1 - var.log/cur.var)) + mean(log.lib)
	new.lib.size <- 2^new.log.lib

	# Adding spike-in variability back in, randomly. 
	fc <- 2^rnorm(length(new.lib.size), 0, sd=sqrt(var.log))
    new.lib.size <- new.lib.size * fc

    # Computing new fitted values by scaling. No need to change the dispersion, as we're not
    # systematically shifting the mean; there should be no need to consider the mean-variance relationship.
	new.fitted <- t(t(param$fitted) * new.lib.size/param$totals)
	new.counts <- edgeR::q2qnbinom(param$counts, input=param$fitted, output=new.fitted, dispersion=param$dispersion)

    new.counts <- round(pmax(new.counts, 0)) # Making it reasonably count-like.
    dim(new.counts) <- dim(param$counts)
    dimnames(new.counts) <- dimnames(param$counts)
	return(new.counts)
}

