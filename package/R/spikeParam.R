spikeParam <- function(spike.counts, design=NULL) 
# This computes the parameters of the spike-in counts, including the
# fitted values and the dispersions used for rescaling later.
#
# written by Aaron Lun
# created 14 October 2015
# last modified 23 June 2017
{
	spike.counts <- as.matrix(spike.counts)
    if (is.null(design)) { design <- .make_intercept(ncol(spike.counts)) }
	totals <- colSums(spike.counts)
    offsets <- log(totals)

	dispersion <- estimateDisp(spike.counts, offset=offsets, design=design, prior.df=0)$tagwise.dispersion
    fitted <- glmFit(spike.counts, design=design, offset=offsets, dispersion=dispersion, prior.count=0)$fitted.values
	return(list(counts=spike.counts, totals=totals, fitted=fitted, dispersion=dispersion))
}

