decomposeVariance <- function(y, design) 
# This decomposes the variances to their relevant components, given the count matrix in 'counts'.
# We need an indication of which genes are in each spike-in set (spike1, spike2).
# We also need to know which samples correspond to premixed and separate additions.
#
# written by Aaron Lun
# created 26 January 2016
# last modified 23 June 2017
{
    ratios <- y$samples$ratio
    separate <- y$samples$separate
    premixed <- y$samples$premixed

    separate.design <- .restore_rank(design[separate,,drop=FALSE])
    premixed.design <- .restore_rank(design[premixed,,drop=FALSE])
    total.var <- estimateVariance(ratios=ratios[separate], design=separate.design)
    premix.var <- estimateVariance(ratios=ratios[premixed], design=premixed.design)
    volume.var <- 0.5*(total.var - premix.var)

#    combined.design <- rbind(do.call(cbind, c(list(separate.design), rep(list(0), ncol(premixed.design)))),
#                             do.call(cbind, c(rep(list(0), ncol(separate.design)), list(premixed.design))))
    if (volume.var < 0) {
        # Following the REML definition when the subtraction is negative (reconstructing the matrix, just in case it didn't have premixed in it).
        # But maybe it's a bit too confusing to have to redefine everything: I think we'll just put up with negative values.
#        premix.var <- total.var <- estimateVariance(ratios=c(ratios[separate], ratios[premixed]), design=combined.design) 
        volume.var <- 0
        volume.sig <- 1
    } else {
        volume.sig <- testVariance(var1=total.var, var2=premix.var, type="one-sided")
    }
    
    # Computing standard error of volume variance.
    attributes(volume.var)$standard.error <- sqrt(attributes(total.var)$standard.error^2 + attributes(premix.var)$standard.error^2)/2

    # Returning all results.
    return(list(total=total.var, premixed=premix.var, volume=volume.var, pval=volume.sig))
}

