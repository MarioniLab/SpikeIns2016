estimateVariance <- function(y, design, ..., ratios, getfit=FALSE) 
# Estimates the variance of the ratios, given the design matrix.
#
# written by Aaron Lun
# created 26 January 2016
# last modified 8 February 2017
{
    if (missing(ratios)) { ratios <- y$samples$ratio }
    if (missing(design)) { design <- .make_intercept(length(ratios)) }
    else { design <- .restore_rank(design) }
    fit <- lm.fit(y=ratios, x=design, ...)
    if (getfit) { return(fit) }

    var.est <- sum(fit$effects[-seq_len(fit$rank)]^2)/fit$df.residual
    df <- nrow(design) - ncol(design)
    attributes(var.est)$standard.error <- var.est/sqrt(df) * sqrt(2) # std.err for variance est.
    attributes(var.est)$df <- df # also returning d.f.
    return(var.est)
}

.make_intercept <- function(n) { as.matrix(rep(1, n)) }

.restore_rank <- function(X) {
    QR <- qr(X)
    X[,QR$pivot[seq_len(QR$rank)],drop=FALSE]
}

splitSpikes <- function(spike.counts) 
# Function to partition spike-in transcripts into two halves with roughly similar abundance distributions.
#
# written by Aaron Lun
# created 26 January 2016
{
    o <- order(rowSums(spike.counts))
    chosen <- rep(c(TRUE, FALSE), length.out=length(o))
    subsum1 <- colSums(spike.counts[o[chosen],,drop=FALSE])
    subsum2 <- colSums(spike.counts[o[!chosen],,drop=FALSE])  
    return(list(first=subsum1, second=subsum2))
}

diagnoseVariance <- function(y, groups, design) 
# This estimates the variance for each level of 'groups'.
# It then does pairwise comparisons between groups to check for significant differences.
#
# written by Aaron Lun
# created 26 January 2016
{
    if (missing(design)) { design <- .make_intercept(length(groups)) }
    ratios <- y$samples$ratio
    all.index <- split(seq_len(length(ratios)), groups)
    all.var <- list()
    all.df <- list()

    for (g in names(all.index)) {
        current <- all.index[[g]]
        cur.ratios <- ratios[current]
        cur.design <- .restore_rank(design[current,,drop=FALSE])
        all.var[[g]] <- estimateVariance(ratios=cur.ratios, design=cur.design)
        all.df[[g]] <- nrow(cur.design) - ncol(cur.design)
    }
    
    pairwise <- matrix(NA, length(all.var), length(all.var))
    for (g1 in seq_along(all.var)) { 
        for (g2 in seq_len(g1-1L)) {
            out <- testVariance(var1=all.var[[g1]], var2=all.var[[g2]],
                                df1=all.df[[g1]], df2=all.df[[g2]], type="two-sided")
            pairwise[g1,g2] <- out
            pairwise[g2,g1] <- out
        }
    }
    rownames(pairwise) <- colnames(pairwise) <- names(all.var)
    return(list(var=all.var, pval=pairwise))
}

decomposeVariance <- function(y, design) 
# This decomposes the variances to their relevant components, given the count matrix in 'counts'.
# We need an indication of which genes are in each spike-in set (spike1, spike2).
# We also need to know which samples correspond to premixed and separate additions.
#
# written by Aaron Lun
# created 26 January 2016
# last modified 23 November 2016
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
        volume.sig <- testVariance(var1=total.var, var2=premix.var, design1=separate.design, 
                                   design2=premixed.design, type="one-sided")
    }
    
    # Computing standard error of volume variance.
    attributes(volume.var)$standard.error <- sqrt(attributes(total.var)$standard.error^2 + attributes(premix.var)$standard.error^2)/2

    # Returning all results.
    return(list(total=total.var, premixed=premix.var, volume=volume.var, pval=volume.sig))
}

