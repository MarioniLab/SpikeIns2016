setupSpikes <- function(counts, spike1, spike2, separate, premixed, ...) 
# Returns a DGEList with the information formatted in a reasonably pretty manner.
# Avoids having to recompute things many times.
#
# written by Aaron Lun
# created 26 January 2016
{
    y <- DGEList(counts)
    y$samples$sum1 <- colSums(counts[spike1,,drop=FALSE])
    y$samples$sum2 <- colSums(counts[spike2,,drop=FALSE])
    y$samples$ratio <- log2(y$samples$sum1/y$samples$sum2)

    # Remembering which samples/genes do what.
    y$samples$separate <- .indexToLogical(counts, separate)
    y$samples$premixed <- .indexToLogical(counts, premixed)
    y$samples <- cbind(y$samples, ...)

    spike.data <- data.frame(spike1=.indexToLogical(counts, spike1, byrow=TRUE), 
                             spike2=.indexToLogical(counts, spike2, byrow=TRUE))
    if (!is.null(y$genes)) { 
        y$genes <- cbind(y$genes, spike.data)
    } else {
        y$genes <- spike.data
    }
    return(y)
}

.indexToLogical <- function(counts, index, byrow=FALSE) {
    if (byrow) { 
        afun <- rownames
        nfun <- nrow
    } else {
        afun <- colnames 
        nfun <- ncol
    }
    if (is.character(index)) { 
        index <- match(index, afun(counts))
    } 
    if (is.numeric(index)) {
        out <- logical(nfun(counts))
        out[index] <- TRUE
        index <- out
    }
    return(index)
}

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

testVariance <- function(var1, var2, df1, df2, ratios1, ratios2, design1, design2, type=c("one-sided", "two-sided"), ...) 
# Tests if the first variance is significantly larger than the second.
# If two-sided, it tests whether there are any significance differences between the variances.
#
# written by Aaron Lun
# created 26 January 2016
{
    if (missing(var1)) { 
        if (missing(design1)) { design1 <- .make_intercept(length(ratios1)) }
        var1 <- estimateVariance(ratios=ratios1, design=design1, ...) 
    } 
    if (missing(var2)) { 
        if (missing(design2)) { design2 <- .make_intercept(length(ratios2)) }
        var2 <- estimateVariance(ratios=ratios2, design=design2, ...) 
    }

    if (missing(df1)) { 
        if (missing(design1)) { design1 <- .make_intercept(length(ratios1)) }
        df1 <- nrow(design1) - ncol(design1) 
    }
    if (missing(df2)) { 
        if (missing(design2)) { design2 <- .make_intercept(length(ratios2)) }
        df2 <- nrow(design2) - ncol(design2) 
    }
    type <- match.arg(type)
    var.ratio <- var1/var2
    attributes(var.ratio) <- NULL

    upper.tail <- pf(var.ratio, df1, df2, lower=FALSE)
    if (type=="one-sided") {
        return(upper.tail)
    } else {
        lower.tail <- pf(var.ratio, df1, df2, lower=TRUE)
        return(2*pmin(lower.tail, upper.tail)) # effectively Bonferroni correction.
    }
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

