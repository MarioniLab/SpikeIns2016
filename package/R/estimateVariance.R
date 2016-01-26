setupSpikes <- function(counts, spike1, spike2, separate, premixed) 
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
    y$samples$separate <- .indexToLogical(separate)
    y$samples$premixed <- .indexToLogical(premixed)
    spike.data <- data.frame(spike1=.indexToLogical(spike1, byrow=TRUE), 
                             spike2=.indexToLogical(spike2, byrow=TRUE))
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

estimateVariance <- function(ratios, design, ...) 
# Estimates the variance of the ratios, given the design matrix.
#
# written by Aaron Lun
# created 26 January 2016
{
    if (missing(design)) { design <- .make_intercept(length(ratios)) }
    fit <- lm.fit(y=ratios, x=design, ...)
    sum(out$effects[-seq_len(out$rank)]^2)/out$df.residual
}

.make_intercept <- function(n) { as.matrix(rep(1, n)) }

testVariance <- function(var1, var2, df1, df2, ratios1, ratios2, design1, design2, type=c("one-sided", "two-sided"), ...) 
# Tests if the first variance is significantly larger than the second.
# If two-sided, it tests whether there are any significance differences between the variances.
#
# written by Aaron Lun
# created 26 January 2016
{
    if (missing(design1)) { design1 <- .make_intercept(length(ratios1)) }
    if (missing(var1)) { var1 <- estimateVariance(ratios1, design1, ...) } 
    if (missing(design2)) { design2 <- .make_intercept(length(ratios2)) }
    if (missing(var2)) { var2 <- estimateVariance(ratios2, design2, ...) }

    if (missing(df1)) { df1 <- nrow(design1) - ncol(design1) }
    if (missing(df2)) { df2 <- nrow(design2) - ncol(design2) }
    type <- match.arg(type)
    var.ratio <- var1/var2

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

}

diagnoseVariance <- function(y, groups, design) 
# This estimates the variance for each level of 'groups'.
# It then does pairwise comparisons between groups to check for significant differences.
#
# written by Aaron Lun
# created 26 January 2016
{
    if (missing(design)) { design <- .make_intercept(length(groups)) }
    ratios <- y$samples$ratios
    all.index <- split(seq_len(length(ratios)), groups)
    all.var <- list()
    all.df <- list()

    for (g in names(all.index)) {
        current <- all.index[[g]]
        cur.ratios <- ratios[current]
        cur.design <- design[current,,drop=FALSE]
        QR <- qr(cur.design)
        cur.design <- cur.design[,QR$pivot[seq_len(QR$rank)]]
        all.var[[g]] <- estimateVariance(cur.ratios, cur.design)
        all.df[[g]] <- nrow(cur.design) - ncol(cur.design)
    }
    
    pairwise <- matrix(NA, length(all.var), length(all.var))
    for (g1 in seq_along(all.var)) { 
        for (g2 in seq_len(g1-1L)) {
            out <- testVariance(var1=all.var[[g1]], var2=all.var[[g2]],
                                df1=all.df[[g1]], df2=all.df[[g2]])
            pairwise[g1,g2] <- out
            pairwise[g2,g1] <- out
        }
    }
    colnames(pairwise) <- paste0("vs.", names(all.var))
    return(data.frame(group=groups, var=unlist(all.var), pairwise))
}

decomposeVariance <- function(y, separate.design, premixed.design) 
# This decomposes the variances to their relevant components, given the count matrix in 'counts'.
# We need an indication of which genes are in each spike-in set (spike1, spike2).
# We also need to know which samples correspond to premixed and separate additions.
{
    ratios <- y$samples$ratio
    total.var <- estimateVariance(ratios[separate], separate.design)
    premix.var <- estimateVariance(ratios[premixed], premixed.design)
    volume.var <- 0.5*(total.var - premix.var)

    # Following the REML definition when the subtraction is negative.
    combined.design <- rbind(do.call(cbind, c(list(separate.design), rep(list(0), ncol(premixed.design)))),
                             do.call(cbind, c(list(premixed.design), rep(list(0), ncol(separate.design)))))
    if (volume.var < 0) {
        premix.var <- total.var <- estimateVariance(c(ratios[separate], ratios[premixed]), combined.design)
        volume.var <- 0
        volume.sig <- 1
    } else {
        volume.sig <- testVariance(var1=total.var, var2=premix.var, design1=separate.design, 
                                   design2=premixed.design, type="one-sided")
    }

    # Computing the variance in behaviour.
    for (s in seq_len(2)) { 
        if (s==1L) {
            spike <- y$genes$spike1
            alt.sum <- y$samples$sum2
            final.value <- "tech.var1"
        } else {
            spike <- y$genes$spike2
            alt.sum <- y$samples$sum1
            final.value <- "tech.var2"
        }
        
        splitted <- splitSpikes(y$counts[spike,,drop=FALSE])
        subsumA <- splitted$first
        subsumB <- splitted$second
        self.ratios <- log2(subsumA/subsumB) 

        # Should not have behavioural variability between spike-ins from the same population.
        self.tech.var <- estimateVariance(c(self.ratios[separate], self.ratios[premixed]), combined.design)
        against.ratiosA <- log2(subsumA/alt.sum)
        premix.varA <- estimateVariance(against.ratiosA[premixed], premixed.desgn)
        against.ratiosB <- log2(subsumB/alt.sum)
        premix.varB <- estimateVariance(against.ratiosB[premixed], premixed.desgn)

        tech.var <- premix.var - 0.5*(premix.varA + premix.varB - self.tech.var)
        tech.var <- pmax(tech.var, 0) # slightly biased, but avoid silly results.
        assign(final.value, tech.var)
    }
    behave.var <- premix.var - tech.var1 - tech.var2
    
    # Reporting the final breakdown of the variance components.
    return(list(direct=data.frame(Total=total.var, Premixed=premix.var, Volume=volume.var, Significance=volume.sig),
                split=data.frame(Tech1=tech.var1, Tech2=tech.var2, Behaviour=behave.var)))
}
