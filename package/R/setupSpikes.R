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

