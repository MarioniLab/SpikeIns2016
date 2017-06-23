estimateVariance <- function(y, design=NULL, ..., ratios)
# Estimates the variance of the ratios, given the design matrix.
#
# written by Aaron Lun
# created 26 January 2016
# last modified 23 June 2017
{
    if (missing(ratios)) { ratios <- y$samples$ratio }
    if (is.null(design)) { design <- .make_intercept(length(ratios)) }
    else { design <- .restore_rank(design) }
    fit <- lm.fit(y=ratios, x=design, ...)

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

