diagnoseVariance <- function(y, groups, design=NULL) 
# This estimates the variance for each level of 'groups'.
# It then does pairwise comparisons between groups to check for significant differences.
#
# written by Aaron Lun
# created 26 January 2016
# last modified 23 June 2017
{
    if (is.null(design)) { design <- .make_intercept(length(groups)) }
    ratios <- y$samples$ratio
    all.index <- split(seq_len(length(ratios)), groups)
    all.var <- vector("list", length(all.index))
    names(all.var) <- names(all.index)

    for (g in names(all.index)) {
        current <- all.index[[g]]
        cur.ratios <- ratios[current]
        cur.design <- .restore_rank(design[current,,drop=FALSE])
        all.var[[g]] <- estimateVariance(ratios=cur.ratios, design=cur.design)
    }
    
    pairwise <- matrix(NA, length(all.var), length(all.var))
    for (g1 in seq_along(all.var)) { 
        for (g2 in seq_len(g1-1L)) {
            out <- testVariance(var1=all.var[[g1]], var2=all.var[[g2]], type="two-sided")
            pairwise[g1,g2] <- out
            pairwise[g2,g1] <- out
        }
    }
    rownames(pairwise) <- colnames(pairwise) <- names(all.var)
    return(list(var=all.var, pval=pairwise))
}

