testVariance <- function(var1, var2, type=c("one-sided", "two-sided")) 
# Tests if the first variance is significantly larger than the second.
# If two-sided, it tests whether there are any significance differences between the variances.
#
# written by Aaron Lun
# created 26 January 2016
# last modified 23 June 2017
{
    df1 <- attributes(var1)$df
    df2 <- attributes(var2)$df

    type <- match.arg(type)
    var.ratio <- var1/var2
    attributes(var.ratio) <- NULL

    upper.tail <- pf(var.ratio, df1, df2, lower.tail=FALSE)
    if (type=="one-sided") {
        return(upper.tail)
    } else {
        lower.tail <- pf(var.ratio, df1, df2, lower.tail=TRUE)
        return(2*pmin(lower.tail, upper.tail)) # effectively Bonferroni correction.
    }
}

