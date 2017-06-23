\name{decomposeVariance}
\alias{decomposeVariance}

\title{Decompose variance to its components}
\description{Decompose the variance of the separate and premixed addition experiments to their components.}

\usage{
decomposeVariance(y, design=NULL) 
}

\arguments{
\item{y}{A DGEList object produced by \code{\link{setupSpikes}}, containing the log-ratios as a \code{ratio} field in the \code{samples} table.}
\item{design}{A design matrix containing uninteresting factors of variation. 
If not supplied, it defaults to an all-intercept matrix.}
}

\details{
This function estimates the variances among the cells in which spike-ins were added separately or premixed, using the \code{\link{estimateVariance}} function.
It tests for a significant difference between these estimates with \code{\link{testVariance}}, and computes the variance of volume addition.
The standard error of the variance of volume addition is also computed.

If the premixed estimate is greater than the separate-addition estimate, the p-value is automatically set to 1 and the variance of volume addition is set to zero.
Technically, when the variance of volume addition is zero, the premixed and separate-addition estimates are equal; 
their estimates could be adjusted to that of the combined pool of all ratios.
However, this is not performed here as it is rather counterintuitive.
}

\value{
A named list with variance estimates for the cells in which spike-ins were added separately (\code{total});
for the cells in which premixed spike-ins were added (\code{premixed});
and for volume addition (\code{volume}).
}

\seealso{
\code{\link{estimateVariance}},
\code{\link{testVariance}}
}

\author{
Aaron Lun
}

\examples{
example(setupSpikes)
decomposeVariance(y)
}

