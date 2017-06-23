\name{diagnoseVariance}
\alias{diagnoseVariance}

\title{Diagnose differences in variance}
\description{Test for differences in the variance of the log-ratios between groups of cells.}

\usage{
diagnoseVariance(y, groups, design=NULL) 
}

\arguments{
\item{y}{A DGEList object produced by \code{\link{setupSpikes}}, containing the log-ratios as a \code{ratio} field in the \code{samples} table.}
\item{groups}{A factor specifying the assigned group for each cell in \code{y}.}
\item{design}{A design matrix containing uninteresting factors of variation. 
If not supplied, it defaults to an all-intercept matrix.}
}

\details{
This function will use \code{\link{estimateVariance}} to estimate the variance of the log-ratios for each group of cells in \code{groups}.
It will then perform two-sided tests between each pair of groups to determine if there is any difference in the log-ratios.
The aim is to quickly spot any differences in the variance estimates for specific groups (e.g., spike-in addition order, treated or untreated cells).

The design matrix is used to block on any uninteresting factors of variation within each group.
It will automatically be converted to full rank if the \code{group} is confounded with any factors in \code{design}.
}

\value{
A named list with two elements:
\describe{
\item{\code{var}:}{A named list of variance estimates for each group, produced by \code{\link{estimateVariance}}.}
\item{\code{pval}:}{A symmetric numeric matrix of p-values.
Each row and column corresponds to a group and the entry contains the p-value for the corresponding pairwise comparison.}
}
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
diagnoseVariance(y, y$samples$separate)
}

