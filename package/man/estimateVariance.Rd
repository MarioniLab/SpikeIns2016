\name{estimateVariance}
\alias{estimateVariance}

\title{Estimates the variance}
\description{Estimates the sample variance of the log-ratios between spike-in sets.} 

\usage{
estimateVariance(y, design=NULL, ..., ratios)
}

\arguments{
\item{y}{A DGEList object produced by \code{\link{setupSpikes}}, containing the log-ratios as a \code{ratio} field in the \code{samples} table.}
\item{design}{A design matrix containing uninteresting factors of variation. 
If not supplied, it defaults to an all-intercept matrix.}
\item{...}{Other arguments to pass to \code{\link{lm.fit}}.}
\item{ratios}{A numeric vector of log-ratios. 
If supplied, it is used instead of \code{y$samples$ratio}.}
}

\details{
This function computes the variance of the log-ratios, after fitting a linear model to regress out any systematic effects on the mean.
The number of degrees of freedom is stored in the attributes, along with the standard error of the variance estimate.
}

\value{
A numeric scalar containing the variance estimate, along with the attributes \code{df} and \code{standard.error}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{setupSpikes}},
\code{\link{lm.fit}}
}

\examples{
example(setupSpikes)
estimateVariance(y[,y$samples$separate])
estimateVariance(y[,y$samples$premixed])
}

% Checking the standard error calculation under normality:
% s^2 ~ sig^2 X_k^2 / k, where X^2_k is chi-squared distributed under k degrees of freedom.
% var(s^2) = sig^4/k^2 * (2k), as the variance of X^2_k is just 2k.
% sd(s^2) = sig^2 * sqrt(2/k)
%
% Confirming with simulations:
% sig <- 2.5
% a <- matrix(rnorm(20000, sd=sig), ncol=100)
% lv <- apply(a, 1, var)
% sqrt(var(lv))  # Observed
% sig^2 * sqrt(2/(ncol(a)-1)) # Expected
