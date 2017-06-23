\name{testVariance}
\alias{testVariance}

\title{Test for differences in variances.}
\description{Test for differences in the variance estimates for the log-ratios.}

\usage{
testVariance(var1, var2, type=c("one-sided", "two-sided")) 
}

\arguments{
\item{var1, var2}{Numeric scalars with attributes, containing variance estimates and associated information. 
These should be the output of \code{\link{estimateVariance}}.}
\item{type}{A string specifying the type of test to perform.}
}

\details{
This function will apply a F-test to test for equality of the true values of \code{var1} and \code{var2}.
This assumes that the two values are independent, e.g., computed from separate observations.
By default, a one-sided test is performed against the null hypothesis that \code{var1} is less than or equal to \code{var2}.
If \code{type="two-sided"}, a two-sided test will be performed for any difference between \code{var1} and \code{var2}.
}

\value{
A numeric scalar containing the p-value.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{estimateVariance}}
}

\examples{
example(setupSpikes)
v1 <- estimateVariance(y[,y$samples$separate])
v2 <- estimateVariance(y[,y$samples$premixed])
testVariance(v1, v2)
}

% # Checking that the test is correct.
% ngenes <- 1000
% a <- matrix(rnorm(ngenes*40), ncol=40)
% b <- matrix(rnorm(ngenes*100), ncol=100)
% output <- numeric(ngenes)
% for (x in seq_len(ngenes)) {
%     v1 <- estimateVariance(ratios=a[x,])
%     v2 <- estimateVariance(ratios=b[x,])
%     output[x] <- testVariance(v1, v2)
% }
% hist(output) # should be uniform. 

