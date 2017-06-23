\name{resampleSpikes}
\alias{resampleSpikes}

\title{Resample spike-in counts}
\description{Use the estimated spike-in parameters to generate simulated spike-in counts.}

\usage{
resampleSpikes(param, var.log)
}

\arguments{
\item{param}{A list of estimates produced by \code{\link{spikeParam}}.}
\item{var.log}{A numeric scalar specifying the variance of the log2-(spike-in totals) due to technical noise.}
}

\details{
This function generates new spike-in totals for each cell based on \code{var.log}.
It does so by computing the variance of the log2-totals; subtracting \code{var.log}, to obtain the variance due to other factors;
and then adding random values from a normal distribution with variance equal to \code{var.log}.
This procedure avoids doubling the contribution of \code{var.log} (which is what would happen if the normal variates were added directly to the log-totals).

The simulated totals are used to scale up the fitted values in \code{param}.
Quantile-quantile mapping is performed to convert the old counts with the old fitted values to new counts based on the new fitted values.
This is done using the \code{\link{q2qnbinom}} function, assuming a NB distribution with the dispersion parameters listed in \code{param}.
Any negative values are set to zero in the returned count matrix.
}

\value{
A numeric matrix of simulated spike-in counts for each transcript in each cell.
} 

\author{
Aaron Lun
}

\seealso{
\code{\link{spikeParam}},
\code{\link{q2qnbinom}}
}

\examples{
example(spikeParam)
new.counts <- resampleSpikes(y, 0.01)
}

% # Testing, testing!
% a <- matrix(rnbinom(10000, mu=1000, size=2), nrow=10)
% p <- spikeParam(a)
% output <- numeric(20)
% for (i in seq_len(20)) {
%     sim <- resampleSpikes(p, 0.1)
%     output[i] <- var(log2(colSums(sim)))
% } 
%
% # Estimates should be quite similar:
% print(mean(output)) 
% var(log2(colSums(a)))

