\name{spikeParam}
\alias{spikeParam}

\title{Sets up spike-in parameters}
\description{Estimates spike-in parameters for simulation of spike-in counts.}

\usage{
spikeParam(spike.counts, design) 
}

\arguments{
\item{spike.counts}{A count matrix containing only rows for spike-in transcripts.}
\item{design}{A design matrix containing uninteresting factors of variation. If not supplied, it defaults to an all-intercept matrix.}
}

\details{
This function estimates the negative binomial dispersion for each row.
It fits a one-way NB GLM to obtain fitted values for each transcript.
The idea is to use these parameters in \code{\link{resampleSpikes}} for simulating new spike-in data.

It is worth mentioning that a related function is \code{\link{setupSpikes}}.
However, that function is used specifically for analyses of real data from the mixture experiments.
\code{spikeParam} is intended for use in estimating parameters for simulation studies.
}

\value{
A named list with components:
\describe{
\item{counts}{The count matrix for a set of spike-in transcripts.}
\item{totals}{A numeric vector of total counts for the spike-in set.}
\item{fitted}{The fitted values of the NB GLM for each transcript in each cell.}
\item{dispersion}{A numeric vector of NB dispersion estimates for each transcript.}
}
} 

\author{
Aaron Lun
}

\seealso{
\code{\link{resampleSpikes}}
}

\examples{
x <- matrix(rpois(100000, lambda=5), ncol=100)
y <- spikeParam(x)
}
