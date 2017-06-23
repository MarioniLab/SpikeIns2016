\name{setupSpikes}
\alias{setupSpikes}

\title{Set up spike-in information}
\description{Extract information about spike-in transcripts from a scRNA-seq count matrix.}

\usage{
setupSpikes(counts, spike1, spike2, separate, premixed, ...) 
} 

\arguments{
\item{counts}{A numeric count matrix, one row per genomic feature (including spike-in transcripts) and one column per cell.}
\item{spike1, spike2}{Subset vectors indicating the rows of \code{counts} corresponding to the first or second spike-in sets.}
\item{separate}{A subset vector specifying the cells for which the two spike-in sets were added separately.}
\item{premixed}{A subset vector specifying the cells for which the two spike-in sets were premixed before addition.}
\item{...}{Other cell-specific metadata fields to be stored.}
}

\details{
For each cell, this function will compute the total coverage of each spike-in set and the log-ratio of the total counts between sets.
The \code{separate} and \code{premixed} vectors are converted to logical vectors and stored as cell metadata.
The identities of the spike-in transcripts (\code{spike1} and \code{spike2}) are also stored as logical vectors in the gene metadata.
}

\value{
A DGEList object containing the count matrix and the metadata for each gene and cell.
}

\author{
Aaron Lun
}

\examples{
x <- matrix(rpois(100000, lambda=5), ncol=100)
spike1 <- 1:10
spike2 <- 11:20
separate <- 1:50
premixed <- 51:100

y <- setupSpikes(x, spike1, spike2, separate, premixed)
}
