# This examines the technical variability of spike-ins in other data sets.

library(scater)
library(edgeR)
library(simpaler)
is.first <- TRUE

for (dataset in c("Wilson", "Islam", "Scialdone", "Grun", "Kolodziejczyk", "Buettner", "Calero", "Liora")) {
    val <- readRDS(file.path("../datasets", paste0(dataset, ".rds")))
    incoming <- val$counts
    spike.in <- val$spikes
    design <- val$design

    # Quality control on the sums of the spike in counts and endogenous transcripts.
    spike.counts <- as.matrix(incoming[spike.in,])
    other.sums <- colSums(incoming) - colSums(spike.counts)
    f <- designAsFactor(design)
    discard <- isOutlier(colSums(spike.counts), log=TRUE, nmads=3, type="lower", batch=f) |
               isOutlier(other.sums, log=TRUE, nmads=3, type="lower", batch=f)

    # Estimating the mean and dispersion.
    d <- DGEList(spike.counts[,!discard])
    redesign <- design[!discard,,drop=FALSE]
    d <- estimateDisp(d, redesign, prior.df=0, trend.method="none")
    fit <- glmFit(d, redesign)

    # Calculating new fitted values with a common offset for all cells.
    fitted <- exp(fit$unshrunk.coefficients %*% t(redesign) + mean(fit$offset))
    fitted[is.na(fitted)] <- 0

    # Simulating data and calculating the variance.
    collected <- vector("list", 20)
    for (it in seq_len(20)) { 
        sim.spikes <- matrix(rnbinom(length(fitted), mu=fitted, size=1/d$tagwise.dispersion),
                             nrow=nrow(fitted), ncol=ncol(fitted))
        collected[[it]] <- estimateVariance(ratios=log2(colSums(sim.spikes)), design=redesign)
    }
    collected <- unlist(collected)
    original <- estimateVariance(ratios=log2(colSums(d$counts)), design=redesign)

    # Saving to file.    
    write.table(file="collated.txt", append=!is.first, row.names=FALSE, col.names=is.first, sep="\t", quote=FALSE,
                data.frame(Dataset=dataset, Original=original, Mean=mean(unlist(collected)), 
                           SE=sqrt(var(collected)/length(collected))))
    is.first <- FALSE 
}



