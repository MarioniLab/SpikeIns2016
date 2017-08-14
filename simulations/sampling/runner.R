# This examines the technical variability of spike-ins in other data sets.

library(scater)
library(edgeR)
library(simpaler)
is.first <- TRUE
set.seed(9999)

all.values <- list()
for (dataset in c("Islam2",
                  "Buettner", 
                  "Grun", 
                  "Kolodziejczyk", 
#                  "Islam", 
                  "Hashimshony", 
                  "Liora",
#                  "Wilson", 
                  "Calero", 
                  "Scialdone", 
                  "Zeisel")) {
    val <- readRDS(file.path("../datasets", paste0(dataset, ".rds")))
    incoming <- val$counts
    spike.in <- val$spikes
    design <- val$design
    name <- val$name

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
    fit <- glmFit(d, redesign, prior.count=0)
    fitted <- fit$fitted.values

    # Simulating data and calculating the variance.
    collected <- vector("list", 20)
    for (it in seq_len(20)) { 
        sim.spikes <- matrix(rnbinom(length(fitted), mu=fitted, size=1/d$tagwise.dispersion),
                             nrow=nrow(fitted), ncol=ncol(fitted))
        spike.totals <- colSums(sim.spikes)
        spike.totals <- spike.totals/mean(spike.totals)
        collected[[it]] <- spike.totals
    }

    spike.totals <- colSums(d$counts)
    spike.totals <- spike.totals/mean(spike.totals)
    collected <- do.call(cbind, collected)/spike.totals
    all.values[[name]] <- apply(collected, 1, sd) * 100
}

pdf("collated.pdf", width=10, height=6)
par(mar=c(9.1, 4.1, 4.1, 2.1))
boxplot(all.values, ylab="Size factor estimation error (%)", cex.axis=1.2, cex.names=1.2, 
        cex.lab=1.4, col="grey80", las=2)
dev.off()
