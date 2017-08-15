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

    # Regressing out any blocking effects.
    spike.totals <- log2(colSums(spike.counts[,!discard]))
    fit <- lm.fit(y=spike.totals, x=design[!discard,,drop=FALSE])
    all.values[[name]] <- (2^abs(fit$residuals) - 1)*100
}

pdf("collated_sizefac.pdf", width=10, height=6)
par(mar=c(9.1, 4.1, 4.1, 2.1))
boxplot(all.values, ylab="Absolute deviation from mean (%)", cex.axis=1.2, cex.names=1.2, 
        cex.lab=1.4, col="grey80", las=2, ylim=c(0, 100))
dev.off()
