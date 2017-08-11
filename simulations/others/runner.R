# This examines the technical variability of spike-ins in other data sets.

library(scater)
library(edgeR)
library(simpaler)
is.first <- TRUE

for (dataset in c("Wilson", "Islam", "Scialdone", "Islam (2)", "Grun", "Hashimshony", "Buettner")) {
    if (dataset=="Wilson") {
        incoming <- read.table('GSE61533_HTSEQ_count_results.tsv.gz', header=TRUE, row.names=1, colClasses=c('character', rep('integer', 96)))
        spike.in <- grepl('^ERCC', rownames(incoming))
        design <- cbind(rep(1, ncol(incoming)))

    } else if (dataset=="Islam") {
        incoming <- read.table("GSE29087_L139_expression_tab.txt.gz", 
                               colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
        spike.in <- grepl("SPIKE", rownames(incoming))
        grouping <- rep(c("ESC", "MEF", "Neg"), c(48, 44, 4))
        incoming <- incoming[,grouping!="Neg"]
        grouping <- grouping[grouping!="Neg"]
        design <- model.matrix(~grouping)

    } else if (dataset=="Scialdone") {
        incoming <- read.table("Marioni_lab_1_Jul_2015_gene_counts_table.txt", skip=1, row.names=1, colClasses=c("character", rep("integer", 96)))
        incoming <- incoming[!grepl("^__", rownames(incoming)),]
        spike.in <- grepl('^ERCC', rownames(incoming))
        design <- cbind(rep(1, ncol(incoming)))

    } else if (dataset=="Islam (2)") {
        incoming <- read.table("GSE46980_CombinedMoleculeCounts.tab.gz", skip=7, row.names=1, sep="\t", 
                               colClasses=c("character", vector("list", 6), rep("integer", 96)))
        spike.in <- grepl("_SPIKE_", rownames(incoming))
        design <- cbind(rep(1, ncol(incoming)))

    } else if (dataset=="Grun") {
        incoming <- read.table("GSE54695_data_transcript_counts.txt.gz", header=TRUE, row.names=1, sep="\t", 
                               colClasses=c("character", rep("numeric", 80), vector("list", 80), 
                                            rep("numeric", 80), vector("list", 80)))
        spike.in <- grepl("^ERCC", rownames(incoming))
        grouping <- sub("SC_([^_]+)_[0-9]+", "\\1", colnames(incoming))
        design <- model.matrix(~grouping)

    } else if (dataset=="Hashimshony") {
        incoming.1 <- read.table("GSE78779/GSE78779_CS1_manual.txt.gz", row.names=1, colClasses=c("character", rep("integer", 24)))
        incoming.2 <- read.table("GSE78779/GSE78779_CS2_manual.txt.gz", row.names=1, colClasses=c("character", rep("integer", 20)))
        stopifnot(identical(rownames(incoming.1), rownames(incoming.2)))
        incoming <- cbind(incoming.1, incoming.2)
        colnames(incoming) <- NULL

        incoming <- incoming[!grepl("^__", rownames(incoming)),]
        spike.in <- grepl('^ERCC', rownames(incoming))
        grouping <- factor(rep(1:2, c(ncol(incoming.1), ncol(incoming.2))))
        design <- model.matrix(~grouping)

    } else if (dataset=="Buettner") {
        incoming.G1 <- read.table("E-MTAB-2805/G1_singlecells_counts.txt", row.names=1, header=TRUE,
                                  colClasses=c("character", vector("list", 3), rep("integer", 96)))
        incoming.G2M <- read.table("E-MTAB-2805/G2M_singlecells_counts.txt", row.names=1, header=TRUE,
                                  colClasses=c("character", vector("list", 3), rep("integer", 96)))
        incoming.S <- read.table("E-MTAB-2805/S_singlecells_counts.txt", row.names=1, header=TRUE,
                                  colClasses=c("character", vector("list", 3), rep("integer", 96)))
        stopifnot(identical(rownames(incoming.G1), rownames(incoming.G2M)))
        stopifnot(identical(rownames(incoming.G1), rownames(incoming.S)))
        incoming <- cbind(incoming.G1, incoming.G2M, incoming.S)

        incoming <- incoming[seq_len(max(grep("^ERCC", rownames(incoming)))),]
        spike.in <- grepl("^ERCC", rownames(incoming))
        grouping <- factor(rep(c("G1", "G2M", "S"), c(ncol(incoming.G1), ncol(incoming.G2M), ncol(incoming.S))))
        design <- model.matrix(~grouping)        
        
    }

    # Quality control on the sums of the spike in counts and endogenous transcripts.
    spike.counts <- as.matrix(incoming[spike.in,])
    other.sums <- colSums(incoming) - colSums(spike.counts)
    discard <- isOutlier(colSums(spike.counts), log=TRUE, nmads=3, type="lower") |
               isOutlier(other.sums, log=TRUE, nmads=3, type="lower")

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



