###########################################################################
# Identifying highly variable genes using the QL strategy.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(matrixStats))
suppressPackageStartupMessages(require(statmod))

temp <- "temp.txt"
if (file.exists(temp))  { unlink(temp) }
filled <- FALSE

top.hits <- c(20, 200, 2000)

###########################################################################
# Running across all data types (saving diagnostics along the way)

for (datatype in c("Wilson", "Scialdone", "Kolodziejczyk", "Calero", "Liora")) {
    val <- readRDS(file.path("../datasets", paste0(datatype, ".rds")))
    incoming <- val$counts
    spike.in <- val$spikes
    design <- val$design
    use.parametric <- !datatype %in% "Kolodziejczyk" # parametric mode doesn't work very well, for some reason.

    # Quality control on cells.
    totals <- colSums(incoming)
    f <- designAsFactor(design)
    okay.libs <- !isOutlier(totals, nmad=3, log=TRUE, type="lower", batch=f) & 
                 !isOutlier(colSums(incoming!=0), nmad=3, log=TRUE, type="lower", batch=f) &
                 !isOutlier(colSums(incoming[spike.in,])/totals, nmad=3, type="higher", batch=f)
	incoming <- incoming[,okay.libs]
    design <- design[okay.libs,,drop=FALSE]

    filter.keep <- calcAverage(incoming) >= 0.1
	countsCell <- incoming[filter.keep & !spike.in,]
	spike.param <- spikeParam(incoming[filter.keep & spike.in,], design=design)
    block <- designAsFactor(design)

    #########################################################################
    # Running through our various methods.

    for (method in c("VarLog", "TechCV")) { 
        for (my.var in c(0.015)) { 
            set.seed(34271)
            results <- top.res <- sig.res <- list()
            
            for (i in seq(0,20)) {
                # Running original data (once for each dataset), then resampled versions.
                if (i) { countsSpike <- resampleSpikes(spike.param, var.log=my.var) }
                else { countsSpike <- spike.param$counts }
                
                sce <- SingleCellExperiment(list(counts=rbind(countsCell, countsSpike)))
                isSpike(sce, "Spike") <- nrow(countsCell) + seq_len(nrow(countsSpike))
                sce <- computeSpikeFactors(sce)
                sce <- normalize(sce)
                
                if (method=="VarLog"){ 
                    # Fitting the trend to the spike-in variances, and computing the biological component.
                    out <- trendVar(sce, parametric=use.parametric, design=design)
                    out2 <- decomposeVar(sce, out)
                    out2 <- out2[!isSpike(sce),]
                    my.rank <- rank(out2$p.value)
                    is.sig <- out2$FDR <= 0.05 & !is.na(out2$FDR)

                    # Some diagnostics for this data set.
        			if (i==0) { 
                        pdf(paste0("diagnostics_", datatype, ".pdf"))
                        plot(out$mean, out$var)
                        curve(out$trend(x), col="red", lwd=2, add=TRUE)
                        dev.off()
                    }
                    if (i<=1) {
                        output <- out2[is.sig,]
                        write.table(file=paste0("log_", datatype, "_", i, ".tsv"), output[order(output$bio, decreasing=TRUE),], 
                                    quote=FALSE, sep="\t", col.names=NA)
                    }

                } else {
                    # Alternatively, fitting a trend to the CV2 using the Brennecke method.
                    sce2 <- sce
                    if (!is.null(block)) {   
                        leftovers <- removeBatchEffect(exprs(sce2), block=block[okay.libs])
                        of.value <- t(t(2^leftovers - 1) * sizeFactors(sce2))
                        of.value[of.value < 0] <- 0
                        counts(sce2) <- of.value
                    }
                    outt <- technicalCV2(sce2, min.bio.disp=0)
                    outt <- outt[!isSpike(sce2),]
                    my.rank <- rank(outt$p.value)
                    is.sig <- outt$FDR <= 0.05 & !is.na(outt$FDR)

                    if (i<=1) {
                        output <- outt[is.sig,]
                        write.table(file=paste0("cv2_", datatype, "_", i, ".tsv"), output[order(output$cv2, decreasing=TRUE),], 
                                    quote=FALSE, sep="\t", col.names=NA)
                    } 
                }
                    
                top.ranked <- lapply(top.hits, function(x) { my.rank <= x })       
                if (i) { 
                    top.res[[i]] <- sapply(seq_along(top.hits), function(j) { sum(best[[j]] & !top.ranked[[j]])/top.hits[j] })
                    sig.res[[i]] <- sum(sig!=is.sig)/sum(sig)
                } else {
                    best <- top.ranked
                    sig <- is.sig
                }
            }

    		top.lost <- do.call(rbind, top.res)
            colnames(top.lost) <- paste0("Top", top.hits)
    		write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method=method, top.lost, Sig=unlist(sig.res)),
    			row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
    		filled <- TRUE
        }
	}
}

###########################################################################

file.rename(temp, "results.txt")
sessionInfo()

###########################################################################
# End.
