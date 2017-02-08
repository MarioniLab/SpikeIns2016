###########################################################################
# Identifying highly variable genes using the QL strategy.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(matrixStats))
suppressPackageStartupMessages(require(statmod))

temp <- "temp.txt"
if (file.exists(temp))  { unlink(temp) }
filled <- FALSE

top.hits <- c(20, 200, 2000)

###########################################################################
# Running across all data types (saving diagnostics along the way)

for (datatype in c("wilson", "calero", "liora")) { 
    block <- NULL

	if (datatype=="wilson") { 
		incoming <- read.table('GSE61533_HTSEQ_count_results.tsv.gz', header=TRUE, row.names=1, colClasses=c('character', rep('integer', 96)))
		spike.in <- grepl('^ERCC', rownames(incoming))
		is.mito <- grepl("^mt-", rownames(incoming)) & !spike.in

	} else if (datatype=="calero" || datatype=="liora") { 
        if (datatype=="calero") {
            incoming.1 <- readRDS("../../real/Calero/trial_20160113/analysis/full.rds")
            incoming.2 <- readRDS("../../real/Calero/trial_20160325/analysis/full.rds")
            stopifnot(identical(rownames(incoming.1), rownames(incoming.2)))

            incoming <- cbind(incoming.1$counts, incoming.2$counts)
            gdata <- incoming.1$genes
            block <- paste0(rep(LETTERS[1:2], c(ncol(incoming.1), ncol(incoming.2))), ".",
                            ifelse(c(incoming.1$samples$induced, incoming.2$samples$induced), "Induced", "Control"))
        } else { 
            incoming.1 <- readRDS("../../real/Liora/test_20160906/analysis/full.rds")
            incoming.1 <- incoming.1[,is.na(incoming.1$samples$control.well)]
            incoming.2 <- readRDS("../../real/Liora/test_20170201/analysis/full.rds")
            incoming.2 <- incoming.2[,is.na(incoming.2$samples$control.well)]
            incoming <- cbind(incoming.1$counts, incoming.2$counts)
            gdata <- incoming.1$genes
            block <- NULL
            block <- rep(LETTERS[1:2], c(ncol(incoming.1), ncol(incoming.2)))
        }

        # Getting rid of SIRVS (spike2), and setting ERCCs (spike1) as the spike-ins.
        keep <- !gdata$spike2 
        incoming <- incoming[keep,]
        gdata <- gdata[keep,]
        spike.in <- gdata$spike1
        is.mito <- gdata$is.mito
    }

    # Quality control on cells.
    totals <- colSums(incoming)
    okay.libs <- !isOutlier(totals, nmad=3, log=TRUE, type="lower") & 
                 !isOutlier(colSums(incoming!=0), nmad=3, log=TRUE, type="lower") &
                 !isOutlier(colSums(incoming[is.mito,])/totals, nmad=3, type="higher") & 
                 !isOutlier(colSums(incoming[spike.in,])/totals, nmad=3, type="higher")
	incoming <- incoming[,okay.libs]

    # Filtering out crappy spikes.
    high.ab <- rowMeans(incoming) >= 1
	countsCell <- as.matrix(incoming[high.ab & !spike.in,])
	spike.param <- spikeParam(incoming[high.ab & spike.in,])
	diag.done <- FALSE

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
                
                sce <- newSCESet(countData=rbind(countsCell, countsSpike))
                sce <- calculateQCMetrics(sce, feature_controls=list(Spike=rep(c(FALSE, TRUE), c(nrow(countsCell), nrow(countsSpike)))))
                isSpike(sce) <- "Spike"
                sce <- computeSpikeFactors(sce)
                sce <- normalize(sce)
                
                if (method=="VarLog"){ 
                    if (is.null(block)) { design <- NULL }
                    else { design <- model.matrix(~block[okay.libs]) }

                    # Fitting the trend to the spike-in variances.
                    out <- trendVar(sce, trend="semiloess", design=design)

                    # Computing the biological component.
                    out2 <- decomposeVar(sce, out)
                    out2 <- out2[!isSpike(sce),]
                    my.rank <- rank(out2$p.value)
                    is.sig <- out2$FDR <= 0.05 

                    # Some diagnostics for this data set.
        			if (!diag.done) { 
                        pdf(paste0("diagnostics_", datatype, ".pdf"))
                        plot(out$mean, out$var)
                        curve(out$trend(x), col="red", lwd=2, add=TRUE)
                        dev.off()

                        output <- out2[is.sig,]
                        write.table(file=paste0("out_", datatype, ".tsv"), output[order(output$bio, decreasing=TRUE),], 
                                    quote=FALSE, sep="\t", col.names=NA)
                        diag.done <- TRUE
                    }
                } else {
                    sce2 <- sce
                    if (!is.null(block)) {   
                        leftovers <- removeBatchEffect(exprs(sce2), block=block[okay.libs])
                        of.value <- t(t(2^leftovers - sce2@logExprsOffset) * sizeFactors(sce2))
                        of.value[of.value < 0] <- 0
                        counts(sce2) <- of.value
                    }
                    outt <- technicalCV2(sce2, min.bio.disp=0)
                    outt <- outt[!isSpike(sce2),]
                    my.rank <- rank(outt$p.value)
                    is.sig <- outt$FDR <= 0.05
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
