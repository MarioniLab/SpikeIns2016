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

	if (datatype=="wilson") { 
		incoming <- read.table('GSE61533_HTSEQ_count_results.tsv.gz', header=TRUE, row.names=1, colClasses=c('character', rep('integer', 96)))
		spike.in <- grepl('^ERCC', rownames(incoming))

		# Getting metrics.
        totals <- colSums(incoming)
		is.mito <- grepl("^mt-", rownames(incoming)) & !spike.in
	} else if (datatype=="calero" || datatype=="liora") { 
        if (datatype=="calero") {
            incoming.1 <- read.table("../../real/Calero/trial_20160113/analysis/genic_counts.tsv", header=TRUE, row.names=1, colClasses=c('character', rep('integer', 97)))
            incoming.1 <- incoming.1[,-1]
            incoming.2 <- read.table("../../real/Calero/trial_20160325/analysis/genic_counts.tsv", header=TRUE, row.names=1, colClasses=c('character', rep('integer', 97)))
            incoming.2 <- incoming.2[,-1]
            incoming <- cbind(incoming.1, incoming.2)
        } else { 
            incoming <- read.table("../../real/Liora/test_20160906/analysis/genic_counts.tsv", header=TRUE, row.names=1, colClasses=c('character', rep('integer', 98)))
            incoming <- incoming[,-1]
            metadata <- read.table("../../real/Liora/test_20160906/analysis/wellID.tsv", header=TRUE)
            incoming <- incoming[,metadata$Colnames[! metadata$Well %in% c("A01", "H12") | is.na(metadata$Well)]]
        }
        incoming <- incoming[!grepl("SIRV", rownames(incoming)),] # Getting rid of SIRVs at this point.
        spike.in <- grepl("^ERCC", rownames(incoming))

		# Getting metrics.
		totals <- colSums(incoming)
        library(TxDb.Mmusculus.UCSC.mm10.ensGene)
        chr.loc <- select(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rownames(incoming), keytype="GENEID", column="CDSCHROM")
        chr.loc <- chr.loc$CDSCHROM[match(rownames(incoming), chr.loc$GENEID)]
		is.mito <- chr.loc=="chrM" & !is.na(chr.loc)
    }

    # Quality control on cells.
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
        for (my.var in c(0.01)) { 
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
                    # Fitting the trend to the spike-in variances.
                    out <- trendVar(sce, trend="semiloess")

                    # Computing the biological component.
                    out2 <- decomposeVar(sce, out)
                    out2 <- out2[!isSpike(sce),]
                    my.rank <- rank(out2$p.value)
                    is.sig <- out2$FDR <= 0.05 & out2$bio > 0.5

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
                    outt <- technicalCV2(sce, min.bio.disp=0)
                    outt <- outt[!isSpike(sce),]
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
