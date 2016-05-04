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

for (datatype in c("wilson", "islam", "brennecke")) { 

	if (datatype=="wilson") { 
		incoming <- read.table('GSE61533_HTSEQ_count_results.tsv.gz', header=TRUE, row.names=1, colClasses=c('character', rep('integer', 96)))
		spike.in <- grepl('^ERCC', rownames(incoming))

		# Quality control on cells.
		totals <- colSums(incoming)
		is.mito <- grepl("^mt-", rownames(incoming)) & !spike.in
		okay.libs <- !isOutlier(totals, n=3, log=TRUE, type="lower") & 
                     !isOutlier(colSums(incoming[is.mito,])/totals, n=3, type="higher") & 
                     !isOutlier(colSums(incoming[spike.in,])/totals, n=3, type="higher")
		incoming <- incoming[,okay.libs]

	} else if (datatype=="islam") { 
		incoming <- read.table('GSE46980_CombinedMoleculeCounts.tab.gz', skip=7, row.names=1, sep="\t", fill=TRUE,
		        colClasses=c(list('character', NULL, NULL, NULL, NULL, NULL, NULL), rep('integer', 96)))
		spike.in <- grepl('SPIKE', rownames(incoming))

		# Quality control on cells.
		totals <- colSums(incoming[!spike.in,])
		is.mito <- grepl("^mt-", rownames(incoming)) & !spike.in
		okay.libs <- !isOutlier(totals, n=3, log=TRUE, type="lower") & 
                     !isOutlier(colSums(incoming[is.mito,])/totals, n=3, type="higher") & 
                     !isOutlier(colSums(incoming[spike.in,])/totals, n=3, type="higher")
		incoming <- incoming[,okay.libs]

	} else if (datatype=="brennecke") {
        # Note that cell-level quality control has already been performed.
        incoming <- read.table('nmeth.2645-S7.csv.gz', header=TRUE, row.names=1, sep=',', colClasses=c('character', rep('integer', 92)))
        gene.length <- as.vector(incoming[,1])
        incoming <- incoming[,-1]
        spike.in <- grepl('^ERCC', rownames(incoming))

    }

    # Filtering out crappy spikes.
    high.ab <- rowMeans(incoming) >= 1
	countsCell <- as.matrix(incoming[high.ab & !spike.in,])
	spike.param <- spikeParam(incoming[high.ab & spike.in,])
	diag.done <- FALSE

    #########################################################################
    # Using our custom method.

	for (my.var in c(0.001, 0.01, 0.1)) { 
		set.seed(34271)
		results <- top.res <- list()

		for (i in seq(0,20)) {
			# Running original data (once for each dataset), then resampled versions.
			if (i) { countsSpike <- resampleSpikes(spike.param, var.log=my.var) }
			else { countsSpike <- spike.param$counts }

            sce <- newSCESet(countData=rbind(countsCell, countsSpike))
            isSpike(sce) <- rep(c(FALSE, TRUE), c(nrow(countsCell), nrow(countsSpike)))
            sce <- computeSpikeFactors(sce)
            sce <- normalize(sce)

            # Fitting the trend to the spike-in variances.
            out <- trendVar(sce, trend="loess")
			if (!diag.done) { 
				pdf(paste0("diagnostics_", datatype, ".pdf"))
                plot(out$mean, out$var)
                sort.ab <- sort(out$mean)
                lines(sort.ab, out$trend(sort.ab), col="red", lwd=2)
				dev.off()
			}

            # Computing the biological component.
            out2 <- decomposeVar(sce, out)
            my.rank <- rank(-out2$bio)
            top.ranked <- lapply(top.hits, function(x) { my.rank <= x })            
			if (i) { 
				top.res[[i]] <- sapply(seq_along(top.hits), function(j) { sum(best[[j]] & !top.ranked[[j]]) })
			} else {
				best <- top.ranked
			}
		
			if (!diag.done) { 
				output <- out2
                colnames(output) <- c("Mean", "Total", "Bio", "Tech")
				write.table(file=paste0("out_custom_", datatype, ".tsv"), output[order(output$Bio, decreasing=TRUE),], 
                            quote=FALSE, sep="\t", col.names=NA)
				diag.done <- TRUE
			}
		}
	
		top.lost <- do.call(rbind, top.res)
        colnames(top.lost) <- paste0("Top", top.hits)
		write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="VarLog", top.lost), 
			row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}

    #########################################################################
    # Using Brennecke's method on the original counts.

    for (my.var in c(0.001, 0.01, 0.1)) { 
		set.seed(3427)
		results <- top.res <- list()

        for (i in seq(0,20)) {
			# Running original data (once for each dataset), then resampled versions.
			if (i) { countsSpike <- resampleSpikes(spike.param, var.log=my.var) }
			else { countsSpike <- spike.param$counts }

            spike.sums <- colSums(countsSpike)
            spike.sums <- spike.sums/exp(mean(log(spike.sums)))
            sfSpike <- sfCell <- spike.sums

#           # Equivalent to original normalization strategy
#           sfCell <- estimateSizeFactorsForMatrix(countsCell)
#           sfSpike <- estimateSizeFactorsForMatrix(countsSpike)

            nCountsSpike <- t(t(countsSpike)/sfSpike)
            nCountsCell <- t(t(countsCell)/sfCell)

            # Copied straight from the supplementary materials.
            meansSpike <- rowMeans( nCountsSpike )
            varsSpike <- rowVars( nCountsSpike )
            cv2Spike <- varsSpike / meansSpike^2
            
            meansCell <- rowMeans( nCountsCell )
            varsCell <- rowVars( nCountsCell )
            cv2Cell <- varsCell / meansCell^2

            minMeanForFitA <- unname( quantile( meansSpike[ which( cv2Spike > .3 ) ], .8 ) )
            useForFitA <- meansSpike >= minMeanForFitA
            fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansSpike[useForFitA] ), cv2Spike[useForFitA] )

            minBiolDisp <- .5^2
            xi <- mean( 1 / sfSpike )
            m <- ncol(countsCell)
            psia1thetaA <- mean( 1 / sfSpike ) + ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfSpike / sfCell )
            cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
            testDenomA <- ( meansCell * psia1thetaA + meansCell^2 * cv2thA ) / ( 1 + cv2thA/m )
#            pA <- 1 - pchisq( varsCell * (m-1) / testDenomA, m-1 )
#            padjA <- p.adjust( pA, "BH" )
            lpA <- pchisq( varsCell * (m-1) / testDenomA, m-1, lower=FALSE, log=TRUE)
   
            # Getting the top-ranked hits.    
            my.rank <- rank(lpA)
            top.ranked <- lapply(top.hits, function(x) { my.rank <= x })            
			if (i) { 
				top.res[[i]] <- sapply(seq_along(top.hits), function(j) { sum(best[[j]] & !top.ranked[[j]]) })
			} else {
				best <- top.ranked
			}
		
			if (!diag.done) { 
				output <- data.frame(GeneID=rownames(countsCell), CV2=cv2Cell, logPValue=lpA)
				write.table(file=paste0("out_brennecke_", datatype, ".tsv"), output[order(output$Bio, decreasing=TRUE),], 
                            quote=FALSE, sep="\t", col.names=NA)
				diag.done <- TRUE
			}
		}
	
		top.lost <- do.call(rbind, top.res)
        colnames(top.lost) <- paste0("Top", top.hits)
		write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="CV2", top.lost), 
			row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}
}

###########################################################################

file.rename(temp, "results.txt")
sessionInfo()

###########################################################################
# End.
