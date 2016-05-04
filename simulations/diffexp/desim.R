#################################################################################
# This performs a differential expression analysis of various data sets, 
# exploiting the presence of spike-ins for normalization.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(monocle))
suppressPackageStartupMessages(library(MAST))

# Note that Kharchenko's SCDE package doesn't seem to support spike-in normalization.

fdr.threshold <- 0.05
top.hits <- 200

temp <- "temp.txt"
if (file.exists(temp))  { unlink(temp) }
filled <- FALSE

#################################################################################

for (datatype in c("kolod", "islam")) {
	if (datatype=="islam") {
		counts <- read.table("GSE29087_L139_expression_tab.txt.gz", 
			colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
		is.spike <- grepl("SPIKE", rownames(counts))
		grouping <- factor(c(rep(c("ESC", "MEF", "Neg"), c(48, 44, 4))))
		
		# Quality control on individual cells.
		totals <- colSums(counts)
		is.mito <- grepl("^mt-", rownames(counts)) & !is.spike
		okay.libs <- !isOutlier(totals, n=3, type="lower", log=TRUE) &
                     !isOutlier(colSums(counts[is.mito,])/totals, n=3, type="higher") & 
                     !isOutlier(colSums(counts[is.spike,])/totals, n=3, type="higher")
		counts <- counts[,okay.libs]
		grouping <- droplevels(grouping[okay.libs])
		
	} else if (datatype=="kolod") {
        all.counts <- read.table("ESpresso/counttable_es.csv", header=TRUE, row.names=1, colClasses=c("character", rep("integer", 704)))
        serum <- sub("ola_mES_([^_]+)_.*", "\\1", colnames(all.counts))
        batch <- sub("ola_mES_[^_]+_([^_]+)_.*", "\\1", colnames(all.counts))
        targets <- data.frame(Serum=serum, Batch=batch)
                                            
        # Only using data from the batch with spike-ins.
        keep <- targets$Batch=="3" & targets$Serum %in% c("lif", "2i")
        counts <- all.counts[,keep]
        is.spike <- grepl("^ERCC", rownames(counts))
        is.mouse <- grepl("^ENSMUSG", rownames(counts))

        # Getting rid of mapping statistics stored in the same file.
        counts <- counts[is.spike|is.mouse,] 
        is.spike <- is.spike[is.spike|is.mouse]
        grouping <- factor(targets$Serum[keep])
	}

	design <- model.matrix(~grouping)
	filter.keep <- rowMeans(counts) >= 1
    spike.param <- spikeParam(counts[is.spike & filter.keep,])

	#################################################################################
	# Running edgeR first.

	# Filtering out crappy genes.
	y.ref <- DGEList(counts[!is.spike & filter.keep,])
	diag.done <- FALSE

    for (my.var in c(0.001, 0.01, 0.1)) {
		set.seed(231234)
		results <- top.set <- list()
		
		for (i in seq(0,20)) {
			# Running original data, then resampled versions.
			if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
			else { spike.data <- spike.param$counts }
		
			# Normalizing on the spike-in counts.
			y <- y.ref
			y$samples$norm.factors <- colSums(spike.data)/y$samples$lib.size

			# Running edgeR.
			y <- estimateDisp(y, design)
			fit <- glmFit(y, design, robust=TRUE)
			res <- glmLRT(fit)

			# Check if diagnostics look okay:
			if (!diag.done) {
				pdf(paste0("diagnostics_", datatype, ".pdf"))
				plotBCV(y) # Nice downward trend
				limma::plotMDS(edgeR::cpm(y, log=TRUE), col=c("red", "blue")[(as.integer(grouping)==1)+1L]) # Mostly separated
				dev.off()
			}

			# Comparing.
			chosen <- p.adjust(res$table$PValue, method="BH") <= fdr.threshold
		    top.ranked <- rank(res$table$PValue) <= top.hits
			if (i) { 
				results[[i]] <- sum(original!=chosen)/sum(original)
				top.set[[i]] <- sum(best!=top.ranked)/top.hits
			} else {
				original <- chosen
				best <- top.ranked
			}

			if (!diag.done) {
				# Top genes also look like the ones listed in Islam's paper, which is good, e.g. Sparc, S100a6, Vim, Fn1.
				de.out <- topTags(res, n=Inf)
				write.table(file=paste0("outEB_", datatype, ".tsv"), de.out, row.names=TRUE, 
					quote=FALSE, sep="\t", col.names=NA)
				diag.done <- TRUE
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- "Total"
		top.lost <- do.call(rbind, top.set)
        colnames(top.lost) <- paste0("Top", top.hits)

        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="edgeR", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}

	#################################################################################
	# Next, monocle.

	for (my.var in c(0.001, 0.01, 0.1)) {
		set.seed(3742)
		results <- top.set <- list()
		
		for (i in seq(0,20)) {
			# Running original data, then resampled versions.
			if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
			else { spike.data <- spike.param$counts }

			# Normalizing by spike-in total counts (scaled to the total library size).
			spike.totals <- colSums(spike.data)
			keep <- !is.spike & filter.keep
            spike.totals <- spike.totals/mean(spike.totals) * mean(colSums(counts[keep,]))
			normalized <- edgeR::cpm(counts[keep,], lib.size=spike.totals)

			pdat <- AnnotatedDataFrame(data=data.frame(grouping=grouping))
			sampleNames(pdat) <- colnames(normalized)
			HSMM <- newCellDataSet(cellData=normalized, phenoData=pdat)

            HSMM <- HSMM[1:2000,] # for convenience, otherwise we'll be here for days. 
			out <- differentialGeneTest(HSMM, fullModelFormulaStr="expression~grouping", cores=10) 

			# Choosing stuff. 
			chosen <- p.adjust(out$pval, method="BH") <= fdr.threshold
		    top.ranked <- rank(out$pval) <= top.hits
			if (i) { 
				results[[i]] <- sum(original!=chosen)/sum(original)
				top.set[[i]] <- sum(best!=top.ranked)/top.hits
			} else {
				original <- chosen
				best <- top.ranked
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- "Total"
		top.lost <- do.call(rbind, top.set)
        colnames(top.lost) <- paste0("Top", top.hits)

        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="monocle", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
    }

    #################################################################################
    # Finally, MAST.

    for (my.var in c(0.001, 0.01, 0.1)) {
        set.seed(3742)
        results <- top.set <- list()
    
        for (i in seq(0,20)) {
            # Running original data, then resampled versions.
            if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
            else { spike.data <- spike.param$counts }

            # Normalizing by spike-in total counts (again, scaled to the total library size).
            spike.totals <- colSums(spike.data)
            keep <- !is.spike & filter.keep
            spike.totals <- spike.totals/mean(spike.totals) * mean(colSums(counts[keep,]))
            lcpms <- log2(t(counts[keep,]+1)/spike.totals*1e6)

            suppressMessages({
                sca <- FromMatrix('SingleCellAssay', lcpms, data.frame(wellKey=rownames(lcpms)), data.frame(primerid=colnames(lcpms)))
                cData(sca) <- cbind(cData(sca), grouping)
                cData(sca)$cngeneson <- colMeans(counts>0)
                                                          
                mfit <- zlm.SingleCellAssay(~ grouping + cngeneson, sca, method="bayesglm", ebayes=TRUE, ebayesControl=list(method="MLE", model="H1"))
                mlrt <- lrTest(mfit, Hypothesis(colnames(design)[ncol(design)], colnames(coef(mfit, "D")))) # cngeneson is added at the end, so it shouldn't affect terms before it.
                MAST.p <- mlrt[, "hurdle", "Pr(>Chisq)"]
            })                                                                            
        
            # Choosing stuff. 
			chosen <- p.adjust(MAST.p, method="BH") <= fdr.threshold
		    top.ranked <- rank(MAST.p) <= top.hits
			if (i) { 
				results[[i]] <- sum(original!=chosen)/sum(original)
				top.set[[i]] <- sum(best!=top.ranked)/top.hits
			} else {
				original <- chosen
				best <- top.ranked
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- "Total"
		top.lost <- do.call(rbind, top.set)
        colnames(top.lost) <- paste0("Top", top.hits)

        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="MAST", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
    }
}

#################################################################################

file.rename(temp, "results.txt")
sessionInfo()

###########################################################################

