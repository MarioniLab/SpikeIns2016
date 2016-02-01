#################################################################################
# This performs a differential expression analysis of various data sets, 
# exploiting the presence of spike-ins for normalization.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(SAMstrt))
suppressPackageStartupMessages(require(monocle))

# Note that Kharchenko's SCDE package doesn't seem to support spike-in normalization.

fdr.threshold <- 0.05
top.hits <- c(20, 200, 2000)

temp <- "temp.txt"
if (file.exists(temp))  { unlink(temp) }
filled <- FALSE

#################################################################################

for (datatype in c("brennecke", "islam")) {
	if (datatype=="islam") {
		counts <- read.table("GSE29087_L139_expression_tab.txt.gz", 
			colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
		is.spike <- grepl("SPIKE", rownames(counts))
		grouping <- factor(c(rep("ESC", 48), rep("MEF", 48)))
		
		# Quality control on individual cells.
		totals <- colSums(counts[!is.spike,])
		is.mito <- grepl("^mt-", rownames(counts)) & !is.spike
		okay.libs <- totals >= 1e5 & colSums(counts[is.mito,])/totals < 0.1 
		counts <- counts[,okay.libs]
		grouping <- grouping[okay.libs]
		
	} else if (datatype=="brennecke") {
		countsAll <- read.table('nmeth.2645-S10.csv.gz', header=TRUE, row.names=1, sep=',', colClasses=c('character', rep('integer', 13)))
		grouping <- factor(substr(colnames(countsAll), 0, 2))
		
		geneTypes <- factor( c( AT="At", pG="pGIBS", EN="HeLa", ER="ERCC" )[ substr( rownames(countsAll), 1, 2 ) ] )
		is.HeLa <- which( geneTypes=="HeLa" )
		countsHeLa <- countsAll[ is.HeLa, ]
		countsAt <- countsAll[ which( geneTypes=="At" ), ]
		
		rownames(countsHeLa) <- paste0("RNA_SPIKE_", rownames(countsHeLa)) # for SAMstrt, which fails unless we give it names.
		counts <- as.matrix(rbind(countsAt, countsHeLa))
		is.spike <- c(logical(nrow(countsAt)), !logical(nrow(countsHeLa)))
	}

	design <- model.matrix(~grouping)
	colnames(design) <- levels(grouping)
    spike.param <- spikeParam(counts[is.spike,])
	filter.keep <- rowSums(counts) >= ncol(counts)

	#################################################################################
	# Running edgeR first.

	# Filtering out crappy genes.
	y.ref <- DGEList(counts[!is.spike & filter.keep,])
	diag.done <- FALSE

    for (my.var in c(0.001, 0.01, 0.1)) {
		set.seed(231234)
		results <- top.res <- list()
		
		for (i in seq(0,20)) {
			# Running original data, then resampled versions.
			if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
			else { spike.data <- spike.param$counts }
		
			# Normalizing on the spike-in counts.
			y <- y.ref
			y$samples$norm.factors <- colSums(spike.data)/y$samples$lib.size

			# Running QL edgeR.
			y <- estimateDisp(y, design)
			fit <- glmQLFit(y, design, robust=TRUE)
			res <- glmQLFTest(fit)

			# Check if diagnostics look okay:
			if (!diag.done) {
				pdf(paste0("diagnostics_", datatype, ".pdf"))
				plotBCV(y) # Nice downward trend
				plotMDS(cpm(y, log=TRUE), col=c("red", "blue")[(as.integer(grouping)==1)+1L]) # Mostly separated
				plotQLDisp(fit) # Funny parabolic shape, but typical of single-cell data where there's loads of zeros.
				dev.off()
			}

			# Comparing.
			chosen <- p.adjust(res$table$PValue, method="BH") <= fdr.threshold
		    my.rank <- rank(res$table$PValue) 
            top.ranked <- lapply(top.hits, function(x) { my.rank <= x })
			if (i) { 
				lost <- sum(original & !chosen)/sum(original)
				gained <- sum(!original & chosen)/sum(chosen)
				results[[i]] <- c(lost, gained)
				top.res[[i]] <- sapply(seq_along(top.hits), function(j) { sum(best[[j]] & !top.ranked[[j]])/top.hits[j] })
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
        colnames(final) <- c("Lost", "Gained")
		top.lost <- do.call(rbind, top.res)
        colnames(top.lost) <- paste0("Top", top.hits)

        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="edgeR", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}

	#################################################################################
	# Also running NB models without any EB shrinkage, only using the tagwise dispersions.

	diag.done <- FALSE
    for (my.var in c(0.001, 0.01, 0.1)) {
		set.seed(231234)
		results <- top.res <- list()
		
		for (i in seq(0,20)) {
			# Running original data, then resampled versions.
			if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
			else { spike.data <- spike.param$counts }
		
			# Normalizing on the spike-in counts.
			y <- y.ref
			y$samples$norm.factors <- colSums(spike.data)/y$samples$lib.size

			# Running no-shrinkage edgeR.
			y <- estimateDisp(y, design, prior.df=0, trend.method="none")
			fit <- glmFit(y, design)
			res <- glmLRT(fit)

			# Comparing.
			chosen <- p.adjust(res$table$PValue, method="BH") <= fdr.threshold
            my.rank <- rank(res$table$PValue) 
            top.ranked <- lapply(top.hits, function(x) { my.rank <= x })
			if (i) { 
				lost <- sum(original & !chosen)/sum(original)
				gained <- sum(!original & chosen)/sum(chosen)
				results[[i]] <- c(lost, gained)
				top.res[[i]] <- sapply(seq_along(top.hits), function(j) { sum(best[[j]] & !top.ranked[[j]])/top.hits[j] })
			} else {
				original <- chosen
				best <- top.ranked
			}

			if (!diag.done) {
				# Top genes also look like the ones listed in Islam's paper, which is good, e.g. Sparc, S100a6, Vim, Fn1.
				de.out <- topTags(res, n=Inf)
				write.table(file=paste0("outNB_", datatype, ".tsv"), de.out, row.names=TRUE, 
					quote=FALSE, sep="\t", col.names=NA)
				diag.done <- TRUE
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- c("Lost", "Gained")
		top.lost <- do.call(rbind, top.res)
        colnames(top.lost) <- paste0("Top", top.hits)
   
        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="edgeR (LRT)", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}

	#################################################################################
	# Finally, monocle.

	for (my.var in c(0.001, 0.01, 0.1)) {
		set.seed(3742)
		results <- top.set <- list()
		
		for (i in seq(0,20)) {
			# Running original data, then resampled versions.
			if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
			else { spike.data <- spike.param$counts }

			# Normalizing by spike-in total counts.
			spike.totals <- colSums(spike.data)
			keep <- !is.spike & filter.keep
			normalized <- cpm(counts[keep,], lib.size=spike.totals)

			pdat <- AnnotatedDataFrame(data=data.frame(grouping=grouping))
			sampleNames(pdat) <- colnames(normalized)
			HSMM <- new("CellDataSet", exprs=normalized, phenoData=pdat)

            if (datatype=="islam") { # for convenience, otherwise we'll be here for days. 
                HSMM <- HSMM[1:2000,]
            }
			out <- differentialGeneTest(HSMM, fullModelFormulaStr="expression~grouping", cores=10) 

			# Choosing stuff. 
			chosen <- out$qval <= fdr.threshold
            my.rank <- rank(out$pval) 
            top.ranked <- lapply(top.hits, function(x) { my.rank <= x })
			if (i) { 
				lost <- sum(original & !chosen)/sum(original)
				gained <- sum(!original & chosen)/sum(chosen)
				results[[i]] <- c(lost, gained)
                top.res[[i]] <- sapply(seq_along(top.hits), function(j) { sum(best[[j]] & !top.ranked[[j]])/top.hits[j] })
            } else {
				original <- chosen
                best <- top.ranked
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- c("Lost", "Gained")
		top.lost <- do.call(rbind, top.res)
        colnames(top.lost) <- paste0("Top", top.hits)
   
        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="monocle", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}
}

#################################################################################

file.rename(temp, "results.txt")
sessionInfo()

###########################################################################

