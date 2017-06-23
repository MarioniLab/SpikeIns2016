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

for (datatype in c("calero", "islam")) {
	if (datatype=="islam") {
		incoming <- read.table("GSE29087_L139_expression_tab.txt.gz", 
			colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
		spike.in <- grepl("SPIKE", rownames(incoming))
		grouping <- rep(c("ESC", "MEF", "Neg"), c(48, 44, 4))
		
        is.mito <- grepl("^mt-", rownames(incoming)) & !spike.in
        
        expvar <- data.frame(Group=grouping, stringsAsFactors=FALSE)
		
	} else if (datatype=="calero") {
        incoming.1 <- readRDS("../../real/Calero/trial_20160113/analysis/full.rds")
        incoming.2 <- readRDS("../../real/Calero/trial_20160325/analysis/full.rds")
        stopifnot(identical(rownames(incoming.1), rownames(incoming.2)))
        incoming <- cbind(incoming.1$counts, incoming.2$counts)
        gdata <- incoming.1$genes

        keep <- !gdata$spike2 # Using only the ERCCs as the spike-ins.
        incoming <- incoming[keep,]
        gdata <- gdata[keep,]
        spike.in <- gdata$spike1 
        is.mito <- gdata$is.mito       

        # Constructing the experimental design.        
        expvar <- data.frame(Batch=rep(LETTERS[1:2], c(ncol(incoming.1), ncol(incoming.2))),
                             Group=ifelse(c(incoming.1$samples$induced, incoming.2$samples$induced), 
                                          "Induced", "Control"), 
                             stringsAsFactors=FALSE)
    }
    
    # Quality control on cells.
    totals <- colSums(incoming)
    okay.libs <- !isOutlier(totals, nmad=3, log=TRUE, type="lower") & 
        !isOutlier(colSums(incoming!=0), nmad=3, log=TRUE, type="lower") &
        !isOutlier(colSums(incoming[is.mito,])/totals, nmad=3, type="higher") & 
        !isOutlier(colSums(incoming[spike.in,])/totals, nmad=3, type="higher")
    incoming <- incoming[,okay.libs]
    expvar <- expvar[okay.libs,,drop=FALSE]

	filter.keep <- rowMeans(incoming) >= 1
    block <- do.call(paste, c(expvar, sep="."))
    spike.param <- spikeParam(incoming[spike.in & filter.keep,], design=model.matrix(~block))
    cell.counts <- incoming[!spike.in & filter.keep,]

	#################################################################################
	# Running edgeR first.

	y.ref <- DGEList(cell.counts)
    if (datatype=="calero") {
        design <- model.matrix(~Batch + Group, expvar)
    } else if (datatype=="islam") {
        design <- model.matrix(~Group, expvar)
    }

    for (my.var in c(0.015)) {
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
			if (i==0) {
				pdf(paste0("diagnostics_", datatype, ".pdf"))
				plotBCV(y) # Nice downward trend
				limma::plotMDS(edgeR::cpm(y, log=TRUE), pch=16,
                               col=c("red", "blue")[(as.integer(factor(expvar$Group))==1)+1L]) # Mostly separated
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

			if (i<=1) {
				# Top genes also look like the ones listed in Islam's paper, which is good, e.g. Sparc, S100a6, Vim, Fn1.
				de.out <- topTags(res, n=Inf)
				write.table(file=paste0("edgeR_", datatype, "_", i, ".tsv"), de.out, row.names=TRUE, 
					quote=FALSE, sep="\t", col.names=NA)
				diag.done <- TRUE
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- "AllDE"
		top.lost <- do.call(rbind, top.set)
        colnames(top.lost) <- paste0("Top", top.hits)

        write.table(file=temp, data.frame(Dataset=datatype, Variance=my.var, Method="edgeR", Total=length(original), Detected=sum(original),
            final, top.lost), row.names=FALSE, col.names=!filled, quote=FALSE, sep="\t", append=filled)
		filled <- TRUE
	}

    #################################################################################
    # Also, MAST.

    if (datatype=="calero") {
        form <- ~cngeneson + Batch + Group
    } else if (datatype=="islam") {
        form <- ~cngeneson + Group
    }

    for (my.var in c(0.015)) {
        set.seed(3742)
        results <- top.set <- list()
    
        for (i in seq(0,20)) {
            # Running original data, then resampled versions.
            if (i) { spike.data <- resampleSpikes(spike.param, var.log=my.var) }
            else { spike.data <- spike.param$counts }

            # Computing effective library sizes that reflect the relative spike-in coverage.
            spike.totals <- colSums(spike.data)
            eff.lib.size <- spike.totals/mean(spike.totals) * mean(colSums(cell.counts))
            
            # Computing log2-(CPM+1), as suggested in the vignette.
            lcpms <- log2(t(t(cell.counts)/eff.lib.size)*1e6 + 1)
            suppressMessages({
                sca <- FromMatrix(lcpms, cData=data.frame(wellKey=colnames(lcpms)), fData=data.frame(primerid=rownames(lcpms)))
                colData(sca) <- cbind(colData(sca), expvar)
                colData(sca)$cngeneson <- scale(colMeans(cell.counts>0))
                
                mfit <- zlm(form, sca)
                mlrt <- lrTest(mfit, "Group")
                MAST.p <- mlrt[, "hurdle", "Pr(>Chisq)"]
            })

            # Choosing significant genes.
            fdr <- p.adjust(MAST.p, method="BH")
			chosen <- fdr <= fdr.threshold
		    top.ranked <- rank(MAST.p) <= top.hits
			if (i) { 
				results[[i]] <- sum(original!=chosen)/sum(original)
				top.set[[i]] <- sum(best!=top.ranked)/top.hits
			} else {
				original <- chosen
				best <- top.ranked
			}

            # Storing diagnostics.
			if (i<=1) {
				write.table(file=paste0("MAST_", datatype, "_", i, ".tsv"), cbind(mlrt, FDR=fdr), 
                            row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
			}
		}

		final <- do.call(rbind, results)
        colnames(final) <- "AllDE"
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

