#################################################################################
# This performs a differential expression analysis of various data sets, 
# exploiting the presence of spike-ins for normalization.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(library(MAST))

# Note that Kharchenko's SCDE package doesn't seem to support spike-in normalization.

fdr.threshold <- 0.05
top.hits <- 200

temp <- "temp.txt"
if (file.exists(temp))  { unlink(temp) }
filled <- FALSE 

#################################################################################

for (datatype in c("Calero", "Islam", "Buettner", "Grun", "Kolodziejczyk")) {
    val <- readRDS(file.path("../datasets", paste0(datatype, ".rds")))
    incoming <- val$counts
    spike.in <- val$spikes
    design <- val$design

	if (datatype=="Islam") {
        test.coef <- "groupingMEF"
	} else if (datatype=="Calero") {
        test.coef <- "GroupInduced"
    } else if (datatype=="Buettner") {
        test.coef <- "groupingG2M"
    } else if (datatype=="Grun") {
        test.coef <- "groupingserum"
    } else if (datatype=="Kolodziejczyk") {
        test.coef <- "conditiona2i"
    }   

    # Quality control on cells.
    totals <- colSums(incoming)
    f <- designAsFactor(design)
    okay.libs <- !isOutlier(totals, nmad=3, log=TRUE, type="lower", batch=f) & 
        !isOutlier(colSums(incoming!=0), nmad=3, log=TRUE, type="lower", batch=f) &
        !isOutlier(colSums(incoming[spike.in,])/totals, nmad=3, type="higher", batch=f)
    incoming <- incoming[,okay.libs]
    design <- design[okay.libs,,drop=FALSE]

	filter.keep <- calcAverage(incoming) >= 0.1
    spike.param <- spikeParam(incoming[spike.in & filter.keep,], design=design)
    cell.counts <- incoming[!spike.in & filter.keep,]

	#################################################################################
	# Running edgeR first.

    y.ref <- DGEList(cell.counts)
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
			y <- estimateDisp(y, design, robust=TRUE)
			fit <- glmFit(y, design)
			res <- glmLRT(fit, coef=test.coef)

			# Check if diagnostics look okay:
			if (i==0) {
				pdf(paste0("diagnostics_", datatype, ".pdf"))
				plotBCV(y) # Nice downward trend
				limma::plotMDS(edgeR::cpm(y, log=TRUE), pch=16,
                               col=c("red", "blue")[(as.integer(design[,test.coef])==1)+1L]) # Mostly separated
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

    redesign <- design
    colnames(redesign) <- gsub("[()]", "", colnames(design)) # Killing the braces around the intercept.
    new.form <- paste(c("~0 + cngeneson", colnames(redesign)), collapse="+")
    form <- eval(parse(text=new.form))

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
                colData(sca) <- cbind(colData(sca), DataFrame(redesign))
                colData(sca)$cngeneson <- scale(colMeans(cell.counts>0))
                
                mfit <- zlm(form, sca)
                mlrt <- lrTest(mfit, test.coef)
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

