#################################################################################
# This script tests out the different clustering methods.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(dendextend))

# Visualization command.

visualize <- function(ref, incoming, col, mult=5, ...) {
	dists <- list()
	for (i in seq_along(incoming)) {
		dists[[i]] <- sqrt(rowSums((incoming[[i]]-ref)^2))
	}
	all.dists <- do.call(cbind, dists)
	r <- apply(all.dists, 1, FUN=quantile, probs=0.95)
#	r <- apply(all.dists, 1, FUN=median) * mult # Effectively MAD; more robust to oddities.

	plot(0,0, xlim=range(ref[,1]), ylim=range(ref[,2]), type="n", ...)
	symbols(ref[,1], ref[,2], circles=r, add=TRUE, inches=FALSE, fg=col)
	points(ref[,1], ref[,2], col=col, pch=16, cex=0.5)
	invisible(NULL)
}

# Estimate support for each cluster.

support.tree <- function(ref, incoming) {
    dendrapply(ref, FUN=function(x) { 
        labs <- labels(x)
        out <- which_node(incoming, labs)
        if (is.null(attributes(x)$support)) { attributes(x)$support <- 0 }
        if (is.null(attributes(x)$support.95)) { attributes(x)$support.95 <- 0 }
        nl <- nleaves(x)
        nm <- get_nodes_attr(incoming, "members", out)
        if (nl==nm) { attributes(x)$support <- attributes(x)$support + 1 }
        if (nl >= nm * 0.95) { attributes(x)$support.95 <- attributes(x)$support.95 + 1 }
        return(x)
    })
}

summarize.support <- function(ref, support="support") {
    new.env <- environment()
    new.env$whee <- list()
    invisible(dendrapply(ref, FUN=function(x) {
        new.env$whee[[length(new.env$whee)+1]] <- c(nleaves(x), attributes(x)[[support]])
    }))
    out <- do.call(rbind, new.env$whee)
    keep <- out[,1]!=1L & out[,1]!=ncol(subcounts)
    return(list(size=out[keep,1], support=out[keep,2]))
}

#################################################################################

for (datatype in c("ola")) { # "scialdone")) { 
	if (datatype=="scialdone") { 
		counts <- read.table("data_2_4_15947.txt.gz", header=TRUE, row.names=1, colClasses=c("character", rep("integer", 42)))
		is.spike <- grepl("^ERCC", rownames(counts))
		is.mouse <- grepl("^ENSMUSG", rownames(counts))

		# Quality control on cells (ignored; all mitochondrial proportions are fairly high, here).
#		chr.loc <- findChr(rownames(counts[is.mouse,]), TxDb.Mmusculus.UCSC.mm10.ensGene)
#		lib.sizes <- colSums(counts[is.mouse,])
#		okay.libs <- lib.sizes > 1e5 & colSums(counts[is.mouse,][!is.na(chr.loc) & chr.loc=="chrM",])/lib.sizes < 0.1
#		counts <- counts[,okay.libs]

		# Only using higher-abundance genes, to avoid problems with technical noise.
		subcounts <- counts[is.mouse & rowSums(counts) > 100,]

		colors <- c(X2="grey", Early="red", EE="darkgreen", ME="blue")[sub("_.*", "", colnames(subcounts))]
		desired.dim <- c(1,2)

	} else if (datatype=="ola") {
        all.counts <- read.table("ESpresso/counttable_es.csv", header=TRUE, row.names=1, colClasses=c("character", rep("integer", 704)))
        serum <- sub("ola_mES_([^_]+)_.*", "\\1", colnames(all.counts))
        batch <- sub("ola_mES_[^_]+_([^_]+)_.*", "\\1", colnames(all.counts))
        targets <- data.frame(Serum=serum, Batch=batch)
        
        # Only using data from the batch with spike-ins.
        keep <- targets$Batch=="3"
		counts <- all.counts[,keep]
	
		is.spike <- grepl("^ERCC", rownames(counts))
		is.mouse <- grepl("^ENSMUSG", rownames(counts))
	
	    # Only using higher-abundance genes, to avoid problems with technical noise (QC on cells is already done).
		subcounts <- counts[is.mouse & rowSums(counts) > ncol(counts),]

		colors <- c("2i"="blue", a2i="red", lif="grey50")[as.character(targets$Serum[keep])]
		desired.dim <- c(1,2)
	}
	spike.param <- spikeParam(counts[is.spike,])

	#################################################################################
	# Performing PCA to reduce dimensions.

	set.seed(578354)
	collected <- list()

	for (i in seq(0,20)) {
		# Running original data, then resampled versions.
		if (i) { spike.data <- resampleSpikes(spike.param, var.log=0.01) }
		else { spike.data <- spike.param$counts }

		# Running PCA with log-CPMs (unadjusted, to keep things simple).
        cur.data <- cpm(subcounts, lib.size=colSums(spike.data), prior.count=1, log=TRUE)
		whee <- prcomp(t(cur.data), scale.=TRUE, center=TRUE) 

		current.co <- whee$x[,desired.dim] 
		if (!i) { 
			ref <- current.co
		} else {
			# Computing the least-squares transformation (already column-mean-centered at zero)
            rotator <- function(par) {
                angle <- par[1]
                scalex <- par[2]
                scaley <- par[3]
                transformation <- rbind(c(cos(angle)*scalex, -sin(angle)*scaley), c(sin(angle)*scalex, cos(angle)*scaley))
                current.co %*% transformation
            }
            angled <- optim(par=c(angle=0, scalex=1, scaley=1), fn=function(par) { sum((rotator(par) - ref)^2) })
			fitted <- rotator(angled$par)
			collected[[i]] <- fitted
		}
	}

	pdf(sprintf("pca_effect_%s.pdf", datatype))
	visualize(ref, collected, col=colors, xlab=paste0("PC", desired.dim[1]), ylab=paste0("PC", desired.dim[2]))
	dev.off()

	#################################################################################
    # Euclidean clustering (identify bootstrap intervals, basically).

    total.its <- 20
	for (i in seq(0,total.its)) {
		# Running original data, then resampled versions.
		if (i) { spike.data <- resampleSpikes(spike.param, var.log=0.01) }
		else { spike.data <- spike.param$counts }

        # Euclidean distance and correlation-based clustering.
        cur.data <- cpm(subcounts, lib.size=colSums(spike.data), prior.count=1, log=TRUE)
        euclid.tree <- hclust(dist(t(cur.data)))
        euclid.tree <- as.dendrogram(euclid.tree)      
        distM <- as.dist(1 - cor(cur.data, method = "spearman"))
        cor.tree <- as.dendrogram(hclust(distM))

        # Comparison of bootstrap support.
        if (!i) { 
			ref.euclid <- euclid.tree
			ref.cor <- cor.tree
		} else {
            ref.euclid <- support.tree(ref.euclid, euclid.tree)
            ref.cor <- support.tree(ref.cor, cor.tree)
		}
	}

    # Plotting cluster size against support.
    out.e <- summarize.support(ref.euclid, "support")
    out.c <- summarize.support(ref.cor, "support")
    pdf(sprintf("clusters_%s_pure.pdf", datatype))
    plot(out.e$size, out.e$support/total.its, pch=16, log="x", ylim=c(0, 1),
         xlab="Cluster size", ylab="Proportion of occurrences", cex.axis=1.2, cex.lab=1.4, main="Perfect match")
    points(out.c$size, out.c$support/total.its, pch=17, col="grey70")
    dev.off()

    out.e <- summarize.support(ref.euclid, "support.95")
    out.c <- summarize.support(ref.cor, "support.95")
    pdf(sprintf("clusters_%s_95.pdf", datatype))
    plot(out.e$size, out.e$support/total.its, pch=16, log="x", ylim=c(0, 1),
         xlab="Cluster size", ylab="Proportion of occurrences", cex.axis=1.2, cex.lab=1.4, main="95% match")
    points(out.c$size, out.c$support/total.its, pch=17, col="grey70")
    dev.off()
}

#################################################################################

sessionInfo()

#################################################################################
