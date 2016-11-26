#################################################################################
# This script tests out dimensionality reduction and clustering.

suppressPackageStartupMessages(require(simpaler))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(dendextend))

# Estimate support for each cluster.

compute.max.jaccard <- function(ref.clust, incoming.clust) {
    collected <- list()
    for (kchar in names(ref.clust)) { # Iterating through all choices of 'k'
        cur.ref.clust <- ref.clust[[kchar]]
        cur.collected <- list()
        
        for (refc in seq_along(cur.ref.clust)) { # Iterating through all reference clusters.
            max.jaccard <- 0 
            my.ref.clust <- cur.ref.clust[[refc]]
            
            for (xi in incoming.clust[[kchar]]) { # Iterating through all observed clusters.
                max.jaccard <- max(max.jaccard, length(intersect(xi, my.ref.clust))/length(union(xi, my.ref.clust)))
            }
            cur.collected[[refc]] <- max.jaccard
        }
        collected[[kchar]] <- unlist(cur.collected)
    }
    return(collected)
}

#################################################################################
# Loading the pancreas data.

headers <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt.gz", nrow=1, comment="", stringsAsFactors=FALSE, sep="\t", fill=TRUE)
metadata <- read.table("E-MTAB-5061.sdrf.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
chosen <- metadata[,"Characteristics[individual]"]=="HP1502401"
cell.type <- metadata[chosen,"Characteristics[cell type]"]
source.name <- metadata[chosen,"Source Name"]

keep <- rep(list(NULL), nrow(metadata) * 2 + 2)
keep[length(headers) + which(headers %in% source.name)] <- list("integer")
keep[1] <- "character"
keep[2] <- "character"
my.cells <- cell.type[match(headers[headers %in% source.name], source.name)]

incoming <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt.gz", colClasses=keep, stringsAsFactors=FALSE, sep="\t", fill=TRUE)
discard <- duplicated(incoming[,2])
incoming <- incoming[!discard,]
rownames(incoming) <- incoming[,2]
incoming <- incoming[,-(1:2)]

spike.in <- grepl("ERCC", rownames(incoming))
totals <- colSums(incoming)

# Quality control on cells (couldn't find mitochondrial genes).
okay.libs <- !isOutlier(totals, nmad=3, log=TRUE, type="lower") & 
             !isOutlier(colSums(incoming!=0), nmad=3, log=TRUE, type="lower") &
             !isOutlier(colSums(incoming[spike.in,])/totals, nmad=3, type="higher")
incoming <- incoming[,okay.libs]
my.cells <- my.cells[okay.libs]

# Filtering out crappy spikes.
high.ab <- rowMeans(incoming) >= 1
countsCell <- as.matrix(incoming[high.ab & !spike.in,])
spike.param <- spikeParam(incoming[high.ab & spike.in,])

#################################################################################

set.seed(578354)
pca.collected <- list(list(), list())
clust.collected <- list()

for (i in seq(0,20)) {
    # Running original data, then resampled versions.
    if (i) { spike.data <- resampleSpikes(spike.param, var.log=0.01) }
    else { spike.data <- spike.param$counts }
  
    sce <- newSCESet(countData=rbind(countsCell, spike.data))
    sce <- calculateQCMetrics(sce, feature_controls=list(Spike=rep(c(FALSE, TRUE), c(nrow(countsCell), nrow(spike.data)))))
    isSpike(sce) <- "Spike"
    sce <- computeSpikeFactors(sce)
    sce <- normalize(sce)

    # Picking HVGs. 
    out <- trendVar(sce)
    out2 <- decomposeVar(sce, out)
    to.use <- out2$FDR <= 0.05 & out2$bio > 0.5 & !is.na(out2$p.value)

    #################################################################################
    # Running PCA with log-normalized expression values.

    cur.data <- exprs(sce)[to.use,]
    whee <- prcomp(t(cur.data), scale.=TRUE, center=TRUE) 
    
    current.co <- list(whee$x[,c(1,2)], whee$x[,c(1,3)])
    if (!i) { 
        ref <- current.co
    } else {
        # Computing the least-squares transformation (already column-mean-centered at zero)
        rotator <- function(par, coords) {
            angle <- par[1]
            scalex <- par[2]
            scaley <- par[3]
            transformation <- rbind(c(cos(angle)*scalex, -sin(angle)*scaley), c(sin(angle)*scalex, cos(angle)*scaley))
            coords %*% transformation
        }
        for (x in seq_along(current.co)) {
            cur.coords <- current.co[[x]]
            ref.coords <- ref[[x]]
            angled <- optim(par=c(angle=0, scalex=1, scaley=1), 
                            fn=function(par) { sum((rotator(par, cur.coords) - ref.coords)^2) })
            fitted <- rotator(angled$par, cur.coords)
            pca.collected[[x]][[i]] <- fitted
        }
    }

    #################################################################################
    # Clustering on the log-expression.

    euclid.tree <- hclust(dist(t(cur.data)), method = "ward.D2")
    all.out.clust <- list()
    for (k in c(2, 5, 10)) { 
        out.clust <- cutree(euclid.tree, k=k)
        all.out.clust[[as.character(k)]] <- split(names(out.clust), out.clust)
    }
   
    # Comparison of bootstrap support.
    if (!i) { 
        alt.tree <- hclust(dist(t(cur.data)))
        alt.clust <- list()
        for (k in c(2, 5, 10)) { 
            out.clust <- cutree(alt.tree, k=k)
            alt.clust[[as.character(k)]] <- split(names(out.clust), out.clust)
        }
        ref.clust <- all.out.clust
    } else {
        clust.collected[[i]] <- compute.max.jaccard(ref.clust, all.out.clust)
    }
}

#################################################################################
# Generating pretty plots.

all.cells <- factor(my.cells)
all.colors <- rainbow(nlevels(all.cells))
colors <- all.colors[as.integer(all.cells)]

# Visualization command.

visualize <- function(ref, incoming, col, mult=5, ...) {
	dists <- list()
	for (i in seq_along(incoming)) {
		dists[[i]] <- sqrt(rowSums((incoming[[i]]-ref)^2))
	}
	all.dists <- do.call(cbind, dists)
	r <- apply(all.dists, 1, FUN=quantile, probs=0.95)

	plot(0,0, xlim=range(ref[,1]), ylim=range(ref[,2]), type="n", ...)
	symbols(ref[,1], ref[,2], circles=r, add=TRUE, inches=FALSE, fg=col)
	points(ref[,1], ref[,2], col=col, pch=16, cex=0.5)
	invisible(NULL)
}

pdf(sprintf("pca_effect.pdf"), width=11, height=5)
layout(cbind(1,2,3), widths=c(2, 2, 1))
visualize(ref[[1]], pca.collected[[1]], col=colors, xlab="PC1", ylab="PC2", cex.axis=1.2, cex.lab=1.4)
visualize(ref[[2]], pca.collected[[2]], col=colors, xlab="PC1", ylab="PC3", cex.axis=1.2, cex.lab=1.4)
par(mar=c(5.1, 0.5, 4.1, 1.1))
plot(0,0,type="n", axes=FALSE, xlab="", ylab="")
legend("topleft", col=all.colors, pch=16, legend=levels(all.cells), cex=1.2)
dev.off()

#################################################################################
# Euclidean clustering (identify bootstrap intervals, basically).

all.means <- all.stder <- all.sizes <- list()
for (k in names(clust.collected[[1]])) {
    curk <- lapply(clust.collected, "[[", i=k)
    curk <- do.call(rbind, curk)
    all.means[[k]] <- colMeans(curk)    
    all.stder[[k]] <- sqrt(apply(curk, 2, var)/nrow(curk))
    all.sizes[[k]] <- unname(lengths(ref.clust[[k]]))
}

origin <- factor(rep(lengths(all.sizes), lengths(all.sizes)))
all.sizes <- unlist(all.sizes)
all.means <- unlist(all.means)
all.stder <- unlist(all.stder)
colors <- c("black", "red", "orange") 

pdf("clusters.pdf")
all.colors <- colors[as.integer(origin)]
plot(all.sizes, all.means, pch=16, col=all.colors, ylim=c(0, 1), ylab="Maximum Jaccard index", 
     xlab="Cluster size", cex.axis=1.2, cex.lab=1.4, cex=1.4, main="Effect of spike-in variability", cex.main=1.4)
upper <- all.means + all.stder 
segments(all.sizes, all.means, all.sizes, upper, col=all.colors)
segments(all.sizes-2, upper, all.sizes+2, upper, col=all.colors)
legend(max(all.sizes), 0, xjust=1, yjust=0, col=colors, legend=sprintf("k = %s", levels(origin)), pch=16, cex=1.4)

# Comparing to an alternative clustering method.
alt.jaccard <- compute.max.jaccard(ref.clust, alt.clust)
plot(all.sizes, unlist(alt.jaccard), col=all.colors, ylim=c(0, 1), ylab="Maximum Jaccard index", 
     xlab="Cluster size", cex.axis=1.2, cex.lab=1.4, cex=1.4, pch=16, main="Effect of clustering algorithm", cex.main=1.4)
dev.off()

#################################################################################

sessionInfo()

#################################################################################
