# This code estimates the increase in technical variance with decreasing depth.

library(simpaler)
pdf("depth_effect.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (person in c("Liora", "Calero")) { 

    if (person=="Liora") {
        loadin1 <- readRDS("../../real/Liora/test_20160906/analysis/object.rds")
        loadin2 <- readRDS("../../real/Liora/test_20170201/analysis/object.rds")
        keep1 <- loadin1$samples$premixed
        keep2 <- loadin2$samples$premixed
        my.counts <- cbind(loadin1$counts[,keep1], loadin2$counts[,keep2])
        batch <- factor(rep(1:2, c(sum(keep1), sum(keep2))))
        design <- model.matrix(~batch)
        top <- "TSC"

    } else if (person=="Calero") {
        loadin1 <- readRDS("../../real/Calero/trial_20160113/analysis/object.rds")
        loadin2 <- readRDS("../../real/Calero/trial_20160325/analysis/object.rds")
        keep1 <- loadin1$samples$premixed
        keep2 <- loadin2$samples$premixed
        my.counts <- cbind(loadin1$counts[,keep1], loadin2$counts[,keep2])
        batch <- factor(rep(1:2, c(sum(keep1), sum(keep2))))
        induction <- factor(c(loadin1$samples$induced[keep1], loadin2$samples$induced[keep2]))
        design <- model.matrix(~batch+induction)
        top <- "416B"

    }

    library(edgeR)
    ercc.only <- DGEList(my.counts[grepl("ERCC", rownames(my.counts)),])
    ercc.only <- estimateDisp(ercc.only, design, prior.df=0)
    ercc.fit <- glmFit(ercc.only, design)
    ercc.means <- exp(ercc.fit$unshrunk.coefficients %*% t(design) + mean(ercc.fit$offset))
    ercc.means[is.na(ercc.means)] <- 0
    
    sirv.only <- DGEList(my.counts[grepl("SIRV", rownames(my.counts)),])
    sirv.only <- estimateDisp(sirv.only, design, prior.df=0)
    sirv.fit <- glmFit(sirv.only, design)
    sirv.means <- exp(sirv.fit$unshrunk.coefficients %*% t(design) + mean(sirv.fit$offset))
    sirv.means[is.na(sirv.means)] <- 0
    
    # Simulating to estimate the effect of depth on the variability of ratios.
    # Note that we assume the offsets are true when estimating the dispersions.
    # This means dispersions may be underestimated for high-abundance genes (and overestimated for low-abundance genes). 
    
    set.seed(100)
    collected <- list()
    all.props <- seq_len(10)/10
    for (p in seq_along(all.props)) { 
        prop <- all.props[p]
        to.collect <- list()
    
        for (it in seq_len(20)) { 
            sim.ercc <- matrix(rnbinom(length(ercc.fit$fitted.values), mu=ercc.means*prop, size=1/ercc.fit$dispersion),
                               nrow=nrow(ercc.only), ncol=ncol(ercc.only))
            sim.sirv <- matrix(rnbinom(length(sirv.fit$fitted.values), mu=sirv.means*prop, size=1/sirv.fit$dispersion),
                               nrow=nrow(sirv.only), ncol=ncol(sirv.only))
    
            ratios <- log2(colSums(sim.ercc)/colSums(sim.sirv))
            to.collect[[it]] <- estimateVariance(ratios=ratios, design=design)
        }
    
        to.collect <- unlist(to.collect)
        collected[[p]] <- to.collect
    }
    
    var.est <- sapply(collected, mean)
    var.sd <- sapply(collected, sd)/sqrt(lengths(collected))
    depth <- all.props*mean(ercc.fit$samples$lib.size)
    upper <- var.est + var.sd
    
    plot(depth, var.est, xlab="Coverage of ERCC transcripts (reads/well)", 
         ylab=expression("Variance of"~log[2]~"[ERCC/SIRV]"),
         ylim=c(min(var.est), max(upper)), pch=16, 
         main=top, cex.axis=1.2, cex.lab=1.4, cex.main=1.5)
    lines(depth, var.est)
    segments(depth, var.est, y1=upper)
    segments(depth-1000, upper, depth+1000, upper)
}
  
dev.off()
