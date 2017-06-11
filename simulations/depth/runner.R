# This code estimates the increase in technical variance with decreasing depth.

library(simpaler)
getRescalingFunction <- function(d, design) {
    d <- estimateDisp(d, design, prior.df=0)
    fit <- glmFit(d, design)
    fitted <- exp(fit$unshrunk.coefficients %*% t(design) + mean(fit$offset))
    fitted[is.na(fitted)] <- 0

    means <- rowMeans(fitted)
    disp <- fit$dispersion
    lmean <- log(means)
    ldisp <- log(disp)
    keep <- is.finite(lmean)
    lmean <- lmean[keep]
    ldisp <- ldisp[keep]
    
    trend <- loess(ldisp ~ lmean, degree=1)
    function(scale) {
        new.fitted <- fitted*scale
        new.means <- rowMeans(new.fitted)
        lnew.means <- log(new.means)[keep]
        
        new.disp <- disp
        new.disp[keep] <- exp(predict(trend, data.frame(lmean = lnew.means)) + residuals(trend))
        upper.bound <- exp(max(fitted(trend))) # Setting to the upper bound at low counts.
        new.disp <- pmin(upper.bound, new.disp)
        new.disp[is.na(new.disp)] <- upper.bound
        return(list(fitted=new.fitted, disp=new.disp))
    }
}

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
    ercc.FUN <- getRescalingFunction(ercc.only, design)
    
    sirv.only <- DGEList(my.counts[grepl("SIRV", rownames(my.counts)),])
    sirv.FUN <- getRescalingFunction(sirv.only, design)
    
    # Simulating to estimate the effect of depth on the variability of ratios.
    # Note that we assume the offsets are true when estimating the dispersions.
    # This means dispersions may be underestimated for high-abundance genes (and overestimated for low-abundance genes). 
    
    set.seed(100)
    collected <- list()
    all.props <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
    for (p in seq_along(all.props)) { 
        prop <- all.props[p]
        to.collect <- list()
    
        for (it in seq_len(20)) { 
            new.ercc <- ercc.FUN(prop)
            sim.ercc <- matrix(rnbinom(length(new.ercc$fitted), mu=new.ercc$fitted, size=1/new.ercc$disp),
                               nrow=nrow(ercc.only), ncol=ncol(ercc.only))

            new.sirv <- sirv.FUN(prop)
            sim.sirv <- matrix(rnbinom(length(new.sirv$fitted), mu=new.sirv$fitted, size=1/new.sirv$disp),
                               nrow=nrow(sirv.only), ncol=ncol(sirv.only))
    
            ratios <- log2(colSums(sim.ercc)/colSums(sim.sirv))
            to.collect[[it]] <- estimateVariance(ratios=ratios, design=design)
        }
    
        to.collect <- unlist(to.collect)
        collected[[p]] <- to.collect
    }
    
    var.est <- sapply(collected, mean)
    var.sd <- sapply(collected, sd)/sqrt(lengths(collected))
    depth <- all.props*100
    upper <- var.est + var.sd
    
    plot(depth, var.est, xlab="Percentage of spike-in coverage",
         ylab=expression("Variance of"~log[2]~"[ERCC/SIRV]"),
         ylim=c(min(var.est), max(upper)), pch=16, log="y", 
         main=top, cex.axis=1.2, cex.lab=1.4, cex.main=1.5)
    lines(depth, var.est)
    segments(depth, var.est, y1=upper)
    segments(depth-1, upper, depth+1, upper)
}
  
dev.off()
