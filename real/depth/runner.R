# This code estimates the increase in technical variance with decreasing depth.

library(simpaler)
getRescalingFunction <- function(d, design) {
    # Estimating the NB dispersion, fitting a GLM.
    d <- estimateDisp(d, design, prior.df=0, trend.method="none")
    fit <- glmFit(d, design)

    # Calculating new fitted values with a common offset for all cells.
    fitted <- exp(fit$unshrunk.coefficients %*% t(design) + mean(fit$offset))
    fitted[is.na(fitted)] <- 0

    # Obtaining non-NA log-means and log-dispersions for trend fitting.
    means <- rowMeans(fitted)
    disp <- fit$dispersion
    lmean <- log(means)
    ldisp <- log(disp)
    keep <- is.finite(lmean)
    lmean <- lmean[keep]
    ldisp <- ldisp[keep]
    
    trend <- loess(ldisp ~ lmean, degree=1)

    # Returning a function to compute statistics given a downscaling factor.
    function(scale) {
        new.fitted <- fitted*scale
        new.means <- rowMeans(new.fitted)
        lnew.means <- log(new.means)[keep]
        
        new.disp <- rep(NA_real_, length(disp))
        new.disp[keep] <- exp(predict(trend, data.frame(lmean = lnew.means)) + residuals(trend))

        # Defining the upper bound as the largest observed variance
        # (or the left edge of the trend, whichever is largest).
        upper.bound <- exp(max(c(ldisp, fitted(trend)))) 
        new.disp <- pmin(upper.bound, new.disp)
        
        # predict.loess doesn't extrapolate, so we set NA values to the upper bound.
        new.disp[is.na(new.disp)] <- upper.bound
        return(list(fitted=new.fitted, disp=new.disp))
    }
}

# pdf("diagnostic.pdf")
# plot(lmean, ldisp, xlab="Log-mean", ylab="Log-Dispersion")
# o <- order(lmean)
# lines(lmean[o], fitted(trend)[o], col="red")
# out <- FUN(0.1)
# points(log(rowMeans(out$fitted)), log(out$disp), col="blue", pch=16, cex=0.5)
# dev.off()

pdf("depth_effect.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (operator in c("Calero", "Liora")) {
    if (operator=="Calero") {
        datasets <- c("trial_20160113", "trial_20160325")
    } else if (operator=="Liora") {
        datasets <- c("test_20160906", "test_20170201")
    }

    collected.cov <- collected.est <- collected.se <- list()
    all.props <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.02, 0.01)
    
    for (dataset in datasets) {
        loadin <- readRDS(file.path("../", operator, dataset, "analysis/object.rds"))
        keep <- loadin$samples$premixed
        my.counts <- loadin$counts[,keep]

        if (operator=="Liora") {
            design <- matrix(1, nrow=ncol(my.counts))
            top <- "TSC"
        } else {
            induction <- factor(loadin$samples$induced[keep])
            design <- model.matrix(~induction)
            top <- "416B"
        }

        library(edgeR)
        d <- DGEList(my.counts[grepl("ERCC|SIRV", rownames(my.counts)),])

        # Setting offsets for each transcript, based on the set of origin.
        is.ercc <- grepl("ERCC", rownames(d))
        is.sirv <- !is.ercc
        totals <- list(colSums(d$counts[is.ercc,,drop=FALSE]), colSums(d$counts[is.sirv,,drop=FALSE]))
        offsets <- do.call(rbind, totals[1 + is.sirv])
        d <- scaleOffset(d, log(offsets))

        FUN <- getRescalingFunction(d, design)
      
        # Simulating to estimate the effect of depth on the variability of ratios.
        # Note that we assume the offsets are true when estimating the dispersions.
        # This means dispersions may be underestimated for high-abundance genes (and overestimated for low-abundance genes). 
        
        set.seed(1000)
        collected <- list()
        for (p in seq_along(all.props)) { 
            prop <- all.props[p]
            to.collect <- list()
        
            for (it in seq_len(20)) { 
                new.spikes <- FUN(prop)

                ercc.fitted <- new.spikes$fitted[is.ercc,,drop=FALSE]
                ercc.disp <- new.spikes$disp[is.ercc]
                sim.ercc <- matrix(rnbinom(length(ercc.fitted), mu=ercc.fitted, size=1/ercc.disp),
                                   nrow=nrow(ercc.fitted), ncol=ncol(ercc.fitted))
    
                sirv.fitted <- new.spikes$fitted[is.sirv,,drop=FALSE]
                sirv.disp <- new.spikes$disp[is.sirv]
                sim.sirv <- matrix(rnbinom(length(sirv.fitted), mu=sirv.fitted, size=1/sirv.disp),
                                   nrow=nrow(sirv.fitted), ncol=ncol(sirv.fitted))
        
                ratios <- log2(colSums(sim.ercc)/colSums(sim.sirv))
                if (any(!is.finite(ratios))) {
                    warning("skipping iteration with non-finite ratios")
                    next
                }
                to.collect[[it]] <- estimateVariance(ratios=ratios, design=design)
            }
        
            to.collect <- unlist(to.collect)
            collected[[p]] <- to.collect
        }

        collected.cov[[dataset]] <- mean(colSums(d$counts))
        collected.est[[dataset]] <- sapply(collected, mean)
        collected.se[[dataset]] <- sapply(collected, sd)/sqrt(lengths(collected))
    }

    upper <- mapply("+", collected.est, collected.se, SIMPLIFY=FALSE)
    ylow <- min(sapply(collected.est, min))
    yhigh <- max(sapply(upper, max))
    xlow <- min(sapply(collected.cov, min))
    xhigh <- max(sapply(collected.cov, max))
    collected.cov <- lapply(collected.cov, FUN=function(x) { x * all.props/1000 })

    # Making a plot.   
    cols <- c("black", "grey70")
    plot(collected.cov[[1]], collected.est[[1]], xlab=expression("Total spike-in counts, ERCC + SIRV ("*10^3*")"),
         ylab=expression("Variance of"~log[2]~"[ERCC/SIRV]"),
         ylim=c(ylow, yhigh), pch=16, log="y", col=cols[1],
         cex.axis=1.2, cex.lab=1.4)
    points(collected.cov[[2]], collected.est[[2]], pch=17, col=cols[2])

    for (x in seq_along(collected.est)) {
        lines(collected.cov[[x]], collected.est[[x]], col=cols[x])
        segments(collected.cov[[x]], collected.est[[x]], y1=upper[[x]], col=cols[x])
        segments(collected.cov[[x]]-1, upper[[x]], collected.cov[[x]]+1, upper[[x]], col=cols[x])
    }
    legend("topright", col=cols, lwd=2, pch=16:17, sprintf("%s (%s)", top, c("I", "II")), cex=1.2)

    # Saving results somewhere.
    for (x in seq_along(datasets)) {
        write.table(file=paste0(operator, "_", datasets[x], ".txt"),
                    data.frame(Coverage=collected.cov[[x]], Variance=collected.est[[x]], SE=collected.se[[x]]),
                    sep="\t", quote=FALSE, row.names=FALSE)
    }
}
  
dev.off()
