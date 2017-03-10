###############################################################

dir.create("pics")
names <- c("416B (I)", "416B (II)", "TSC (I)", "TSC (II)")

###############################################################
# Making a barplot of the variance estimates.

require(simpaler)
total <- premixed <- ERCC.first <- SIRV.first <- volume <- ERCC.sf <- SIRV.sf <- list()
total.err <- premixed.err <- ERCC.first.err <- SIRV.first.err <- volume.err <- ERCC.sf.err <- SIRV.sf.err <- list()
total.df <- premixed.df <- ERCC.first.df <- SIRV.first.df <- volume.df <- ERCC.sf.df <- SIRV.sf.df <- list()
index <- 1L

for (operator in c("Calero", "Liora")) {
    if (operator=="Calero") {
        datasets <- c("trial_20160113", "trial_20160325")
    } else if (operator=="Liora") {
        datasets <- c("test_20160906", "test_20170201")
    }
    
    for (dataset in datasets) {
        out <- readRDS(file.path(operator, dataset, "analysis", "results.rds"))
        total[[index]] <- out$ratioSep.var
        total.err[[index]] <- attributes(out$ratioSep.var)$standard.error
        total.df[[index]] <- attributes(out$ratioSep.var)$df
        
        premixed[[index]] <- out$ratioPre.var
        premixed.err[[index]] <- attributes(out$ratioPre.var)$standard.error
        premixed.df[[index]] <- attributes(out$ratioPre.var)$df
        
        volume[[index]] <- out$ratioVol.var
        volume.err[[index]] <- attributes(out$ratioVol.var)$standard.error
        volume.df[[index]] <- attributes(out$ratioVol.var)$df
        
        ERCC.sf[[index]] <- out$sfERCC.var
        ERCC.sf.err[[index]] <- attributes(out$sfERCC.var)$standard.error
        ERCC.sf.df[[index]] <- attributes(out$sfERCC.var)$df
        
        SIRV.sf[[index]] <- out$sfSIRV.var
        SIRV.sf.err[[index]] <- attributes(out$sfSIRV.var)$standard.error
        SIRV.sf.df[[index]] <- attributes(out$sfSIRV.var)$df
        
        ERCC.first[[index]] <- out$ratioERCCfirst.var 
        ERCC.first.err[[index]] <- attributes(out$ratioERCCfirst.var)$standard.error
        ERCC.first.df[[index]] <- attributes(out$ratioERCCfirst.var)$df
        
        SIRV.first[[index]] <- out$ratioERCCsecond.var
        SIRV.first.err[[index]] <- attributes(out$ratioERCCsecond.var)$standard.error
        SIRV.first.df[[index]] <- attributes(out$ratioERCCsecond.var)$df
        index <- index + 1L
    }
}

pdf("pics/variance_exp.pdf", width=14, height=7)
par(mar=c(5.1, 5.1, 2.1, 2.1), mfrow=c(1,2), xpd=TRUE)

# First making a plot of separate vs premixed.

cols <- c("orange", "purple", "grey80")
final <- rbind(Separate=unlist(total), Premixed=unlist(premixed), Volume=unlist(volume))
colnames(final) <- names
final.err <- rbind(Separate=unlist(total.err), Premixed=unlist(premixed.err), Volume=unlist(volume.err))
upper.limit <- final  + final.err
final.df <- rbind(Separate=unlist(total.df), Premixed=unlist(premixed.df), Volume=unlist(volume.df))

ylim <- 0.03
out <- barplot(final, beside=TRUE, ylab=expression("Variance of"~log[2]~"[ERCC/SIRV]"), 
               cex.axis=1.2, cex.lab=1.4, cex.names=1.4, col=cols, ylim=c(0, ylim))
segments(out, final, y1=upper.limit)
segments(out-0.1, upper.limit, out+0.1)
text(out[-3,], upper.limit[-3,], final.df[-3,], pos=3, cex=0.9, offset=0.2)
legend(out[1]-0.5, ylim, fill=cols, rownames(final), cex=1.2)
curcoords <- par()$usr
mtext("a", line=0, cex=1.5, at=curcoords[1] - 0.14*(curcoords[2] - curcoords[1]), font=2)

final <- rbind(ERCC=unlist(ERCC.sf), SIRV=unlist(SIRV.sf))
colnames(final) <- names
final.err <- rbind(ERCC=unlist(ERCC.sf.err), SIRV=unlist(SIRV.sf.err))
upper.limit <- final  + final.err
final.df <- rbind(ERCC=unlist(ERCC.sf.df), SIRV=unlist(SIRV.sf.df))

# Making a plot of the size factors.

my.cols <- c("grey30", "grey70")
out <- barplot(final, beside=TRUE, ylab=expression("Variance of"~log[2]~"size factors"), 
        cex.axis=1.2, cex.lab=1.4, cex.names=1.4, col=my.cols, ylim=c(0, 0.3))
segments(out, final, y1=upper.limit)
segments(out-0.1, upper.limit, out+0.1)
text(colMeans(out), apply(upper.limit, 2, max), final.df[1,], pos=3, cex=0.9, offset=0.2)
legend("topright", fill=my.cols, rownames(final), cex=1.2)
curcoords <- par()$usr
mtext("b", line=0, cex=1.5, at=curcoords[1] - 0.14*(curcoords[2] - curcoords[1]), font=2)
dev.off()

# Ordering the elements.

pdf("pics/variance_order.pdf", width=9, height=7)
par(mar=c(5.1, 5.1, 2.1, 10.1), xpd=TRUE)

final <- rbind("ERCC+SIRV"=unlist(ERCC.first), "SIRV+ERCC"=unlist(SIRV.first))
colnames(final) <- names
final.err <- rbind("ERCC+SIRV"=unlist(ERCC.first.err), "SIRV+ERCC"=unlist(SIRV.first.err))
upper.limit <- final  + final.err
final.df <- rbind("ERCC+SIRV"=unlist(ERCC.first.df), "SIRV+ERCC"=unlist(SIRV.first.df))

out <- barplot(final, beside=TRUE, ylab=expression("Variance of"~log[2]~"[ERCC/SIRV]"), 
               cex.axis=1.2, cex.lab=1.4, cex.names=1.4, col=my.cols, ylim=c(0, max(upper.limit)))
segments(out, final, y1=upper.limit)
segments(out-0.1, upper.limit, out+0.1)
text(out, upper.limit, final.df, pos=3, cex=0.9, offset=0.2)
legend(max(out)+0.5, max(final), fill=my.cols, rownames(final), cex=1.2)
dev.off()

###############################################################
# Checking normality of the ratios.

separate <- premixed <- list()
index <- 1L
for (operator in c("Calero", "Liora")) {
    if (operator=="Calero") {
        datasets <- c("trial_20160113", "trial_20160325")
    } else if (operator=="Liora") {
        datasets <- c("test_20160906", "test_20170201")
    }

    for (dataset in datasets) {
        y <- readRDS(file.path(operator, dataset, "analysis", "object.rds"))

        sep <- y$samples$separate
        sep.groups <- factor(y$samples$group[sep])
        if (nlevels(sep.groups)==1L) { 
            design <- cbind(rep(1, length(sep.groups)))
        } else {
            design <- model.matrix(~sep.groups)
        }
        fit <- lm.fit(y$samples$ratio[sep], x=design)
        separate[[index]] <- sort(fit$residuals)

        pre <- y$samples$premixed
        pre.groups <- factor(y$samples$group[pre])
        if (nlevels(pre.groups)==1L) { 
            design <- cbind(rep(1, length(pre.groups)))
        } else {
            design <- model.matrix(~pre.groups)
        }
        fit <- lm.fit(y$samples$ratio[pre], x=design)
        premixed[[index]] <- sort(fit$residuals)
        
        index <- index + 1L
    }
}

pdf("pics/qq_separate.pdf")
collected <- list()
xlim <- ylim <- c(0, 0)
for (index in seq_along(separate)) {
    out <- qqnorm(separate[[index]], plot.it=FALSE)
    out$y <- out$y - mean(out$y)
    out$y <- out$y/sd(out$y)
    collected[[index]] <- out
    xlim <- pmax(xlim, abs(range(out$x)))
    ylim <- pmax(ylim, abs(range(out$y)))
}

qq.cols <- c("blue", "lightblue", "red", "pink")
plot(0,0,type="n", xlim=c(-xlim[1], xlim[2]), ylim=c(-ylim[1], ylim[2]), main="Separate addition", cex.main=1.4, 
     xlab="Theoretical quantiles", ylab="Sample quantiles", cex.axis=1.2, cex.lab=1.4)
abline(a=0, b=1, col="black", lwd=2, lty=2)
for (index in seq_along(collected)) {
    out <- collected[[index]]
    points(out$x, out$y, col=qq.cols[index], pch=16)
    lines(out$x, out$y, col=qq.cols[index], lwd=2)
}
legend("bottomright", col=qq.cols, lwd=2, pch=16, legend=names, cex=1.2)
dev.off()

pdf("pics/qq_premixed.pdf")
collected <- list()
xlim <- ylim <- c(0, 0)
for (index in seq_along(premixed)) {
    out <- qqnorm(premixed[[index]], plot.it=FALSE)
    out$y <- out$y - mean(out$y)
    out$y <- out$y/sd(out$y)
    collected[[index]] <- out
    xlim <- pmax(xlim, abs(range(out$x)))
    ylim <- pmax(ylim, abs(range(out$y)))
}

plot(0,0,type="n", xlim=c(-xlim[1], xlim[2]), ylim=c(-ylim[1], ylim[2]), main="Premixed addition", cex.main=1.4, 
     xlab="Theoretical quantiles", ylab="Sample quantiles", cex.axis=1.2, cex.lab=1.4)
abline(a=0, b=1, col="black", lwd=2, lty=2)
for (index in seq_along(collected)) {
    out <- collected[[index]]
    points(out$x, out$y, col=qq.cols[index], pch=16)
    lines(out$x, out$y, col=qq.cols[index], lwd=2)
}
dev.off()

###############################################################
# Checking that the sums for the spike-ins are the same between premixed and separate additions.

erccs <- sirvs <- list()
index <- 1L
for (operator in c("Calero", "Liora")) {
    if (operator=="Calero") {
        datasets <- c("trial_20160113", "trial_20160325")
    } else if (operator=="Liora") {
        datasets <- c("test_20160906", "test_20170201")
    }

    for (dataset in datasets) {
        y <- readRDS(file.path(operator, dataset, "analysis", "object.rds"))
        erccs[[index]] <- y$samples$sum1[y$samples$separate]
        sirvs[[index]] <- y$samples$sum2[y$samples$separate]
        index <- index + 1L
        erccs[[index]] <- y$samples$sum1[y$samples$premixed]
        sirvs[[index]] <- y$samples$sum2[y$samples$premixed]
        index <- index + 1L
    }
}

pdf("pics/total_ercc.pdf")
my.cols <- cols[1:2]
boxplot(erccs, at=rep(c(-0.5, 0.5), 4) + rep(1:4*3, each=2), xaxt="n", log="y", ylab="Total ERCC read count",
        col=rep(my.cols, 4), cex.axis=1.2, cex.lab=1.4)
axis(1, at=1:4*3, names, cex.axis=1.1)
legend("topright", fill=my.cols, legend=c("Separate", "Premixed"), cex=1.2)
dev.off()

pdf("pics/total_sirv.pdf")
boxplot(erccs, at=rep(c(-0.5, 0.5), 4) + rep(1:4*3, each=2), xaxt="n", log="y", ylab="Total SIRV read count",
        col=rep(my.cols, 4), cex.axis=1.2, cex.lab=1.4)
axis(1, at=1:4*3, names, cex.axis=1.1)
dev.off()

