pdf("setplot.pdf", width=14, height=14)
par(mar=c(3.1, 5.1, 2.1, 2.1), mfrow=c(2,1))

###############################################################
# Making plots of how DE results change.

de.results <- read.table("~/AaronDocs/Research/SpikeIns/simulations/diffexp/results.txt", header=TRUE, sep="\t", stringsAsFactor=FALSE)    
de.results <- de.results[de.results$Variance==0.015,]

l.means <- r.means <- l.stder <- r.stder <-  list()
for (dataset in unique(de.results$Dataset)) { 
    current <- de.results[de.results$Dataset==dataset,]
    c.lost <- split(current$AllDE*100, current$Method)
    c.rank <- split(current$Top200*100, current$Method)
    l.means[[dataset]] <- sapply(c.lost, FUN=mean)
    l.stder[[dataset]] <- sapply(c.lost, FUN=function(x) sqrt(var(x)/length(x)))
    r.means[[dataset]] <- sapply(c.rank, FUN=mean)
    r.stder[[dataset]] <- sapply(c.rank, FUN=function(x) sqrt(var(x)/length(x)))
}

l.means <- do.call(rbind, l.means)
r.means <- do.call(rbind, r.means)
l.stder <- do.call(rbind, l.stder)
r.stder <- do.call(rbind, r.stder)

combined.means <- rbind(l.means[,1], r.means[,1], l.means[,2], r.means[,2])
combined.stder <- rbind(l.stder[,1], r.stder[,1], l.stder[,2], r.stder[,2])
upper <- combined.means + combined.stder

renamed <- c("416B\n(induced/control)", "mESC\n(G2M/G1)", "mESC\n(serum/2i)", "mESC\n(a2i/2i)")[match(colnames(upper), c("Calero", "Buettner", "Grun", "Kolodziejczyk"))]
all.cols <- c("forestgreen", "lightgreen", "darkgoldenrod", "yellow2")
out <- barplot(combined.means, beside=TRUE, col=all.cols, ylab='Change in the DEG set (%)', 
               ylim=c(0, 10), cex.lab=1.4, cex.axis=1.2, names.arg=character(ncol(upper)))
mtext(at=colMeans(out), text=renamed, side=1, line=2, cex=1.4)

segments(out, combined.means, out, upper)
segments(out-0.2, upper, out+0.2, upper)

legend("topright", fill=all.cols, cex=1.4,
       legend=c(expression("edgeR, all"), 
                expression("edgeR, top 200"), 
                "MAST, all",
                "MAST, top 200"))

curcoords <- par()$usr
mtext("a", line=0, cex=1.5, at=curcoords[1]-1, font=2)

###############################################################
# Making plots of how HVG results change.

hvg.results <- read.table("~/AaronDocs/Research/SpikeIns/simulations/variance/results.txt", header=TRUE, sep="\t", stringsAsFactor=FALSE)    
hvg.results <- hvg.results[hvg.results$Variance==0.015,]

l.means <- r.means <- l.stder <- r.stder <-  list()
for (dataset in unique(hvg.results$Dataset)) { 
    current <- hvg.results[hvg.results$Dataset==dataset,]
    c.lost <- split(current$Sig*100, current$Method)
    c.rank <- split(current$Top200*100, current$Method)
    l.means[[dataset]] <- sapply(c.lost, FUN=mean)
    l.stder[[dataset]] <- sapply(c.lost, FUN=function(x) sqrt(var(x)/length(x)))
    r.means[[dataset]] <- sapply(c.rank, FUN=mean)
    r.stder[[dataset]] <- sapply(c.rank, FUN=function(x) sqrt(var(x)/length(x)))
}

l.means <- do.call(rbind, l.means)
r.means <- do.call(rbind, r.means)
l.stder <- do.call(rbind, l.stder)
r.stder <- do.call(rbind, r.stder)

combined.means <- rbind(l.means[,1], r.means[,1], l.means[,2], r.means[,2])
combined.means[combined.means < 0.05] <- 0.05
combined.stder <- rbind(l.stder[,1], r.stder[,1], l.stder[,2], r.stder[,2])
upper <- combined.means + combined.stder

colnames(combined.means) <- sub(" \\(Kolod\\)", "", colnames(combined.means))
all.cols <- c("blue", "lightblue", "red", "pink")
out <- barplot(combined.means, beside=TRUE, col=all.cols, ylab='Change in the HVG set (%)', 
               ylim=c(0, 10), cex.lab=1.4, cex.axis=1.2, cex.names=1.4)

segments(out, combined.means, out, upper)
segments(out-0.2, upper, out+0.2, upper)

legend("topright", fill=all.cols, cex=1.4,
       legend=c(expression("Brennecke et al. (CV"^2*"), all"), 
                expression("Brennecke et al. (CV"^2*"), top 200"), 
                "Variance of log-expression, all",
                "Variance of log-expression, top 200"))

curcoords <- par()$usr
mtext("b", line=0, cex=1.5, at=curcoords[1]-1, font=2)

###############################################################
# End.

dev.off()
