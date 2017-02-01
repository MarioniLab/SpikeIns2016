###############################################################
# Making plots of how DE results change.

de.results <- read.table("~/AaronDocs/Research/SpikeIns/simulations/diffexp/results.txt", header=TRUE, sep="\t", stringsAsFactor=FALSE)    
de.results <- de.results[de.results$Variance==0.01,]

calero.de <- de.results[de.results$Dataset=="calero",]
c.lost <- split(calero.de$AllDE*100, calero.de$Method)
c.rank <- split(calero.de$Top200*100, calero.de$Method)
cl.means <- sapply(c.lost, FUN=mean)
cl.stder <- sapply(c.lost, FUN=function(x) sqrt(var(x)/length(x)))
cr.means <- sapply(c.rank, FUN=mean)
cr.stder <- sapply(c.rank, FUN=function(x) sqrt(var(x)/length(x)))

islam.de <- de.results[de.results$Dataset=="islam",]
i.lost <- split(islam.de$AllDE*100, islam.de$Method)
i.rank <- split(islam.de$Top200*100, islam.de$Method)
il.means <- sapply(i.lost, FUN=mean)
il.stder <- sapply(i.lost, FUN=function(x) sqrt(var(x)/length(x)))
ir.means <- sapply(i.rank, FUN=mean)
ir.stder <- sapply(i.rank, FUN=function(x) sqrt(var(x)/length(x)))

pdf("setplot.pdf", width=14, height=7)
par(mar=c(3.1, 5.1, 2.1, 2.1), mfrow=c(1,2))
means <- rbind(cl.means, cr.means, il.means, ir.means)
spacing <- matrix(rep(c(0.2, 0), length(means)/2), ncol=2)
spacing[1,1] <- 0
spacing[1,2] <- 2

all.cols <- c("blue", "lightblue", "forestgreen", "lightgreen")
out <- barplot(means, beside=TRUE, col=all.cols, ylab='Change in the DEG set (%)', 
               ylim=c(0, 8), cex.lab=1.4, cex.axis=1.2, cex.names=1.4, space=spacing)

stderrs <- rbind(cl.stder, cr.stder, il.stder, ir.stder)
segments(out, means, out, means+stderrs)
segments(out-0.2, means+stderrs, out+0.2, means+stderrs)

legend(0, 8, fill=all.cols, cex=1.4, 
        legend=c("416B (all DEGs)", "416B (top 200)", "mESC/MEF (all DE)", "mESC/MEF (top 200)")) 

curcoords <- par()$usr
mtext("a", line=0, cex=1.5, at=curcoords[1] - 0.14*(curcoords[2] - curcoords[1]), font=2)

###############################################################
# Making plots of how HVG results change.

hvg.results <- read.table("~/AaronDocs/Research/SpikeIns/simulations/variance/results.txt", header=TRUE, sep="\t", stringsAsFactor=FALSE)    
hvg.results <- hvg.results[hvg.results$Variance==0.01,]

calero.hvg <- hvg.results[hvg.results$Dataset=="calero",]
c.lost <- split(calero.hvg$Sig*100, calero.hvg$Method)
c.rank <- split(calero.hvg$Top200*100, calero.hvg$Method)
cl.means <- sapply(c.lost, FUN=mean)
cl.stder <- sapply(c.lost, FUN=function(x) sqrt(var(x)/length(x)))
cr.means <- sapply(c.rank, FUN=mean)
cr.stder <- sapply(c.rank, FUN=function(x) sqrt(var(x)/length(x)))

liora.hvg <- hvg.results[hvg.results$Dataset=="liora",]
l.lost <- split(liora.hvg$Sig*100, liora.hvg$Method)
l.rank <- split(liora.hvg$Top200*100, liora.hvg$Method)
ll.means <- sapply(l.lost, FUN=mean)
ll.stder <- sapply(l.lost, FUN=function(x) sqrt(var(x)/length(x)))
lr.means <- sapply(l.rank, FUN=mean)
lr.stder <- sapply(l.rank, FUN=function(x) sqrt(var(x)/length(x)))

wilson.hvg <- hvg.results[hvg.results$Dataset=="wilson",]
w.lost <- split(wilson.hvg$Sig*100, wilson.hvg$Method)
w.rank <- split(wilson.hvg$Top200*100, wilson.hvg$Method)
wl.means <- sapply(w.lost, FUN=mean)
wl.stder <- sapply(w.lost, FUN=function(x) sqrt(var(x)/length(x)))
wr.means <- sapply(w.rank, FUN=mean)
wr.stder <- sapply(w.rank, FUN=function(x) sqrt(var(x)/length(x)))

means <- rbind(cl.means, cr.means, ll.means, lr.means, wl.means, wr.means)
spacing <- matrix(rep(c(0.2, 0), length(means)/2), ncol=2)
spacing[1,1] <- 0
spacing[1,2] <- 2

all.cols <- c("blue", "lightblue", "red", "pink", "darkgoldenrod1", "darkgoldenrod4")
out <- barplot(means, beside=TRUE, col=all.cols, space=spacing,
    names.arg=c(expression("Brennecke et al. (CV"^2*")"), "Variance of log-expression"),
    ylab='Change in the HVG set (%)', ylim=c(0, 15), cex.lab=1.4, cex.axis=1.2, cex.names=1.4)

stderrs <- rbind(cl.stder, cr.stder, ll.stder, lr.stder, wl.stder, wr.stder)
segments(out, means, out, means+stderrs)
segments(out-0.2, means+stderrs, out+0.2, means+stderrs)

legend("topright", fill=all.cols, legend=c("416B (all HVGs)", "416B (top 200)", 
                                          "Trophoblast (all HVGs)", "Trophoblast (top 200)", 
                                          "HSC (all HVGs)", "HSC (top 200"), cex=1.4)

curcoords <- par()$usr
mtext("b", line=0, cex=1.5, at=curcoords[1] - 0.14*(curcoords[2] - curcoords[1]), font=2)
dev.off()

###############################################################
# End.

