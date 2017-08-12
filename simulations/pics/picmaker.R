pdf("setplot.pdf", width=14, height=7)
par(mar=c(3.1, 5.1, 2.1, 2.1), mfrow=c(2,1))

###############################################################
# Making plots of how DE results change.

de.results <- read.table("~/AaronDocs/Research/SpikeIns/simulations/diffexp/results.txt", header=TRUE, sep="\t", stringsAsFactor=FALSE)    
de.results <- de.results[de.results$Variance==0.015,]

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

means <- rbind(cl.means, cr.means, il.means, ir.means)
spacing <- matrix(rep(c(0.2, 0), length(means)/2), ncol=2)
spacing[1,1] <- 0
spacing[1,2] <- 2

ylim <- 10
all.cols <- c("blue", "lightblue", "forestgreen", "lightgreen")
out <- barplot(means, beside=TRUE, col=all.cols, ylab='Change in the DEG set (%)', 
               ylim=c(0, ylim), cex.lab=1.4, cex.axis=1.2, cex.names=1.4, space=spacing)

stderrs <- rbind(cl.stder, cr.stder, il.stder, ir.stder)
segments(out, means, out, means+stderrs)
segments(out-0.2, means+stderrs, out+0.2, means+stderrs)

legend(0, ylim, fill=all.cols, cex=1.4, 
        legend=c("416B (all DEGs)", "416B (top 200)", "mESC (all DE)", "mESC (top 200)")) 

curcoords <- par()$usr
mtext("a", line=0, cex=1.5, at=curcoords[1] - 0.14*(curcoords[2] - curcoords[1]), font=2)

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
combined.stder <- rbind(l.stder[,1], r.stder[,1], l.stder[,2], r.stder[,2])
upper <- combined.means + combined.stder

all.cols <- c("blue", "lightblue", "red", "pink")
out <- barplot(combined.means, beside=TRUE, col=all.cols, ylab='Change in the HVG set (%)', 
               ylim=c(0, max(upper)), cex.lab=1.4, cex.axis=1.2, cex.names=1.4)

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
