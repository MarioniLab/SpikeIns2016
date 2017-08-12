# Making a plot.

results <- read.table("collated.txt", header=TRUE)
dat <- rbind(results$Original, results$Mean)
colnames(dat) <- results$Dataset

pdf("collated.pdf", width=10, height=6)
par(mar=c(7.1, 4.1, 2.1, 7.5), xpd=TRUE)
col <- c("grey80", "grey20")
x <- barplot(dat, beside=TRUE, las=2, ylab="Variance estimate", col=col, 
             cex.lab=1.4, cex.axis=1.2, cex.names=1.2)
legend(max(x), max(dat), legend=c("Size factors", "Sampling noise"), 
       fill=col, cex=1.2)
dev.off()
