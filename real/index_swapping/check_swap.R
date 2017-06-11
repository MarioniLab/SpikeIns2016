# This checks the swapping in the Calero data set by examining the expression of the CBFB-MYH11 oncogene.
# If swapping was occurring, the log-fold change should be closer to zero in the 4000 data set.

pdf("Calero_check.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (dataset in c("trial_20160113", "trial_20160325")) {
    type <- ifelse(dataset=="trial_20160113", "HiSeq 2500", "HiSeq 4000")
    spike.data <- readRDS(file.path("..", "Calero", dataset, "analysis", "object.rds"))
    oncogene.cpm <- edgeR::cpm(spike.data, log=TRUE, prior.count=3)["CBFB-MYH11-mcherry",]
    by.cherry <- split(oncogene.cpm, ifelse(spike.data$samples$induced, "Induced", "Control"))
    boxplot(by.cherry, ylab=expression("Log"[2]~"CPM (CBFB-MYH11-mCherry)"), main=type,
            col="grey80", cex.axis=1.2, cex.names=1.2, cex.main=1.4, cex.lab=1.4)
    legend("bottomright", legend=sprintf("Log-fold change: %.3f", mean(by.cherry[[2]])-mean(by.cherry[[1]])),
           bty="n", cex=1.2)
}
dev.off()

# This checks the swapping in the Liora data set by examining the negative control wells.
# Swapping should be at least partly responsible for any mouse counts in the negative control.

pdf("Liora_check.pdf")
collected <- is.neg <- list()
for (dataset in c("test_20160906", "test_20170201")) {
    full.data <- readRDS(file.path("..", "Liora", dataset, "analysis", "full.rds"))
    mouse.sums <- colSums(full.data$counts[full.data$genes$mouse,])
    log.counts <- log10(mouse.sums)    

    out <- hist(log.counts, breaks=20, xlab=expression("Log"[10]~"total mouse counts"),
         cex.axis=1.2, cex.lab=1.4, cex.main=1.4, col="grey80",
         main=ifelse(dataset=="test_20160906", "TSC (I)", "TSC (II)"))
    is.neg <- which(full.data$samples$control.well=="-")
    i <- findInterval(log.counts[is.neg], out$breaks)
    points(out$mids[i], out$counts[i]+1, col="red", pch=25, cex=1.8, bg="red")
}
dev.off()

# Consider the proportion of mouse counts in the negative control to that of other cells.
# This represents the proportion of counts in each well due to switching (and other factors, e.g., cell-free RNA).
# Note that there's no need to consider cell-specific biases here, as swapping occurs after cell-specific capture processes.

