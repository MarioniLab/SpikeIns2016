# This checks the swapping in the Calero data set by examining the
# expression of the CBFB-MYH11 oncogene.

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

collected <- is.pos <- is.neg <- list()
for (dataset in c("test_20160906", "test_20170201")) {
    full.data <- readRDS(file.path("..", "Liora", dataset, "analysis", "full.rds"))
    ercc.sums <- colSums(full.data$counts[full.data$genes$spike1,])
#    sirv.sums <- log10(colSums(full.data$counts[full.data$genes$spike2,]))
    mouse.sums <- colSums(full.data$counts[full.data$genes$mouse,])
            
    dname <- ifelse(dataset=="test_20160906", "TSC (I)", "TSC (II)")
    collected[[dname]] <- log10(mouse.sums/ercc.sums)
    is.neg[[dname]] <- which(full.data$samples$control.well=="-")
    is.pos[[dname]] <- which(full.data$samples$control.well=="+")
#    points(1:3, c(ercc.sums[is.neg], sirv.sums[is.neg], mouse.sums[is.neg]), col="red", pch=16, cex=1.8)
#    points(1:3, c(ercc.sums[is.pos], sirv.sums[is.pos], mouse.sums[is.pos]), col="blue", pch=17, cex=1.8)
}

pdf("Liora_check.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
boxplot(collected, ylab=expression("Log"[10]~"(Mouse/ERCC)"),
        col="grey80", cex.axis=1.2, cex.names=1.2, cex.main=1.4, cex.lab=1.4)
for (x in seq_along(is.pos)) {
    points(x, collected[[x]][[is.neg[[x]]]], col="red", pch=16, cex=1.8)
    points(x, collected[[x]][[is.pos[[x]]]], col="blue", pch=17, cex=1.8)
}
dev.off()
