# This looks at the properties of the HVGs or DE genes that
# were missed between the original and one simulation iteration.

# For HVGs with the log-variance method.

pdf("comparison_log_HVG.pdf", width=4, height=6)
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (x in c("calero", "liora")) {
    original <- sprintf("../variance/log_%s_0.tsv",x)
    modified <- sprintf("../variance/log_%s_1.tsv",x)
    odata <- read.table(original, header=TRUE)
    mdata <- read.table(modified, header=TRUE)
    
    kept <- rownames(odata) %in% rownames(mdata)
    lost <- ! rownames(odata) %in% rownames(mdata)
    dname <- ifelse(x=="calero", "416B", "TSC")

    boxplot(list(Detected=odata$mean[kept], Lost=odata$mean[lost]), ylab=expression("Log"[2]~"Mean count"),
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
    boxplot(list(Detected=odata$bio[kept], Lost=odata$bio[lost]), ylab="Biological component of variance",
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
}
dev.off()   

# For HVGs with the CV2 method.

pdf("comparison_CV2_HVG.pdf", width=4, height=6)
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (x in c("calero", "liora")) {
    original <- sprintf("../variance/cv2_%s_0.tsv",x)
    modified <- sprintf("../variance/cv2_%s_1.tsv",x)
    odata <- read.table(original, header=TRUE)
    mdata <- read.table(modified, header=TRUE)
    
    kept <- rownames(odata) %in% rownames(mdata)
    lost <- ! rownames(odata) %in% rownames(mdata)
    dname <- ifelse(x=="calero", "416B", "TSC")

    log.mean <- log2(odata$mean)
    boxplot(list(Detected=log.mean[kept], Lost=log.mean[lost]), ylab=expression("Log"[2]~"Mean count"),
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
    ratio <- odata$cv2/odata$trend
    boxplot(list(Detected=ratio[kept], Lost=ratio[lost]), ylab="Variance ratio (Total:Technical)", 
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, log='y', range=0)
}
dev.off()

# For DEGs with edgeR.

pdf("comparison_edgeR_DEG.pdf", width=4, height=6)
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (x in c("calero", "islam")) {
    original <- sprintf("../diffexp/edgeR_%s_0.tsv",x)
    modified <- sprintf("../diffexp/edgeR_%s_1.tsv",x)
    odata <- read.table(original, header=TRUE)
    odata <- odata[odata$FDR <= 0.05,]
    mdata <- read.table(modified, header=TRUE)
    mdata <- mdata[mdata$FDR <= 0.05,]
    
    kept <- rownames(odata) %in% rownames(mdata)
    lost <- ! rownames(odata) %in% rownames(mdata)
    dname <- ifelse(x=="calero", "416B", "mESC")

    boxplot(list(Detected=odata$logCPM[kept], Lost=odata$logCPM[lost]), ylab=expression("Average Log"[2]~"CPM"),
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
    lfc <- abs(odata$logFC)
    boxplot(list(Detected=lfc[kept], Lost=lfc[lost]), ylab=expression("Absolute Log"[2]~"Fold Change"),
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
}
dev.off()   

# For DEGs with MAST.

pdf("comparison_MAST_DEG.pdf", width=4, height=6)
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (x in c("calero", "islam")) {
    original <- sprintf("../diffexp/MAST_%s_0.tsv",x)
    modified <- sprintf("../diffexp/MAST_%s_1.tsv",x)
    odata <- read.table(original, header=TRUE)
    odata <- odata[odata$FDR <= 0.05,]
    mdata <- read.table(modified, header=TRUE)
    mdata <- mdata[mdata$FDR <= 0.05,]
    
    kept <- rownames(odata) %in% rownames(mdata)
    lost <- ! rownames(odata) %in% rownames(mdata)
    dname <- ifelse(x=="calero", "416B", "mESC")

    # Get statistics from edgeR, because MAST table isn't informative.
    edata <- read.table(sprintf("../diffexp/edgeR_%s_0.tsv",x))
    edata <- edata[rownames(odata),]

    boxplot(list(Detected=edata$logCPM[kept], Lost=edata$logCPM[lost]), ylab=expression("Average Log"[2]~"CPM"),
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
    lfc <- abs(edata$logFC)
    boxplot(list(Detected=lfc[kept], Lost=lfc[lost]), ylab=expression("Absolute Log"[2]~"Fold Change"),
            col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0)
}
dev.off()   


