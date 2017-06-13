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
    
    kept <- odata$bio[rownames(odata) %in% rownames(mdata)]
    diff <- c(odata$bio[! rownames(odata) %in% rownames(mdata)], mdata$bio[! rownames(mdata) %in% rownames(odata)])
    dname <- ifelse(x=="calero", "416B", "TSC")

    out <- boxplot(list(Common=kept, "Lost/Gained"=diff), ylab="Biological component of variance",
                   col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0,
                   ylim=c(0, max(max(kept), max(diff))+1))
    text(1:2, out$stats[5,], c(length(kept), length(diff)), pos=3)
}
dev.off()   

# For HVGs with the CV2 method.

#pdf("comparison_CV2_HVG.pdf", width=4, height=6)
#par(mar=c(5.1, 5.1, 4.1, 2.1))
#for (x in c("calero", "liora")) {
#    original <- sprintf("../variance/cv2_%s_0.tsv",x)
#    modified <- sprintf("../variance/cv2_%s_1.tsv",x)
#    odata <- read.table(original, header=TRUE)
#    mdata <- read.table(modified, header=TRUE)
#    
#    dname <- ifelse(x=="calero", "416B", "TSC")
#    oratio <- odata$cv2/odata$trend
#    mratio <- mdata$cv2/mdata$trend
#    kept <- oratio[rownames(odata) %in% rownames(mdata)]
#    diff <- c(oratio[! rownames(odata) %in% rownames(mdata)], mratio[! rownames(mdata) %in% rownames(odata)])
#
#    out <- boxplot(list(Common=kept, "Lost/Gained"=diff), ylab="Variance ratio (Total:Technical)", 
#                   col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, log='y', range=0,
#                   ylim=c(0, max(max(kept), max(diff))+1))
#    text(1:2, out$stats[5,], c(length(kept), length(diff)), pos=3)
#}
#dev.off()

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
    
    dname <- ifelse(x=="calero", "416B", "mESC")
    olfc <- abs(odata$logFC)
    mlfc <- abs(mdata$logFC)
    kept <- olfc[rownames(odata) %in% rownames(mdata)]
    diff <- c(olfc[! rownames(odata) %in% rownames(mdata)], mlfc[!rownames(mdata) %in% rownames(odata)])

    out <- boxplot(list(Common=kept, "Lost/Gained"=diff), ylab=expression("Absolute Log"[2]~"Fold Change"),
                   col="grey80", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, main=dname, range=0,
                   ylim=c(0, max(max(kept), max(diff))+1))
    text(1:2, out$stats[5,], c(length(kept), length(diff)), pos=3)
}
dev.off()   

