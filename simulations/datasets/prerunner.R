# This converts all files into a standard format for easy parsing.

library(simpaler)
is.first <- TRUE

for (dataset in c("Wilson", "Islam", "Scialdone", "Grun", "Buettner", "Kolodziejczyk", "Calero", "Liora")) {
    if (dataset=="Wilson") {
        incoming <- read.table('GSE61533_HTSEQ_count_results.tsv.gz', header=TRUE, row.names=1, colClasses=c('character', rep('integer', 96)))
        spike.in <- grepl('^ERCC', rownames(incoming))
        design <- cbind(rep(1, ncol(incoming)))

    } else if (dataset=="Islam") {
        incoming <- read.table("GSE29087_L139_expression_tab.txt.gz", 
                               colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
        spike.in <- grepl("SPIKE", rownames(incoming))
        grouping <- rep(c("ESC", "MEF", "Neg"), c(48, 44, 4))
        incoming <- incoming[,grouping!="Neg"]
        grouping <- grouping[grouping!="Neg"]
        design <- model.matrix(~grouping)

    } else if (dataset=="Scialdone") {
        incoming <- read.table("Marioni_lab_1_Jul_2015_gene_counts_table.txt", skip=1, row.names=1, colClasses=c("character", rep("integer", 96)))
        incoming <- incoming[!grepl("^__", rownames(incoming)),]
        spike.in <- grepl('^ERCC', rownames(incoming))
        design <- cbind(rep(1, ncol(incoming)))

    } else if (dataset=="Grun") {
        incoming <- read.table("GSE54695_data_transcript_counts.txt.gz", header=TRUE, row.names=1, sep="\t", 
                               colClasses=c("character", rep("numeric", 80), vector("list", 80), 
                                            rep("numeric", 80), vector("list", 80))) # keep only single cells.
        spike.in <- grepl("^ERCC", rownames(incoming))
        grouping <- sub("SC_([^_]+)_[0-9]+", "\\1", colnames(incoming))
        design <- model.matrix(~grouping)

    } else if (dataset=="Buettner") {
        incoming.G1 <- read.table("E-MTAB-2805/G1_singlecells_counts.txt", row.names=1, header=TRUE,
                                  colClasses=c("character", vector("list", 3), rep("integer", 96)))
        incoming.G2M <- read.table("E-MTAB-2805/G2M_singlecells_counts.txt", row.names=1, header=TRUE,
                                  colClasses=c("character", vector("list", 3), rep("integer", 96)))
        incoming.S <- read.table("E-MTAB-2805/S_singlecells_counts.txt", row.names=1, header=TRUE,
                                  colClasses=c("character", vector("list", 3), rep("integer", 96)))
        stopifnot(identical(rownames(incoming.G1), rownames(incoming.G2M)))
        stopifnot(identical(rownames(incoming.G1), rownames(incoming.S)))
        incoming <- cbind(incoming.G1, incoming.G2M, incoming.S)

        incoming <- incoming[seq_len(max(grep("^ERCC", rownames(incoming)))),]
        spike.in <- grepl("^ERCC", rownames(incoming))
        grouping <- factor(rep(c("G1", "G2M", "S"), c(ncol(incoming.G1), ncol(incoming.G2M), ncol(incoming.S))))
        design <- model.matrix(~grouping)

    } else if (grepl("Kolodziejczyk", dataset)) {
        incoming <- read.table("counttable_es.csv", row.names=1, header=TRUE, colClasses=c("character", rep("integer", 704)))
        incoming <- incoming[!grepl("^__", rownames(incoming)),]
        spike.in <- grepl('^ERCC', rownames(incoming))
        batch <- sub(".*_([0-9]+)_[0-9]+.counts", "\\1", colnames(incoming))
        condition <- sub("ola_mES_([^_]+)_.*", "\\1", colnames(incoming))

        keep <- batch=="3" # Keep only the batch with spike-ins for all three conditions.
        incoming <- incoming[,keep]
        condition <- condition[keep]
        design <- model.matrix(~condition)
        
    } else if (dataset=="Calero") {
        incoming.1 <- readRDS("../../real/Calero/trial_20160113/analysis/full.rds")
        incoming.2 <- readRDS("../../real/Calero/trial_20160325/analysis/full.rds")
        stopifnot(identical(rownames(incoming.1), rownames(incoming.2)))
        incoming <- cbind(incoming.1$counts, incoming.2$counts)
        gdata <- incoming.1$genes
        
        keep <- !gdata$spike2 # Using only the ERCCs as the spike-ins.
        incoming <- incoming[keep,]
        gdata <- gdata[keep,]
        spike.in <- gdata$spike1 

        Batch <- rep(LETTERS[1:2], c(ncol(incoming.1), ncol(incoming.2)))
        Group <- ifelse(c(incoming.1$samples$induced, incoming.2$samples$induced), "Induced", "Control")
        design <- model.matrix(~0 + Batch + Group)

    } else if (dataset=="Liora") {
        incoming.1 <- readRDS("../../real/Liora/test_20160906/analysis/full.rds")
        incoming.1 <- incoming.1[,is.na(incoming.1$samples$control.well)]
        incoming.2 <- readRDS("../../real/Liora/test_20170201/analysis/full.rds")
        incoming.2 <- incoming.2[,is.na(incoming.2$samples$control.well)]
        stopifnot(identical(rownames(incoming.1), rownames(incoming.2)))
        incoming <- cbind(incoming.1$counts, incoming.2$counts)
        gdata <- incoming.1$genes
        
        keep <- !gdata$spike2 # Using only the ERCCs as the spike-ins.
        incoming <- incoming[keep,]
        gdata <- gdata[keep,]
        spike.in <- gdata$spike1

        block <- rep(LETTERS[1:2], c(ncol(incoming.1), ncol(incoming.2)))
        design <- model.matrix(~block)

    } else {
        stop("Unknown data set")
    }

    saveRDS(file=paste0(dataset, ".rds"), list(counts=as.matrix(incoming), spikes=spike.in, design=design))
}



