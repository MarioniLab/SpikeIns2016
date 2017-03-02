########################################################################################
# Print out important bits from the metadata for Liora's lot.

collected <- list()
fpath <- "../Liora"
relink <- "make_links_Liora.sh"
write(file=relink, c("set -e", "set -u", "mkdir fastq_Liora", "cd fastq_Liora"), ncol=1)
library(edgeR)

for (sample in c("test_20160906", "test_20170201")) {
    cpath <- file.path(fpath, sample, "analysis", "genic_counts.tsv")
    all.files <- read.table(cpath, nrows=1, stringsAsFactor=FALSE, comment="")
    prefixes <- as.character(all.files[-c(1:2)])
    batch <- sub(".*_", "", basename(sample))

    # Loading in the corresponding object.
    full.obj <- readRDS(file.path(fpath, sample, "analysis", "full.rds"))
    m <- match(prefixes, colnames(full.obj))
    full.obj <- full.obj[,m]

    # Adding in the metadata.
    addition.mode <- character(length(prefixes))
    addition.mode[full.obj$samples$ercc.first] <- "ERCC+SIRV"
    addition.mode[full.obj$samples$sirv.first] <- "SIRV+ERCC"
    addition.mode[full.obj$samples$premixed] <- "Premixed"
    treatment <- "not applicable"
    well.type <- rep("single cell", length(prefixes))
    well.type[full.obj$samples$control=="+"] <- "50 cells"
    well.type[full.obj$samples$control=="-"] <- "empty"

    # Adding in the MD5 sums.
    md5.sums <- read.table(file.path(fpath, sample, "fastq", "md5.all"),
                           header=FALSE, stringsAsFactor=FALSE, comment="")
    md5.sums <- md5.sums[order(gsub("_", "X", md5.sums[,2])),] # weird sorting order with underscores in UTF-8
    m <- match(sub("_[12].fq.gz$", "", md5.sums[,2]), prefixes)

    # Creating links to files.
    curpath <- file.path(fpath, sample, "fastq")
    chosen <- list.files(curpath, pattern="fq.gz$")
    write(file=relink, paste0("ln -s ", file.path("..", curpath, chosen), " ", chosen), append=TRUE, ncol=1)
    new.count.file <- paste0("counts_Liora_", batch, ".tsv")
    write(file=relink, paste0("ln -s ", normalizePath(cpath), " ", new.count.file), append=TRUE, ncol=1)

    # Setting manual standard deviation values.
    if (sample=="test_20170201") {
        ave.frag <- 461
        sd.frag <- 182
    } else {
        ave.frag <- 402
        sd.frag <- 194       
    }

    out <- data.frame(Sample=prefixes, Batch=batch, Addition=addition.mode, Treatment=treatment, 
                      Well=well.type, Counts=new.count.file, MeanFrag=ave.frag, SDFrag=sd.frag)[m,]
    out$File <- md5.sums[,2]
    out$MD5 <- md5.sums[,1]
    collected[[sample]] <- out 
}

collected <- do.call(rbind, collected)

output <- list()
output[["Source Name"]] <- collected$Sample
output[["Characteristics[organism]"]] <- "Mus musculus"
output[["Characteristics[cell line]"]] <- "Trophoblast stem cell"
output[["Characteristics[single cell well quality]"]] <- collected$Well
output[["Material Type"]] <- "RNA"
output[[paste0(rep(c("Protocol REF", "Performer"), 5), collapse="\t")]] <- paste0(c("Obtaining TSCs", "Liora Vilmovsky",
                                                                                    "Culturing TSCs", "Liora Vilmovsky",
                                                                                    "Reverse transcription", "Liora Vilmovsky",
                                                                                    "Extracting RNA", "Liora Vilmovsky",
                                                                                    "Creating libraries","Liora Vilmovsky"
                                                                                    ), collapse="\t")
output[["Extract Name"]] <- collected$Sample
output[["Comment[LIBRARY_LAYOUT]"]] <- "PAIRED"
output[["Comment[LIBRARY_SELECTION]"]] <- "Oligo-dT"
output[["Comment[LIBRARY_SOURCE]"]] <- "TRANSCRIPTOMIC"
output[["Comment[LIBRARY_STRAND]"]] <- "not applicable"
output[["Comment[LIBRARY_STRATEGY]"]] <- "RNA-seq"
output[["Comment[NOMINAL_LENGTH]"]] <- collected$MeanFrag
output[["Comment[NOMINAL_SDEV]"]] <- collected$SDFrag
output[["Comment[ORIENTATION]"]] <- "5'-3'-3'-5'"
output[["Protocol REF\tPerformer"]] <- "Sequencing libraries\tLiora Vilmovsky"
output[["Assay Name"]] <- collected$Sample
output[["Technology Type"]] <- "sequencing assay"
output[["Array Data File"]] <- collected$File
output[["Protocol REF\tPerformer"]] <- "Assigning reads to genes\tAaron Lun"
output[["Derived Array Data File"]] <- collected$Counts
output[["Comment[MD5]"]] <- collected$MD5
output[["Factor Value[spike-in addition]"]] <- collected$Addition
output[["Factor Value[treatment]"]] <- collected$Treatment
output[["Factor Value[block]"]] <- collected$Batch

# Constructing the sdrf.tsv file.

output$check.names <- FALSE
sdrf <- do.call(data.frame, output)
write.table(file="sdrf_Liora.tsv", sdrf, row.names=FALSE, sep="\t", quote=FALSE)

