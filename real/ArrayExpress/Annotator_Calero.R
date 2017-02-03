########################################################################################
# Print out important bits from the metadata for Fernando's lot.

collected <- list()
fpath <- "/run/user/1753941046/gvfs/smb-share:server=jmlab-data,share=jmlab/group_folders/lun01/Internal/SpikeIns"
relink <- "make_links_Calero.sh"
write(file=relink, c("set -e", "set -u", "mkdir fastq_Calero"), ncol=1)


for (sample in c("Calero/trial_20160113", "Calero/trial_20160325")) {
    cpath <- file.path("..", sample, "analysis", "genic_counts.tsv")
    all.files <- read.table(cpath, nrows=1, stringsAsFactor=FALSE)
    prefixes <- as.character(all.files[-c(1:2)])
    batch <- sub(".*_", "", basename(sample))

    # Loading in the corresponding object.
    full.obj <- readRDS(file.path("..", sample, "analysis", "full.rds"))
    m <- match(prefixes, colnames(full.obj))
    full.obj <- full.obj[,m]

    # Adding in the metadata.
    addition.mode <- character(length(prefixes))
    addition.mode[full.obj$samples$ercc.first] <- "ERCC+SIRV"
    addition.mode[full.obj$samples$sirv.first] <- "SIRV+ERCC"
    addition.mode[full.obj$samples$premixed] <- "Premixed"
    treatment <- ifelse(full.obj$samples$induced, "Induced", "Control")

    # Adding in the MD5 sums.
    fnames <- paste0(sub("SLX\\.", "SLX-", prefixes), ".fq.gz")
    md5.sums <- read.table(file.path(fpath, sample, "fastq", "md5.all"),
                           header=FALSE, stringsAsFactor=FALSE)
    m <- match(fnames, md5.sums[,2])
    md5.sums <- md5.sums[m,1]

    # Creating links to files.
    curpath <- file.path(fpath, sample, "fastq")
    chosen <- list.files(curpath, pattern="fq.gz$")
    write(file=relink, paste0("ln -s ", file.path(curpath, chosen), " ", file.path("fastq", chosen)), append=TRUE, ncol=1)
    new.count.file <- paste0("counts_Calero_", batch, ".tsv")
    write(file=relink, paste0("ln -s ", normalizePath(cpath), " ", file.path("fastq", new.count.file)), append=TRUE, ncol=1)

    collected[[sample]] <- data.frame(Sample=prefixes, Batch=batch, Addition=addition.mode, Treatment=treatment,
                                      File=fnames, MD5=md5.sums, Counts=new.count.file)
}

collected <- do.call(rbind, collected)


output <- list()
output[["Source Name"]] <- collected$Sample
output[["Characteristics[organism]"]] <- "Mus musculus"
output[["Characteristics[cell line]"]] <- "416B"
output[["Material Type"]] <- "RNA"
output[[paste0(rep(c("Protocol REF", "Performer"), 6), collapse="\t")]] <- paste0(c("Obtaining 416B cells", "Fernando Calero",
                                                                                    "Culturing 416B cells", "Fernando Calero",
                                                                                    "Reverse transcription", "Fernando Calero",
                                                                                    "Extracting RNA", "Fernando Calero",
                                                                                    "Creating libraries","Fernando Calero"
                                                                                    ), collapse="\t")
output[["Extract Name"]] <- collected$Sample
output[["Comment[LIBRARY_LAYOUT]"]] <- "SINGLE"
output[["Comment[LIBRARY_SELECTION]"]] <- "Oligo-dT"
output[["Comment[LIBRARY_SOURCE]"]] <- "TRANSCRIPTOMIC"
output[["Comment[LIBRARY_STRAND]"]] <- "not applicable"
output[["Comment[LIBRARY_STRATEGY]"]] <- "RNA-seq"
output[["Comment[NOMINAL_LENGTH]"]] <- "not applicable"
output[["Comment[NOMINAL_SDEV]"]] <- "not applicable"
output[["Comment[ORIENTATION]"]] <- "not applicable"
output[["Protocol REF\tPerformer"]] <- "Sequencing libraries\tFernando Calero"
output[["Assay Name"]] <- collected$Sample
output[["Technology Type"]] <- "sequencing assay"
output[["Comment[experiment batch]"]] <- collected$Batch
output[["Array Data File"]] <- collected$File
output[["Protocol REF"]] <- "Assigning reads to genes"
output[["Derived Array Data File"]] <- collected$Counts
output[["Comment[MD5]"]] <- collected$MD5
output[["Characteristics[single cell well quality]"]] <- "single cell"
output[["Factor Value[spike-in addition]"]] <- collected$Addition
output[["Factor Value[treatment]"]] <- collected$Treatment

# Constructing the sdrf.tsv file.

output$check.names <- FALSE
sdrf <- do.call(data.frame, output)
write.table(file="sdrf_Calero.tsv", sdrf, row.names=FALSE, sep="\t", quote=FALSE)


