
anno.files <- c("~/lustre/annotation/mm10.gtf", "~/lustre/annotation/ERCC92.gtf", "~/lustre/annotation/SIRV_C_150601a.gtf", "../../genomes/sequences/oncogene.gtf")
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
stat.file <- "../all_qual.tsv"
ispet <- FALSE

# Files of interest.
bam.files
anno.files
stat.file

if (!exists("ispet")) { 
    ispet <- FALSE
} 
if (!exists("strandspec")) {
    strandspec <- 0
}

ispet
strandspec

if (length(anno.files)==1L) {
    file.symlink(anno.files, "temp.gtf")
} else {
    system(paste(c("cat", anno.files, "> temp.gtf"), collapse=" "))
}

# Running featureCounts.
require(Rsubread)
out <- featureCounts(bam.files, annot.ext="temp.gtf", isGTFAnnotationFile=TRUE, minMQS=10, nthreads=4, isPairedEnd=ispet, strandSpecific=strandspec)

# Saving counts to file, with gene names.
colnames(out$counts) <- sub("\\.bam$", "", basename(bam.files))
final <- data.frame(GeneID=rownames(out$counts), Length=out$annotation$Length, out$counts)
write.table(file="genic_counts.tsv", final, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Augmenting the stats.
my.stats <- read.table(stat.file, header=TRUE)
m <- match(my.stats$Sample, colnames(out$counts))
my.stats$Genic <- as.integer(out$stat[1,-1][m])
write.table(file="my_qual.tsv", my.stats, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Saving the session information.
unlink("temp.gtf")
sessionInfo()

