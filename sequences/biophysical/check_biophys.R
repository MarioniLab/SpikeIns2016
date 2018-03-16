#############################################################
# Get the length and GC content of the ERCCs.

library(Biostrings)
epath <- "/lustre/jmlab/resources/genomes/sequences/spikes/ERCC92.fa"
erccs <- readDNAStringSet(epath)

elen <- width(erccs)
edead <- DNAStringSet(sub("A+$", "", as.character(erccs)))
tab <- alphabetFrequency(edead)
egc <- rowSums(tab[,c("G", "C")])/rowSums(tab)

# Get the length and GC content of the SIRVs (these need to be processed into transcripts).

spath <- "/lustre/jmlab/resources/genomes/sequences/spikes/SIRV_150601a.fasta"
raw.sirvs <- readDNAStringSet(spath)

library(rtracklayer)
sirv.gtf <- import("/lustre/jmlab/resources/annotation/original/SIRV_C_150601a.gtf")

library(BSgenome)
sub.sirvs <- getSeq(raw.sirvs, sirv.gtf)
by.transcript <- split(sub.sirvs, sirv.gtf$transcript_id)
sirvs <- DNAStringSet(unlist(lapply(by.transcript, paste0, collapse="")))

slen <- width(sirvs)
sdead <- DNAStringSet(sub("A+$", "", as.character(sirvs)))
tab <- alphabetFrequency(sdead)
sgc <- rowSums(tab[,c("G", "C")])/rowSums(tab)

# Get the length and GC content of endogenous transcripts.

library(BSgenome.Mmusculus.UCSC.mm10)
mm10.gtf <- import("/lustre/jmlab/resources/annotation/processed/mm10.gtf")
chosen.mm10.gtf <- mm10.gtf[grep("^chr", seqnames(mm10.gtf))]
sub.mm10 <- getSeq(BSgenome.Mmusculus.UCSC.mm10, chosen.mm10.gtf)
all.transcript <- split(sub.mm10, chosen.mm10.gtf$transcript_id)

set.seed(1000) # Picking 2000 transcripts to show, otherwise this takes too long.
by.transcript <- all.transcript[sample(length(all.transcript), 2000)]
mm10 <- DNAStringSet(unlist(lapply(by.transcript, paste0, collapse="")))

mlen <- width(mm10)
tab <- alphabetFrequency(mm10)
mgc <- rowSums(tab[,c("G", "C")])/rowSums(tab)

# Making a pretty plot of biophysical characteristics.

pdf("comparison.pdf")
cols <- c("dodgerblue", "salmon", "orange")
boxplot(list(elen/1e3, slen/1e3, mlen/1e3), ylab="Length (kbp)", col=cols, log="y", cex.axis=0.9, cex.lab=1.4, names=rep("", 3))
axis(1, at=1:3, labels=c("ERCC", "SIRV", "mm10"), cex.axis=1.2)
boxplot(list(ERCC=egc*100, SIRV=sgc*100, mm10=mgc*100), ylab="GC content (%)", col=cols, cex.axis=1.2, cex.lab=1.4, cex.names=1.4)
dev.off()

