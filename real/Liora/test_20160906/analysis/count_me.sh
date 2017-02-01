echo '
anno.files <- c("/lustre/jmlab/resources/annotation/processed/mm10.gtf", 
	"/lustre/jmlab/resources/annotation/original/ERCC92.gtf", 
	"/lustre/jmlab/resources/annotation/original/SIRV_C_150601a.gtf")
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
stat.file <- "../all_qual.tsv"
ispet <- TRUE
# NOT STRAND-SPECIFIC.
' | cat - ~/Code/mapping/counter.R > count_me.R

R CMD BATCH --no-save count_me.R 
