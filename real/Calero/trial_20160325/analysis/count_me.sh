echo '
anno.files <- c("~/lustre/annotation/mm10.gtf", "~/lustre/annotation/ERCC92.gtf", "~/lustre/annotation/SIRV_C_150601a.gtf", "../../genomes/sequences/oncogene.gtf")
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
stat.file <- "../all_qual.tsv"
ispet <- FALSE
' | cat - ~/Code/mapping/counter.R > count_me.R

R CMD BATCH --no-save count_me.R 
