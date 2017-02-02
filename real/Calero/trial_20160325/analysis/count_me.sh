echo '
anno.files <- file.path("/lustre/jmlab/resources/annotation", 
    c("processed/mm10.gtf", "original/ERCC92.gtf", "original/SIRV_C_150601a.gtf", "original/CBFB-MYH11-mcherry.gtf"))
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
stat.file <- "../all_qual.tsv"
ispet <- FALSE
# NOT STRAND-SPECIFIC.
' | cat - ~/Code/mapping/counter.R > count_me.R

R CMD BATCH --no-save count_me.R 
