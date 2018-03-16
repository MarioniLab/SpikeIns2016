echo '
anno.files <- file.path("../../../../sequences/annotation/",
    c("mm10.gtf", "ERCC92.gtf", "SIRV_C_150601a.gtf"))
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
stat.file <- "../all_qual.tsv"
ispet <- TRUE
# NOT STRAND-SPECIFIC.
' | cat - ~/Code/mapping/counter.R > count_me.R

R CMD BATCH --no-save count_me.R 
