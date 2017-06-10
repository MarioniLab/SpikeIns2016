ln -s ../../ArrayExpress/counts_Liora_20170201.tsv genic_counts.tsv
ln -s ../../ArrayExpress/E-MTAB-5522.sdrf.txt sdrf.tsv
echo "knitr::knit('public_Liora.Rmd')" | R --no-save 
