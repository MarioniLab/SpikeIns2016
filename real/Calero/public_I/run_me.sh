ln -s ../../ArrayExpress/counts_Calero_20160113.tsv genic_counts.tsv
ln -s ../../ArrayExpress/sdrf.tsv sdrf.tsv
echo "knitr::knit('public_Calero.Rmd')" | R --no-save 
