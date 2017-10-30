ln -s ../../../ArrayExpress/counts_Calero_20160325.tsv genic_counts.tsv
ln -s ../../../ArrayExpress/E-MTAB-5522.sdrf.txt sdrf.tsv
echo ".spike.local <- FALSE; knitr::knit('techanal.Rmd')" | R --no-save 

