ln -s ../../ArrayExpress/counts_Liora_20160906.tsv genic_counts.tsv
ln -s ../../ArrayExpress/E-MTAB-5522.sdrf.txt sdrf.tsv
echo "knitr::knit('public_Liora.Rmd')" | R --no-save 

# Reconstructing as if we had run it in the original directory.
curwd=$(pwd)
cd ../test_20160906/analysis/ 
ln -s ${curwd}/*.rds . 
