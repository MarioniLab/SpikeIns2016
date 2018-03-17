# Downloading all the necessary materials.
cd real/ArrayExpress
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.processed.1.zip
unzip E-MTAB-5522.processed.1.zip
rm E-MTAB-5522.processed.1.zip
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.sdrf.txt

# Moving them to the correct locations.
cd ../
mv ArrayExpress/counts_Calero_20160113.tsv Calero/trial_20160113/analysis/genic_counts.tsv
mv ArrayExpress/counts_Calero_20160325.tsv Calero/trial_20160325/analysis/genic_counts.tsv
mv ArrayExpress/counts_Liora_20160906.tsv Liora/test_20160906/analysis/genic_counts.tsv
mv ArrayExpress/counts_Liora_20170201.tsv Liora/test_20170201/analysis/genic_counts.tsv

# Also cloning the tools (using the last commit known to work).
cd ../
git clone https://github.com/LTLA/CRUKTools tools
cd tools
git reset --hard 915706ac016f6f59e74b4ec266f06349293b997f
