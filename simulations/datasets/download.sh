# Kolod et al.
wget http://www.ebi.ac.uk/teichmann-srv/espresso/static/counttable_es.csv

# Buettner et al.
mkdir E-MTAB-2805
cd E-MTAB-2805
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2805/E-MTAB-2805.processed.1.zip 
unzip E-MTAB-2805.processed.1.zip
rm E-MTAB-2805.processed.1.zip
cd -

# Zeisel et al.
wget https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
wget https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt

# Islam et al. (1)
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29087/suppl/GSE29087%5FL139%5Fexpression%5Ftab%2Etxt%2Egz

# Islam et al. (2)
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46980/suppl/GSE46980%5FCombinedMoleculeCounts%2Etab%2Egz

# Grun et al.
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54695/suppl/GSE54695%5Fdata%5Ftranscript%5Fcounts%2Etxt%2Egz

# Wilson et al.
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61533/suppl/GSE61533%5FHTSEQ%5Fcount%5Fresults%2Exls%2Egz

# Hashimshony et al.
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78779/suppl/GSE78779%5FExpression%5FC1%5F96%5Fcells%2Etxt%2Egz

# Scialdone et al.
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3707/E-MTAB-3707.processed.1.zip
unzip E-MTAB-3707.processed.1.zip
rm E-MTAB-3707.processed.1.zip
