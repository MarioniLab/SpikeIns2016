# Assessing the reliability of spike-ins for scRNA-seq data

To run the real data analyses, enter `real` and follow these instructions:

1. Enter `ArrayExpress`, download the count matrices from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.processed.1.zip) and unpack them. 
Also download the [SDRF file](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.sdrf.txt).
2. Install `package` in your R installation, using `R CMD INSTALL --preclean package`.
3. Enter `Calero/public_I` and run `run_me.sh`.
This will produce a Markdown document containing the variance estimates, along with some serialized R objects for further inspection if necessary.
4. Repeat for `Calero/public_II`, `Liora/public_I` and `Liora/public_II`.
5. Run `make_pics.R` to reproduce the figures in the paper.

To run the simulations, enter the `simulations` directory:

1. Enter `variance` and run `varsim.R` to perform the simulations for detecting HVGs.
2. Enter `diffexp` and run `desim.R` to perform the simulations for detecting DEGs.
This requires you to download a count matrix from [GEO](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29087/suppl/GSE29087%5FL139%5Fexpression%5Ftab%2Etxt%2Egz).
3. Enter `clustering` and run `clustsim.R` to perform the simulations for clustering and PCA.
This requires you to obtain `pancreas_refseq_rpkms_counts_3514sc.txt.gz` from the processed files in [E-MAT-5061](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip), along with the corresponding [SDRF file](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt).
4. Enter `pics/` and run `picmaker.R` and `get_properties.R` to reproduce the figures in the paper.

The `sequence_check/biophysical` directory contains a script to examine the differences in biophysical properties between the spike-in and endogenous mouse genes.

The `manuscript` directory contains all LaTeX code used to generate the manuscript.
This can be compiled with `make`.
