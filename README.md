# Assessing the reliability of spike-ins for scRNA-seq data

To run the real data analyses, enter `real` and follow these instructions:

1. Enter `ArrayExpress/`, download the count matrices from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.processed.1.zip) and unpack them. 
Also download the [SDRF file](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.sdrf.txt).
2. Install `package/` in your R installation, using `R CMD INSTALL --preclean package`.
3. Enter `Calero/trial_20160113/` and run `run_me.sh`.
This will produce a Markdown document containing the variance estimates, along with some serialized R objects for further inspection if necessary.
4. Repeat for `Calero/trial_20160325`, `Liora/test_20160906` and `Liora/test_20170201`.
5. Run `make_pics.R` to reproduce the figures in the paper.
6. Enter `depth/` and run `runner.R` to perform the simulations for sampling noise.
7. Enter `index_swapping/` and run `check_swap.R` to generate the figures checking for index swapping.

To run the simulations, enter the `simulations` directory:

1. Enter `datasets/` and run `download.sh` to download the various data files. 
Also run `prerunner.R` to pre-process them into serialized R objects.
2. Enter `sampling/` and run `resampler.R` to simulate sampling noise.
You can also run `sizevar.R` to quantify the variation in size factors across cells.
3. Enter `variance/` and run `varsim.R` to perform the simulations for detecting HVGs.
4. Enter `diffexp/` and run `desim.R` to perform the simulations for detecting DEGs.
5. Enter `clustering/` and run `clustsim.R` to perform the simulations for clustering and PCA.
This requires you to obtain `pancreas_refseq_rpkms_counts_3514sc.txt.gz` from the processed files in [E-MAT-5061](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip), along with the corresponding [SDRF file](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt).
6. Enter `pics/` and run `picmaker.R` and `get_properties.R` to reproduce the figures in the paper.

The `sequence_check/biophysical/` directory contains a script to examine the differences in biophysical properties between the spike-in and endogenous mouse genes.

The `manuscript` directory contains all LaTeX code used to generate the manuscript.
This can be compiled with `make`.
