# Assessing the reliability of spike-ins for scRNA-seq data

## Overview

This repository provides the code (and manuscript files!) for the paper **Assessing the reliability of spike-in normalization for analyses of single-cell RNA sequencing data**
by Lun _et al._, published in _Genome Research_ in 2017.

## Repeating the real data analysis

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

## Repeating the simulations

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

The `sequences/biophysical/` directory contains a script to examine the differences in biophysical properties between the spike-in and endogenous mouse genes.

## Realigning and recounting

Readers interested in regenerating the count matrices from the FASTQ files are advised to:

1. Obtain master scripts from https://github.com/LTLA/CRUKTools (last tested with 915706ac016f6f59e74b4ec266f06349293b997f),
in addition to the _subread_ aligner and _Rsubread_.
2. Follow the instructions in `sequences/genomes/README.md` to build the genome indices.
Similarly, follow the instructions in `sequences/annotation/README.md` to obtain the annotation.
3. Download the FASTQ files from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5522/).
Files corresponding to each batch of data should be placed in `Calero/trial_20160113/fastq`, etc.
4. Run the various `mapme.sh` scripts to execute the master scripts for alignment, and `count_me.sh` for read counting.
Obviously the paths may need to be changed to point to the location of the master scripts.

## Creating the manuscript

The `manuscript` directory contains all LaTeX code used to generate the manuscript.
This can be compiled with `make`.


