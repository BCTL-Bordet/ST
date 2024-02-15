# Scripts used for the analyses in the ST TNBC paper

First, the two packages `STstuff` and `nonSTstuff` must be installed.

Then the `STscript.R` file must be sourced.
The first line of that file specify the directory where the data is (in  `dataDir`).
The data is supposed to correspond to the one on [Zenodo](https://zenodo.org/records/8135722).
The `.tar.gz` should be decompressed in that directory, to obtain the expected directory structure.

After this, it is possible to reproduce the analysis done.

## `install.R`
Extract the data from the tables as RDS/RData, extract the annotation information.
Also perform the batch effect correction between slides from the same patient.

## `morpho.R`
Calculate the distribution of patches sizes for different annotations.

## `regression.R`
Estimates from the gene expression in each spot the prevalence of each annotation.
Does the validation (leave a patient out) and apply on each sample.

## `deconvolution.R`
Uses the annotation prevalence obtained before to deconvolute the gene expression relative to
each annottation in each sample.

## `MC ET.R`
Does the clustering of spots intra patient, then the megaclustering of chose clustering across patients,
and finally determine ecotypes, that is groups of patients that share similar megacluster distributions.

## `ST TNBC start.R`
A script to load all the data (regressions, megaclusters, bulk, pseudo-bulk,...)
to be able to perform final analyses.
Must be sourced before running anything in `ST TNBC figs.R`.


## `ST TNBC figs.R`
Code to generate the figures from the paper.
Note that that the `figDir` variable, at the top of the file, is the directory where the
figures will be plotted.