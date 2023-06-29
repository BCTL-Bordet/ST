# Scripts used for the analyses in the ST TNBC paper

First, the two packages STstuff and nonSTstuff must be installed.

Then the STscript.R file must be sourced. The first line of that file specify where the data 

After this, it is possible to reproduce the analysis done.

## Install.R
Extract the data from the tables as RDS/RData, extract the annotation information.
Also perform the batch effect correction between slides from the same patient.

## Regression.R
Estimates from the gene expression in each spot the prevalence of each annotation.
Does the validation (leave a patient out) and apply on each sample.

## Deconvolution.R
Uses theannotaiton prevalenc eobtained before to deconvolute the gene expression relative to each annottation in each sample.

## MC ET .R
Does the clustering of spots intra patient, then the megaclustering of chose clustering across patients, and finally determine ecotypes om, that is gropus of patients that share similar megacluster disctribution.s

