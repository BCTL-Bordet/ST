# ST
Scripts for the analysis of data from the spatial transcriptomics on triple negative breast cancers.

`STstuff` and `nonSTstuff` are two R packages that contain various functions for the ST analysis
(`STstuff`) and other plotting etc functions (`nonSTstuff`).

The ST TNBC folder contains scripts that were used to analyse the data for the ST TNBC paper (submitted).
The data in this paper are available on [Zenodo]https://zenodo.org/records/8135722.

## Installation
The two packages can be installed directly from `R` -
```
library(devtools)
install_github("BCTL-Bordet/ST/STstuff")
install_github("BCTL-Bordet/ST/nonSTstuff")
```
The packages do not need compilation and so should install quickly.

## Software requirement
This software should work on any machine on which R can be installed. It has been tested on macOS
and Linux.