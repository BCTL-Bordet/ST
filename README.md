# ST
Scripts for the analysis of data from the spatial transcriptomics on triple negative breast cancers.

`STstuff` and `nonSTstuff` are two R packages that contain various functions for the ST analysis
(`STstuff`) and other plotting etc functions (`nonSTstuff`).

The ST TNBC folder contains scripts that were used to analyse the data for the ST TNBC paper (submitted).
The data used in this paper are on [Zenodo](https://zenodo.org/records/8135722) and will be
made publicaly available upon acceptance of the manuscript.

## Usage
Both packages include manual pages.

Detailed instructions on how to use the original methods are available in the vignette from the `STstuff` package.
It is available for example from the `R` prompt using
```
vignette("STstuff")
```
A pre-compiled `pdf` version of the vignette is included.
It would take about 1.5 minutes to recompile on a desktop computer.

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
and Linux with R 4.2.1.

The `STstuff` package depends on `matrixStats`, `FNN` and `xgboost`.
The `nonSTstuff` package depends on `matrixStats`.
The scripts for the paper further depend on `RColorBrewer`, `corrplot`, `EnsDb.Hsapiens.v86`,
  `EBImage`, `RANN`, `mclust`, `umap`, `xCell`, `irlba`, `pROC`, `reldist`, `GSVA`, `lme4`,
  `TeachingDemos`, `xlsx`, `alluvial`, `NbClust`, `jpeg`, `Polychrome`, `glmGamPoi`, `scran`,
  `Matrix`, `rdist`.