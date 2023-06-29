source("~/prog/ngs/R scripts/STscripts.R")

library(NNLM)
#library(EnsDb.Hsapiens.v86)
library(survival)
library(xgboost)
library(irlba);
options(mc.cores=24)
library(pROC)
library(Matrix);

# Use regressor to estimate proportion of each annotation in each spot
ist = readRDS(paste0(dataDir, "/classification/classifAll.RDS"))
ist = lapply(ist, function(ii)
{ ii = ii$pr;
  ii = ii[,c('Fat tissue','in situ','Lactiferous duct','Lymphoid nodule','Necrosis',
    "Tumor", "Stroma", "Lymphocyte", 'Vessels')];
})

for (nm in names(ist))
{ message(nm);
  y = readRDS(paste0(dataDir, "counts/TNBC", nm, ".RDS"));
  x = y$cnts;
  
  #x = cntt;
  
  n = ist[[nm]];
  n = n/rowSums(n);
  rownames(x) = rownames(n) = rownames(y$spots);

  n = n[,colSums(n)>1];
  x = x[,colSums(x)>5];

  fm = expressionInAnnotations(x, n, maxIter=5, mccores=TRUE, quiet=FALSE)
  prot = fm$W; rownames(prot) = colnames(x); colnames(prot) = colnames(n);
  Ntot = colSums(n)
  save(prot, fm, Ntot, file=paste0(dataDir, "deconvolution/TNBC", nm, ".RData"));
}