source("~/prog/ngs/R scripts/STscripts.R")

library(NNLM)
#library(EnsDb.Hsapiens.v86)
library(survival)
library(xgboost)
library(irlba);
options(mc.cores=4)
library(pROC)
library(Matrix);

dataDir = "~/Data/Spatial/TNBC/datadist/"
ids = readRDS(paste0(dataDir, "Clinical/ids.RDS"))

# Install (Note: need a lot of RAM!)
#####################################
nm = ids$id[which(ids$hasAnnot)]
# Read annotations
annots = lapply(nm, function(i) readRDS(paste0(dataDir, "annotsBySpot/TNBC", i, ".RDS"))[[1]]); names(annots)=nm;

# 1. Read counts with annotations
cnt = mclapply(nm, function(i)
{ #message(i);
  x = readRDS(paste0(dataDir, "counts/TNBC", i, ".RDS"));
  w = rownames(ids)[ids$id==i & ids$hasAnnot];
  r = Matrix(t(x$cnts[x$spots$slide==w,]));
}, mc.cores=4); names(cnt) = nm;

g = unique(unlist(lapply(cnt, rownames)));
for (i in seq_along(cnt)) { cnt[[i]] = fullMatSparse(cnt[[i]], g) }

m = do.call(cbind, mclapply(cnt, function(x) apply(x,1,max)))
m = rowMaxs(m)>=5;
for (i in seq_along(cnt)) { cnt[[i]] = cnt[[i]][m,]; }

# 2. Make Big matrix
cnts = do.call(cbind, lapply(cnt, as.matrix))
pts = rep(names(cnt), sapply(cnt, function(i) ncol(i)))

# Full PCA
x = makePR(t(cnts), Ng=4e3);
xx = (x-rep(m<-colMeans(x), each=nrow(x)))*rep(1/colSds(x, center=m), each=nrow(x));
prc = prcomp_irlba(xx, n=500, center=FALSE, scale.=FALSE,  maxit = 10000, work=1e3);
pr = prc$x

rot = prc$rotation; rownames(rot) = colnames(x)

save(pts, pr, annots, file=paste0(dataDir, "classification/LOOstart.RData"))
saveRDS(list(rot=rot, center=colMeans(x), scale=colSds(x), g=rownames(cnts)),
  file=paste0(dataDir, "classification/baseClassif.RDS"))
save(pr, xx, cntsT, file=paste0(dataDir, "classification/classifData.RData"))

# Project all samples into the PC and save
bc = readRDS(paste0(dataDir, "classification/baseClassif.RDS"))
f = mclapply(dir(paste0(dataDir, "counts")), function(nmf)
{ x = readRDS(paste0(dataDir, "counts/", nmf));
  z = rescaleForClassif(bc, t(x$cnts));
  spot = x$spots;
  save(z, spot, file=paste0(dataDir, "classification/projectedSamples/", sub(".RDS", ".RData", nmf)))
}, mc.preschedule=FALSE, mc.cores=12);

# LPO X-validation
######################
load(paste0(dataDir, "classification/LOOstart.RData"))

n = do.call(rbind, annots);
n = cbind(n, Stroma=NA);
n = n[,!(colnames(n) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))]
n[,"Tumor"] = n[,"Tumor"] + n[,"Tumor region"];
n[,"Stroma"] = n[,"Stroma cell"] + n[,"Low TIL stroma"] + n[,"High TIL stroma"]/2 + n[,"Acellular stroma"];
n[,"Lymphocyte"] = n[,"Lymphocyte"] + n[,"High TIL stroma"]/2
n = n[,!(colnames(n) %in% c("Low TIL stroma", "High TIL stroma", "Stroma cell", "Acellular stroma",
  "Tumor region"))];
  
pr2 = pr; pts2 = pts;
w = rowSums(n)>1000; n=n[w,]; pr2=pr2[w,]; pts2=pts2[w];
n = n/rowSums(n);
prw = pr2[,1:250];

loos=list();
for (what in colnames(n))
{ message(what);
  loo = lapply(split(seq_along(pts2), pts2), function(i)
  { fm = xgboost(prw[-i,], n[-i,what], nrounds=25, booster = "gblinear", objective = "reg:logistic", verbose=0,
      nthread=options("mc.cores"))
    r = predict(fm, newdata=prw[i,], nthread=options("mc.cores"))
  })
  loos[[what]] = unsplit(loo, pts2);
}
save(loos, n, file=paste0(dataDir, "classification/loosReg.RData"));

# Make regressors (note: starts with n etc from before)
########################################
for (what in colnames(n))
{ fm = xgboost(prw, n[,what], nrounds=25, booster = "gblinear", objective = "reg:logistic", verbose=0,
        nthread=options("mc.cores"))
  xgb.save(fm, fname=path.expand(paste0(dataDir,"classification/regressors/", what, "Reg.xgb")))
}

# Apply on all
###############
d = dir(paste0(dataDir,"classification/regressors/"), pattern=".xgb");
fmGlob = lapply(d, function(what) xgb.load(path.expand(paste0(dataDir,"classification/regressors/", what))));
names(fmGlob) = sub("Reg.xgb", "", d);

d = dir(paste0(dataDir, "classification/projectedSamples/"))
ist = lapply(d, function(nmf)
{ message(nmf);
  load(paste0(dataDir, "classification/projectedSamples/", nmf));
  ist = do.call(cbind, lapply(fmGlob, function(fm) predict(fm, newdata=z[,1:250])));
  return(list(xy=xy, pr=ist, spot=spot));
}); names(ist) = sub("TNBC([0-9]+).RData", "\\1", d);
saveRDS(ist, file=paste0(dataDir, "classification/classifAll.RDS"))