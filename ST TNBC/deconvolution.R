source("~/prog/ST/ST TNBC/STscripts.R")

# Use regressor to estimate proportion of each annotation in each spot
ist = readRDS(paste0(dataDir, "/classification/classifAll.RDS"))
ist = lapply(ist, function(ii)
{ ii = ii$pr;
  ii = ii[,c('Fat tissue','in situ','Lactiferous duct','Lymphoid nodule','Necrosis',
    "Tumor", "Stroma", "Lymphocyte", 'Vessels')];
})

for (nm in names(ist))
{ message(nm);
  y = readRDS(paste0(dataDir, "Robjects/counts/TNBC", nm, ".RDS"));
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

# TLS Signature
##################
ids = readRDS(paste0(dataDir, "Clinical/ids.RDS"))
cli = readRDS(paste0(dataDir, "Clinical/Clinical.RDS"))
geneMap = readRDS(paste0(dataDir, "misc/geneMap.RDS"))
# Select only the samples with actual TLS
hasTLS = c(12,15,16,19,2,20,27,28,29,30,31,32,33,38,39,42,43,44,45,46,48,5,50,52,53,54,56,57,6,61,62,63,64,65,67,69,70,71,72,
  73,77,78,80,82,85,86,87,88,89,9,92,94)
cli$hasTLS = rownames(cli) %in% hasTLS;

# 1. Load deconvolution data (Same code in "ST TNBC start.R")
allS = mclapply(d<-dir(paste0(dataDir, "deconvolution")), function(i)
{ load(paste0(dataDir, "deconvolution/", i))
  list(prot=rowsum(prot, geneMap$PB[rownames(prot)]), Ntot=Ntot);
}); names(allS) = sub("TNBC([0-9]+).RData", "\\1", d); allS=allS[rownames(cli)];
Ntot = unlist(lapply(allS, function(i) i$Ntot))
allS = lapply(allS, function(i) i$prot);

g = table(unlist(lapply(allS, function(i) rownames(i)))); g = names(g)[g>=max(g)/5]
xCanRaw = fullMat3(allS, g);
idCan = cbind(what=colnames(xCanRaw), id=rep(names(allS), sapply(allS, ncol)))
colnames(xCanRaw) = paste(idCan[,1], idCan[,2]);
Nt = matrix(NA, nrow=length(unique(idCan[,2])), ncol=length(unique(idCan[,1])));
rownames(Nt) = rownames(cli); colnames(Nt) = unique(idCan[,1]);
Nt[idCan[,c(2,1)]] = Ntot;

xCan = normIt(xCanRaw);

idAn = matrix(NA, nrow=length(unique(idCan[,2])), ncol=length(unique(idCan[,1])));
dimnames(idAn) = dimnames(Nt);
idAn[idCan[,c(2,1)]] = 1:nrow(idCan);

# 2. Make comparisons
tt = setdiff(colnames(idAn), "Lymphoid nodule"); # Comparisons to do

pinp = lapply(tt, function(z) # Not paired
{ a = xCan[,idAn[cli$hasTLS,"Lymphoid nodule"]]; b = xCan[, idAn[,z]];
  a = a[,!is.na(a[1,])]; b = b[,!is.na(b[1,])]
  p = unlist(mclapply(1:nrow(a), function(i) wilcox.test(a[i,], b[i,], alt='g')$p.value))
  m = rowMedians(a)-rowMedians(b);
  cbind(m=m, p=p);
}); names(pinp) = tt;

pi = pinp;
ps = do.call(cbind, lapply(pi, function(i) i[,"p"]));
ms = do.call(cbind, lapply(pi, function(i) i[,"m"]));
rownames(ps) = rownames(ms) = rownames(xCan);
p = rowMaxs(ps); m = rowMins(ms); names(p) = names(m) = rownames(ps);

# And finally the signature
gnsTLS = names(p)[p<1e-9 & m>log(2)]