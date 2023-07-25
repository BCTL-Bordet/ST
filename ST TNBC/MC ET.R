source("~/prog/ST/ST TNBC/STscripts.R")

# Some data...
ids = readRDS(paste0(dataDir, "Clinical/ids.RDS"))
cli = readRDS(paste0(dataDir, "Clinical/Clinical.RDS"))
geneMap = readRDS(paste0(dataDir, "misc/geneMap.RDS"))
PBraw = readRDS(paste0(dataDir, "PB_count.RDS"))
PBraw = rowsum(PBraw, geneMap$PB[rownames(PBraw)]);
PBraw = PBraw[,rownames(cli)];
PB = normRNAseq(PBraw, lim=1e3)
genePB = rownames(PB)

cli$barPB = factBareche(TNBCclassif(PB, version="bareche", shortName=TRUE, coef=FALSE))

ist = readRDS(paste0(dataDir, "classification/classifAll.RDS"))
ist = lapply(ist, function(ii)
{ ii = ii$pr;
  ii = ii[,c('Fat tissue','in situ','Lactiferous duct','Lymphoid nodule','Necrosis',
    "Tumor", "Stroma", "Lymphocyte", 'Vessels')];
})
ist = ist[rownames(cli)] # Classification per spot

# Intra-patient clustering
############################
library(FNN)
library(umap)

for (nm in rownames(cli))
{ message(nm);
  x = readRDS(paste0(dataDir, "Robjects/counts/TNBC", nm, ".RDS"));
  cntt = x$cnts; spot = x$spot[,c("pixel_x", "pixel_y")];

  clusts = sp = list();

  c2 = cntt[, colMeans(cntt>2)>.05];
  c2 = c2/rowSums(c2)*1e4
  c2 = log(1+c2)
  c2 = apply(c2,2,scale)

  nei = get.knn(spot, k=10)
  scg = apply(c2, 2, function(i) { mean((i - i[nei$nn.index])^2 ) })

  tst = seq(from=100, to=min(500, ncol(c2)), by=100); o = order(scg);

  kms = mclapply(tst, function(lim)   
  { xx = c2[,o[1:lim]];
    tsn2 = umap(xx)$layout;

    km = lapply(2:10, function(k)
    { km1 = mkmeans(xx, k, minNstart=100, nstart=1e4, iter.max=50)
      km2 = mkmeans(tsn2, k, minNstart=100, nstart=1e4, iter.max=50)
      return(list(x=km1$cluster, tsn=km2$cluster))
    });
    return(list(tsn=tsn2, km=km)) 
  });
  qqt = unlist(lapply(kms, function(i) unlist(i$km, recursive=FALSE)), recursive=FALSE);
  #qv = sapply(qqt, function(i) if (any(table(i)<3)) { return(0) } else { mean(i == i[nei$nn.index])-1/max(i) })
  qv = sapply(qqt, function(i) if (any(table(i)<3)) { return(0) } else { mean(i == i[nei$nn.index]) })
  N = sapply(qqt, function(i) max(i))
  id = tapply(seq_along(N), N, function(i) i[which.max(qv[i])])
  #qv = qv-1/N;
  km = do.call(cbind, qqt[id]);#qqt[[wm<-which.max(qv)]]
  qv = qv[id];
  N = N[id];

  clusts = lapply(1:ncol(km), function(x) do.call(cbind, tapply(1:nrow(km), km[,x], function(i) colSums(cntt[i,]))))

  STgenes = scg[o];
  save(clusts, km, qv, N, scg, file=paste0(dataDir, "/clustering/intraPatientClust/TNBC",nm, ".RData"))
}

# Deconvoluted prototypes
###########################
ers = lapply(rownames(cli), function(nm) 
{ message(nm);
  load(paste0(dataDir, "/clustering/intraPatientClust/TNBC", nm , ".RData"))
  cntt = readRDS(paste0(dataDir, "Robjects/counts/TNBC", nm, ".RDS"))$cnts;
  bst = which.max(qv-1/N)
  k = as.integer(km[,bst]);
  cntt = cntt[,intersect(colnames(cntt), genePB)]; cntt = cntt[,colSums(cntt)>0];
  tmp = deconvoluteClusters(cntt, k, anchor=1e3, rmCutOff=-Inf, mccores=TRUE, Niter=40, clean=FALSE)
  tmp$prevProt = NULL; tmp$k = k;
  saveRDS(tmp, file=paste0(dataDir, "clustering/clustPrototypes/TNBC", nm , ".RDS"))
})

# Load full proto
######################
annots = lapply(rownames(cli), function(i)
{ x = readRDS(paste0(dataDir, "Robjects/annotsBySpot/TNBC", i, ".RDS"))
  r = x$annots;
  rownames(r) = paste(rownames(ids)[ids$hasAnnot & ids$id==i], rownames(r));
  return(r);
}); names(annots) = rownames(cli);

annotsClean = lapply(annots, function(n)
{ n = cbind(n, Stroma=NA);
  n = n[,!(colnames(n) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))]
  n[,"Tumor"] = n[,"Tumor"] + n[,"Tumor region"];
  n[,"Stroma"] = n[,"Stroma cell"] + n[,"Low TIL stroma"] + n[,"High TIL stroma"]/2 + n[,"Acellular stroma"];
  n[,"Lymphocyte"] = n[,"Lymphocyte"] + n[,"High TIL stroma"]/2
  n = n[,!(colnames(n) %in% c("Low TIL stroma", "High TIL stroma", "Stroma cell", "Acellular stroma",
    "Tumor region"))];
  n = n/rowSums(n);
}); # Annotations per spot

allClust = mclapply(d<-dir(paste0(dataDir, "clustering/clustPrototypes/")), function(i)
{ n = sub("TNBC([0-9]+).RDS", "\\1", i)
  x = readRDS(paste0(dataDir, "clustering/clustPrototypes/", i));
  if (!identical(x$k, x$kOrig)) { browser(); }
  colnames(x$proto$proto) = x$fit$cl;
  annot = rowsum(ist[[n]], x$k)/tabulate(x$k);
  annot = annot[x$fit$cl,];
  annot2 = do.call(rbind, tapply(1:nrow(annotsClean[[n]]), factor(x$k)[1:nrow(annotsClean[[n]])],
    function(i) colMeans(annotsClean[[n]][i,], na.rm=TRUE)));
  if (length(a<-setdiff(as.character(x$fit$cl), rownames(annot2)))>0)
  { m = matrix(NA, nrow=length(a), ncol=ncol(annot2)); rownames(m)=a; annot2 = rbind(annot2, m); } 
  annot2 = annot2[as.character(x$fit$cl),];
  return(list(proto=x$proto, annot=annot, annot2=annot2,
    nReads=cbind(orig=x$Nreads, deconv=colSums(x$fit$mix)*colSums(x$fit$proto), origSpot=tabulate(x$k),
    deconvSpot=colSums(x$fit$mix))));
}); names(allClust) = sub("TNBC([0-9]+).RDS", "\\1", d);

g = unique(unlist(lapply(allClust, function(i) rownames(i$proto$proto))));
xC = fullMat3(lapply(allClust, function(i) i$proto$proto), g);
colnames(xC) = paste(rep(names(allClust), sapply(allClust, function(i) ncol(i$proto$proto))), colnames(xC))
idC = as.data.frame(do.call(rbind, strsplit(colnames(xC), " ")));
if (ncol(idC)==2) { colnames(idC) = c("id", "n") } else { colnames(idC) = c("id", "tumor", "n"); }
idC$barPB = cli[idC[,'id'], "barPB"];
sigma = fullMat2(lapply(allClust, function(i) i$proto$sigma), g)+.2;
annot = do.call(rbind, lapply(allClust, function(i) i$annot));
annot2 = do.call(rbind, lapply(allClust, function(i) i$annot2));
Nr = do.call(rbind, lapply(allClust, function(x) x$nReads));

# Merge those below Nmin reads
Nmin = 5e4;
wg = order(-rowSums(xC))[1:5e3];
o = order(Nr[,2]); klink = 1:ncol(xC);
a = xC[wg,]; a[a<.1]=.1;
for (j in seq_along(o))
{ i = o[j];
  if (Nr[i,2]>5e4) { next; }
  w2 = setdiff(which(idC[,"id"]==idC[i,"id"]), o[1:j])
  if (length(w2)==0) { next; }
  wb = which.max(cor(a[,i], a[,w2], method='s'));
  xC[,w2[wb]] = xC[,w2[wb]] + xC[,i]; Nr[w2[wb],] = Nr[w2[wb],] + Nr[i,];
  a = xC[wg,]; a[a<.1]=.1;
  klink[klink==i] = w2[wb]
  o[j]=-o[j];
}
names(klink) = paste(idC$id, idC$n);
w = -o[o<0]
xC = xC[,-w]; idC = idC[-w,];  Nr = Nr[-w,]; annot=annot[-w,]; annot2=annot2[-w,]
a = klink; for (i in sort(w, dec=TRUE)) { a[klink>i]=a[klink>i]-1; }
klink=a;

xCn = log(1+1e4*xC/rep(colSums(xC), each=nrow(xC)))

# Add Bareche subtype per cluster
y2 = rowsum(exp(xCn)-1, geneMap$PB[rownames(xCn)]);

g = intersect(rownames(y2), rownames(PBraw));

pb = log(1e4*PBraw[g,]/rep(colSums(PBraw[g,]), each=nrow(PBraw[g,]))+1)
y3 = log(1e4*y2[g,]/rep(colSums(y2[g,]), each=nrow(y2[g,]))+1)

y3n = (y3-rowMeans(pb))/rowSds(pb);
idC$bar = factor(TNBCclassif(y3n, version='bar', shortName=TRUE, rescale=FALSE), levels=names(colBar));

# Kmeans on prototypes
#######################
for (i in 5:20)
{ message(i);
  km = mkmeans(t(y), i, nstart=1e7, fact=100)
  saveRDS(km, file=paste0(dataDir, "clustering/Kmeans MC/km", i, ".RDS"));
}

# LOO on MC
#############
pFromW = function(W, pr=pres)
{ W = W/rowSums(W); W[W<1e-2] = 0;
  p = sapply(1:nrow(pr), function(i) wilcox.test(W[,i]~pr[i,], alt='l')$p.value);
  au = sapply(1:nrow(pr), function(i) auc(pr[i,],  W[,i], direction="<", quiet=TRUE));
  cbind(p=p, auc=au);
}

kms = lapply(d<-dir(fd<-paste0(dataDir, "clustering/Kmeans MC/")), function(n) readRDS(paste0(fd, n)))
names(kms) = sub("km([0-9]+).RDS", "\\1", d); kms = kms[order(as.integer(names(kms)))];

bu = readRDS(paste0(dataDir, "bulk_count.RDS"))
bu = bu[rowSums(bu)>100,];

x = 1e4*xC/rep(colSums(xC), each=nrow(xC))
w = rowMeans(x>1)>.05;
si = sigma[w,];
y = x[w,];
dataForMC = log(y+1); rm(y);

XF = 5;
for (nk in names(kms))
{ message("Doing ", nk);
  k = assignClust(kms[[nk]], t(y), idC, annot, cli, clLim=0, minPatient=1, rmProbOnly=TRUE)
  #k = tmp$k;
  pres = table(k, idC[,1])>0; pres=pres[,rownames(cli)]
  mk = max(k); tk = tabulate(k, nbins=mk);
  sp = unlist(lapply(1:XF, function(iii)
  { for (iter in 1:10000) # Check if OK
    { a = split(rownames(cli), floor(10*sample(0:(nrow(cli)-1))/nrow(cli)))
      if (any(sapply(a, function(i) { min(tk-tabulate(k[idC[,1] %in% i], nbins=mk)) })==0)) { next; }
      return(a);
    }
    stop("No possible split");
  }), recursive=FALSE)
  loos = lapply(sp, function(i)
  { w = !(idC[,1] %in% i)
    pr = rowsum(t(dataForMC[,w]), k[w])/as.vector(table(k[w]))
    pr = exp(pr)-1;
    proto=t(pr); s = rowMedians(si[,!(colnames(si) %in% i)]); names(s) = rownames(proto);
    tmp = clustersInOtherDS(bu[,i], proto, s, maxIter=3)
    W = tmp$W/rowSums(tmp$W); rownames(W) = i;#rownames(cli);
    W[i,];
  })
  W = do.call(rbind, loos);
  pLoo = pFromW(W, pres[,unlist(sp)])
  save(tmp, W, pres, pLoo, sp, file=paste0(dataDir, "clustering/looMC/k", nk, ".RData"));
}

# Make prototypes of MC
#########################
tmp = assignClust(kms[["15"]], t(dataForMC), idC, annot, cli, minPatient=1, clLim=0)
o = c(1, 2, 5, 6, 3, 4, 7, 8, 11, 12, 9, 10, 13, 14); oo = o; oo[o]=1:14
tmp$k = oo[tmp$k]; tmp$tp = tmp$tp[o,]; tmp$Np = tmp$Np[o]; tmp$desc = tmp$desc[o];
K = tmp$k;

pr = rowsum(t(dataForMC), K)/as.vector(table(K))
pr = exp(pr)-1;
proto=t(pr);

sigma=rowMedians(si); names(sigma)=rownames(proto);
#save(proto, sigma, file="STstuff/data/MC.RData")

# Deconvolve MC at spot level
##############################
ers = lapply(rownames(cli), function(nm)
{ message(nm);

  x = readRDS(paste0(dataDir, "Robjects/counts/TNBC", nm, ".RDS"));
  cnts = t(x$cnts); spots = x$spots; rm(x);

  fm = computeMC(cnts, maxIter=2, rescale=FALSE, mccores=TRUE, quiet=FALSE)
  fm$spots = spots;

  saveRDS(fm, file=paste0(dataDir, "clustering/MC deconv/TNBC", nm, ".RDS"))
})

# Load deconvolution per spot
###############################
ms = lapply(names(ist), function(nm)
{ load(paste0(dataDir, "clustering/MC deconv/TNBC", nm, ".RData"));
  # Correct idSpot...
  w = which(!duplicated(idSpot$slide));
  sl = sub("([A-F][1-9])[0-9]+$","\\1", rownames(idSpot)[w]);
  rownames(m) = paste(sl[idSpot[,"slide"]], idSpot[,"spot"]);
  return(m);
}); names(ms) = names(ist)

m2s = lapply(ms, function(i) { i=i/rowSums(i); i[i<1e-2]=0; i=i/rowSums(i); i; });

fa = t(sapply(m2s, function(i) colMeans(i)));colnames(fa)=1:ncol(fa)
fa[fa<.01]=0

faN = fa; faN[faN<.01] = .01;

# Clustering in ecotypes
##########################
nm='ward.D2';

hc = hclust(dist(log10(100*(faN))), method=nm);
a = cutree(hc, 9); m = unique(a[hc$order]); n = rep(NA, length(m)); n[m]=seq_along(m); 
ecot = n[a];
#save(faN, ecot, file="STstuff/data/ET.RData")

# MC & ET in external datasets
###############################
# Install external datasets
ds = readRDS(paste0(dataDir, "external datasets/otherTNBC.RDS"))
for (i in names(ds))
{ x = ds[[i]];
  cl = x$cli; dn = x$dn;
  
  cl$bar = TNBCclassif(dn, version='bareche', shortName=TRUE)
  cl$leh = TNBCclassif(dn, version='lehmann', shortName=TRUE)
  cl$bur = TNBCclassif(dn, version='burstein', shortName=TRUE)
  cl$bur2 = TNBCclassif(dn, version='burstein2', shortName=TRUE)
  
  a = x$d;
  if (i=="METABRIC") { a=a/100; }
  a = a[rowSums(a)>100,];
  z = computeMC(a)

  H = z$W/rowSums(z$W); H[H<1e-3] = 0;
  colnames(H) = paste0("k", colnames(proto))
  x[["15"]]$H = H;
  
  et = computeET(H);
  x$ecot=et;
  
  x$cli=cl;
  ds[[i]] = x;
}
saveRDS(ds, file=paste0(dataDir, "external datasets/otherTNBC.RDS"))