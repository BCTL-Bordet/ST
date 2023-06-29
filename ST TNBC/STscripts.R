dataDir = "~/Data/Spatial/TNBC/datadist/"

#library(cancerData);
library(nonSTstuff)
library(STstuff)
#library(jpeg);
library(parallel);
options(mc.cores=4)#$detectCores()-1);
library(gplots);
library(matrixStats)
library(RColorBrewer)
require(corrplot)
library(NNLM)
library(EnsDb.Hsapiens.v86)

library(EBImage)
library(RANN);
library(mclust)
library(umap);
library(xCell)

library(survival)
library(xgboost)
#library(e1071)
#library(liquidSVM)
library(irlba);
library(pROC)
library(FNN)
library(grid)
library(reldist)
library(GSVA)
library(lme4);
  
midPoint = function(i, n)
{ s = seq(from=i[1], to=i[2], len=2*n+1);
  s[seq(from=2, to=length(s), by=2)];
}
  
colHeatmap = function(a, annots, lbl="", colAnnots, file=NULL, horizontal=TRUE)
{ #a = a[rowSds(a)>1e-5,]
  a = a-rowMins(a); a=a/rowMaxs(a);
  a = a[names(annots)[names(annots)%in%rownames(a)],];
  if (all(annots%in%names(colAnnots))) { z=colAnnots[annots]; names(z)=names(annots); annots=z; }
  if (!is.null(file)) 
  { if (horizontal) { pdf(file, width=nrow(a)*0.2+2, height=4+max(nchar(rownames(a)))/30); }
    else { pdf(file, height=nrow(a)*0.2+2, width=ncol(a)/4+max(nchar(rownames(a)))/20); }
  }
  if (!is.null(lbl)) { lbl = paste0(lbl, 1:ncol(a)); }
  if (horizontal)
  { par(mar=c(3,3,max(nchar(rownames(a))+4)/3,3)); b=a[,ncol(a):1];
    o=1:4; lbl=rev(lbl);
  }
  else
  { par(mar=c(3,max(nchar(rownames(a))+4)/2.5,.5,.5)); b=t(a); #a = a[nrow(a):1,]; 
    o = c(3, 4, 1, 2);
  }
  image(b, xaxt='n', yaxt='n', col=colorRampPalette(c('blue', 'yellow'))(100), xaxs='i', yaxs='i')
  pU = par('usr')[o];
  if (!is.null(lbl)) { mtext(lbl, side=1+horizontal, line=.2, at=midPoint(pU[3:4], ncol(a)), las=2) }
  if (horizontal)
  { text(midPoint(pU[1:2], nrow(a)), pU[4]+strheight("1"), rownames(a), xpd=NA, srt=45, adj=0,
      col=annots[rownames(a)])
  }
  else
  { mtext(rownames(a), at=midPoint(pU[1:2], nrow(a)), adj=1, side=2, line=.5,
      col=annots[rownames(a)], las=2)
  }
  
  if (horizontal) { legend(x=mean(pU[1:2]), y=pU[3]-5*strheight("1"), bty='n', xpd=NA, #text.width=strwidth(levels(co)),
      legend=names(colAnnots), fill=colAnnots, ncol=ceiling(length(colAnnots)/2), xjust=.5, yjust=0); }
  if (!is.null(file)) { dev.off(); }
}

factBareche = function(i) factor(i, levels=c("IM", "BL", "M", "MSL", "LAR"));

calcAllSig = function(x)
{ a = calcSig(x, sigs);
  b = t(gsva(x, sigH, mx.diff = TRUE, kcdf="Gaussian", parallel.sz=options("mc.cores")));
  cbind(a,b);
}

logAbs = function(x) log10(abs(x))*sign(x);

normIt = function(x)
{ x = x/rep(colSums(x), each=nrow(x)); log(1+x*1e4); }

classNewPrev = function(x)
{ x = x[,!(colnames(x) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))];
  rs = rowSums(x); x = x/rs; x[rs<1e3,] = NA;
  y = x[,colnames(x) != "Lymphocyte"]; y = y/rowSums(y);
  r = cbind(Tumor=rowSums(y[,c("Tumor", "in situ", "Tumor region"), drop=FALSE]), #
    Stroma=rowSums(y[,c("Low TIL stroma", "Vessels",
      "Stroma cell", "Nerve", "Acellular stroma", "Lymphoid nodule", "High TIL stroma")]))
}
classNew = function(x)
{ x = x[,!(colnames(x) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))];
  rs = rowSums(x); x = x/rs; x[rs<1e3,] = NA;
  y = x[,colnames(x) != "Lymphocyte"]; y = y/rowSums(y);
  r = cbind(Tumor=rowSums(y[,c("Tumor", "Tumor region"), drop=FALSE]), #
    Stroma=rowSums(x[,c("Low TIL stroma", "Vessels", "Lymphocyte",
      "Stroma cell", "Nerve", "Acellular stroma", "Lymphoid nodule", "High TIL stroma")]))
}

classOld = function(x)
{ x = x[,!(colnames(x) %in% c("None", "Artefact", "Outside", "Whitespace"))];
  rs = rowSums(x); x = x/rs; x[rs<1e3,] = NA;
  r = cbind(Tumor=rowSums(x[,c("Tumor", "Tumor highly infiltrated by TIL")]),
    Stroma=rowSums(x[,c("Low TIL stroma", "High TIL stroma", "Vessels", "Nodule lymphoid",
      "B cells", "Nerve")]));
}

evalStroma = function(x)
{ if (colnames(x)[2]=="Tumor") # New
  { x = x[,!(colnames(x) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))];
    rs = rowSums(x); x = x/rs; x[rs<1e3,] = NA;
    y = x[,colnames(x) != "Lymphocyte"]; y = y/rowSums(y);
    return(rowSums(y[,c("Low TIL stroma", "Vessels",
        "Stroma cell", "Nerve", "Lymphoid nodule", "High TIL stroma")]))
  }
  else
  { x = x[,!(colnames(x) %in% c("None", "Artefact", "Outside", "Whitespace"))];
    rs = rowSums(x); x = x/rs; x[rs<1e3,] = NA;
    return(rowSums(x[,c("Low TIL stroma", "High TIL stroma", "Vessels", "Nodule lymphoid",
        "B cells", "Nerve")]));
  }
}

mkmeans = function(x, ..., nstart=1000, minNstart=1000, minOK=5, fact=100, doPca=TRUE)
{ if (doPca && ncol(x)>nrow(x)) { xb=x; x = prcomp(x, scale=FALSE, center=FALSE)$x; } else { doPca=FALSE; }
  ns = minNstart;
  while(ns <= nstart)
  { res = mclapply(1:fact, function(i) kmeans(x, ..., nstart=ceiling(ns/fact)) );
    if (any(sapply(res, is, "try-error"))) { return(res); }
    w = sapply(res, function(i) i$tot.withinss);
    r = res[[which.min(w)]];
    r$Nok = sum(w==min(w)); r$nstart=ns; r$Nmax=fact;
    if (r$Nok >= minOK) { break; } else { ns = ns*10; }
  }
  
  if (doPca) { r$centers = rowsum(xb, r$cluster)/tabulate(r$cluster); } 
  return(r);
}

## MC and ET
#####################
cleanDeconv = function(i) { x = as.matrix(i[,-1]); rownames(x) = unlist(i[,1]); x; }

colsCateg = c(`Low TIL stroma`="orange",`High TIL stroma`="#ffff80",`Tumor`="#c80000",`Fat tissue`="#1a3399",
  `Artefact`="#b4b4b4", `Tumor highly infiltrated by TIL`="#338c07", `Canal galactophore`="#669999",
  `Necrosis`="#323232",`Vessels`="#ff6366", `Nodule lymphoid`="#b3b31a", `B cells`='#ff9999',
  `in situ`="#ade148", `Nerve`="#228caf", `None`="#ffffff", `Outside`="#b4b4b4", `Whitespace`="#b4b4b4") # For the colors

xCellPoss = c("aDC","Adipocytes","B-cells","Basophils","CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells",
  "CD4+ Tcm","CD4+ Tem","CD8+ naive T-cells","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","cDC","Class-switched memory B-cells",
  "CLP","CMP","DC","Endothelial cells","Eosinophils","Epithelial cells","Fibroblasts","GMP","iDC","Macrophages",
    "Macrophages M1","Macrophages M2","Mast cells","Memory B-cells","Monocytes","MPP","MSC","Myocytes",
    "naive B-cells","Neurons","Neutrophils","NK cells","NKT","pDC","Pericytes","Plasma cells","Preadipocytes",
    "pro B-cells","Smooth muscle","Tgd cells","Th1 cells","Th2 cells","Tregs")

colLeh = c(BL1="#1976D2", BL2="#4fc3f7", IM="#66BB6A", LAR="#EF5350", M="#FFEE58", MSL="#FFA726", UNS="#8D6E63")
colBar = c(BL="#1976D2", IM="#66BB6A", LAR="#EF5350", M="#FFEE58", MSL="#FFA726")
colBu = colLeh[c("IM", "BL1", "MSL", "LAR")];
names(colBu) = c("BLIA","BLIS", "MES", "LAR");
colBar = colBar[c("IM", "BL", "M", "MSL", "LAR")];

ring.genes = c("FTH1", "EEF2", "BEST1", "LRRC59", "PRDX1", "CD63", "DYNC1H1", "ENO1",
                "PSMB3", "RNF187", "RNASE1", "CFL1", "GRN", "UBC", "TAX1BP3", "COX4I1",
                "CUTA", "NME1", "H3F3B", "AKR7A2", "IMPDH2")

colEcot = palette.colors(palette="Polychrome 36")[-c(1, 5)][1:12]

bubullePlot = function(x, cmps, what="Other", fileName, width=NULL,
  leftMar=NULL, horizontal=TRUE)
{ co=NULL;
  if (what=="Genes")
  { x = x[, colnames(x) %in% rownames(geneList)];
    x = x[,order(-match(colnames(x), rownames(geneList)))];
    #if (v2 == "Stroma")
    #{ x = x[,!(geneList[colnames(x), "Role.of.genes"] %in% c("Methylation", "HRD", "Cycling", "MMR"))] }
    co = factor(geneList[colnames(x), "Role.of.genes"], levels=unique(geneList[, "Role.of.genes"]));
    co = palette.colors(palette="Polychrome 36")[-2][unclass(co)];
    names(co) = colnames(x) = geneList[colnames(x), "List.of.interesting.genes"]
  }
  if (what=="xCell")
  { x = x[,!grepl("Score", colnames(x)),drop=FALSE];
    x = x[,order(-match(colnames(x), names(colXct))),drop=FALSE];
    co = colXct[colnames(x)]; names(co) = colnames(x);
  }
  if (what=="Sigs")
  { x = x[,intersect(names(sigInfo), colnames(x)),drop=FALSE];
    co = colSig[sigInfo[colnames(x)]]; names(co) = colnames(x);
  }
  ok = which(!is.na(cmps));
  ps = do.call(cbind, tapply(seq_along(cmps), cmps, function(w)
  { apply(x,2,function(i) wilcox.test(i[w], i[setdiff(ok, w)])$p.value)
  }, simplify=FALSE))
  p = ps; p[] = p.adjust(ps, method='fdr');
  w = rownames(p)[which(rowAnys(p<=.05))];

  if (length(w)==0) { return(invisible("Nothing significant")); }
  
  p2 = p.adjust(ps[w,], method='fdr');
  cutOff = max(p2[p[w,]<=.05])+1e-10;

  if (is.null(width))
  { width = max(nchar(w))+3 + ncol(ps)*max(nchar(colnames(ps))+3.5) + 30; width=width*.035;
    leftMar = 1+max(nchar(w))/2; 
  }
  
  if (what !="Other") { fileName = sub("pdf$", paste0(what, ".pdf"), fileName); } 
  
  cairo_pdf(fileName, height=(length(w)+7)/7, width=width)
  par(mar=c(3,leftMar,.5,.5), mgp=c(1.5,.5,0), cex=.5);
  dotPlot(x[,w,drop=FALSE], factor(cmps), col.lbl=co[w], maxP=cutOff);
  dev.off()
  
  if (horizontal)
  { sl = pmax(0.06*(max(nchar(w))+3)*.7, 0.06*max(nchar(colnames(ps))+8)); sb=0.06*leftMar*.7+2;  
    cairo_pdf(sub(".pdf$", ".horiz.pdf", fileName), width=2*(length(w)+7)/7+sl,
      height=ncol(ps)*.45+sb)
    #par(mar=c(leftMar*.7+7, leftMar*.7,.5,.5), mgp=c(1.5,.5,0));
    par(mai=c(sb, sl, .01, .01), mgp=c(1.5,.5,0));
    dotPlot(x[,w,drop=FALSE], factor(cmps), col.lbl=co[w], horizontal=TRUE, maxP=cutOff, inMa=.8, oma=NULL);
    legendDotPlot(par('usr')[1], par('usr')[3]-(leftMar+4)*.75*strheight("M"),
      horizontal=TRUE)
    dev.off()
  }
  
  return(invisible(list(p, w)));
}

loadData = function(bd, loadImage=TRUE, removeFold=TRUE, cutOff=500, geneMap=NULL, autoCut)
{nm = sub(".+/(CN[0-9]+)/([CDE][12]).*", "\\1_\\2", bd);
  if (loadImage) { im = readImage(paste0(bd, "/small.jpg")); } else { im=NULL; }
  load(paste0(bd, "/selection.RData"));
  if (!is.null(autoCut[[nm]])) { spots = spots[-autoCut[[nm]],]; cnts = cnts[-autoCut[[nm]],]; }
  w = rowSums(cnts)>cutOff; cnts = cnts[w,]; spots=spots[w,];
  g = sub("\\.[0-9]+$", "", colnames(cnts));
  if (is.null(geneMap))
  { gid = ensembldb::select(EnsDb.Hsapiens.v86, keys=g, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    colnames(cnts) = gid[match(gid[,2], g), 1];
  }
  else { colnames(cnts)  = geneMap[g]; }
  cnts = t(cnts);
  cnts = rowsum(cnts, rownames(cnts));
  spots[,c("pixel_x", "pixel_y")] = spots[,c("pixel_x", "pixel_y")]/8;
  if (file.exists(f<-paste0(bd, "/plis.RDS")))
  { artefacts=readRDS(f)[colnames(cnts),];
    if (removeFold) { w=artefacts[,"OK"]>.99*rowSums(artefacts); cnts = cnts[,w]; spots=spots[w,]; artefacts=NULL; }
  } else { artefacts=NULL; }
  if (file.exists(f<-paste0(bd, "/annotBySpot.RDS"))) { annot=readRDS(f)[colnames(cnts),]; } else { annot=NULL; }
  if (file.exists(f<-paste0(bd, "/annotBySpotNew.RDS")))
  { annotNew=readRDS(f)[colnames(cnts),];
    if (removeFold)
    { w = annotNew[,"Artefacts"]/rowSums(annotNew) <= .01;
      cnts = cnts[,w]; spots=spots[w,]; annotNew=annotNew[w,];
    }
  } else { annotNew=NULL; }
  if (loadImage && ncol(cnts)>0) { imSpot = spotXY(im, spots, diam=6) } else { imSpot=NULL; }
  return(list(im=im, cnts=cnts, imSpot=imSpot, annot=annot, spots=spots, artefacts=artefacts, annotNew=annotNew))
}

# Utilities and stats
#######################
getCoxP = function(x, HR=FALSE)
{ p=summary(x)$logtest[3]
  if (HR) { return(c(p, coef(x))); } else { return(p); }
}
 
moulinette = function(x, y, w=NULL)
{ if (!is.null(w))
  { if (is.factor(w) || is.character(w))
    { w = factor(w);
      return(do.call(cbind, tapply(1:nrow(x), w, function(i) moulinette(x,y,i))))
    }
    x = x[w,];
    y = y[w];
  }
  if (is.data.frame(x))
  { x = x[,sapply(x, is.numeric)];
  }
  
  if (is(y, "Surv"))
  { p = apply(x, 2, function(i) summary(coxph(y~i))$logtest["pvalue"]);
    return(p);
  }
  
  if (length(unique(y)) > 2)
  { p = apply(x, 2, function(i) suppressWarnings(kruskal.test(i~y))$p.value);
    return(p);
  }
  y = unclass(factor(y));
  #p = pt.test(t(x[y==1,]), t(x[y==2,]));
  p = apply(x, 2, function(i) suppressWarnings(wilcox.test(i[y==1], i[y==2]))$p.value);
  return(p);
}

fullMat = function(x, g, fill=0)
{ g2 = intersect(g, rownames(x))
  m = matrix(fill, ncol=ncol(x), nrow=length(g), dimnames=list(g, colnames(x)))
  m[g2,] = x[g2,];
  return(m)
}

fullMatSparse = function(x, g)
{ g2 = setdiff(g, rownames(x))
  m = sparseMatrix(c(), c(), dims=c(length(g2), ncol(x)), dimnames=list(g2, colnames(x)))
  m = rbind(x, m);
  m = m[g,]
  return(m)
}

fullMat2 = function(x, g)
{ do.call(cbind, lapply(x, function(i)
  { r = rep(0, length(g)); names(r) = g; r[names(i)] = i; r;
  }))
}

fullMat3 = function(x, g)
{ do.call(cbind, lapply(x, function(i)
  { r = matrix(0, nrow=length(g), ncol=ncol(i)); rownames(r) = g; colnames(r) = colnames(i);
    i = i[intersect(g, rownames(i)),];
    r[rownames(i),] = i; r;
  }))
}

norm2 = function (d, lim = 20) 
{  if (any(d < 0, na.rm = TRUE)) {
        stop("Negative value(s) present, probably already in log")
    }
    w = rowSums(d, na.rm = TRUE) > lim
    d = 1e4*d/rep(colSums(d, na.rm = TRUE), each = nrow(d))
    d = d[w,];
 
    d = log10(d + 1)
    return(d)
}

rescaleForClassif = function(fm, cnts)
{ cnts = cnts[intersect(rownames(cnts), fm$g), ]
  x = 1e4*cnts/rep(colSums(cnts), each=nrow(cnts));
  x2 = matrix(0, nrow=length(fm$center), ncol=ncol(x)); colnames(x2) = colnames(x); rownames(x2) = names(fm$center);
  i = intersect(rownames(x), names(fm$center)); x2[i,] = x[i,]; x = x2; rm(x2);
  #x = x[names(fm$center), ] # To check
  x = log(x+1);
  x = (x-fm$center)/fm$scale;
 
  xx = t(qr.solve(fm$rot, x));
}

classifTumor = function(fm, sv, cnts)#, what=c("Tumor", "Stroma"))
{ xx = rescaleForClassif(fm, cnts)
  pr = predict(sv, newdata=xx);
  #pr < c(Tumor=.7, Stroma=.5)[what];
}

makePR = function(x, Ng=1000)
# Genes in columns
{ nt = colMeans(x>=5)
  x = 1e4*x/rowSums(x);
  x = x[, colMeans(x!=0)>.05 & nt>.01]
  x = log(x+1)
  if (ncol(x)>Ng)
  { z = colSds(x)/colMeans(x);
    w = rank(-z)<=Ng;
    x = x[,w];
  }
  return(x);
}

maskAndResize = function(img)
{ x = gblur(img, 10)
  x2 = resize(x, h=1000)

  y = channel(x2, "grey")

  b = y < otsu(y);
  b2 = closing(b, makeBrush(7, 'disc'))
  b2 = opening(b2, makeBrush(15, 'disc'))

  lbl = bwlabel(b2);
  t = table(lbl);
  if (any(t)<1000) { lbl = rmObjects(lbl, which(t<1000)-1); }
  b2 = lbl>0;
  return(list(mask=b2, img=x2));
}

centerAndRot = function(mr, f)
{ if (is.list(mr)) { p = centerAndRotParams(mr[[1]], f); }
  else { p = centerAndRotParams(mr, f); }
  transfo(mr, p);
}

centerAndRotParams = function(mr, f)
{ return(c(dim(mr)[1:2]/2-f[1:2], theta=-mean(f[grep("theta", names(f))])*180/base::pi));
}

qual = function(x, y, theta, delta)
{ if (theta!=0) { y2 = rotate(y, theta, output.dim=dim(y)[1:2], output.origin=dim(y)[1:2]/2, filter='none'); }
  else { y2 = y; }
  y2 = translate(y2, delta);
  w = x[,,2]==1 & y2[,,2]==1;
  mean(x[,,1][w] == y2[,,1][w])
}

qualRI = function(x,y,theta,delta)
{ if (theta!=0) { y2 = rotate(y, theta, output.dim=dim(y)[1:2], output.origin=dim(y)[1:2]/2, filter='none'); }
  else { y2 = y; }
  y2 = translate(y2, delta);
  #w = x!=0 & y2!=0;
  adjustedRandIndex(x, y2)
}

makeImage = function(x, lbl, rescale=TRUE)
{ if (rescale) { x = x-min(x); x = x/max(x); }
  r = lbl;
  x = c(0, x);
  r[] = x[lbl+1];
  return(r);
}

pdist = function(x, y)
{ matrix(nrow=nrow(y), colSums((t(x)[,rep(1:nrow(x), each=nrow(y))]-t(y)[,rep(1:nrow(y),nrow(x))])^2))
}


## Clustering functions
#########################
assignClust = function(km, y, idC, annot, cli, clLim=1, minPatient=3, rmProbOnly=FALSE)
{ library(rdist);
  k = km$cluster;
  
  # Reassign cluster with only one/two patient
  t = tapply(idC[,"id"], k, function(i) length(unique(i)))
  prob = which(t<=minPatient)
  if (any(prob))
  { d = cdist(y,km$centers[-prob,,drop=FALSE])
    k = apply(d, 1, which.min)
  }
  
  # Reassign those that go into only one cluster
  w = names(which(colSums(table(k, idC[,1])>0)<=clLim));
  for (i in w)
  { ww = which(idC[,1]==i)
    pr = rowsum(y[-ww,], k[-ww])/as.vector(table(k[-ww]))
    d = pdist(pr, y[ww,]);
    kn = apply(d,1,which.min)
    if (any(kn!=unclass(k)[ww])) { message("Youpi"); k[ww] = kn; } #levels(k)[kn] }
  }
  
  if (rmProbOnly) { return(k); }

  t = table(k, cli[idC[,1], "barPB"])
  t2 = rowsum(annot, k)/tabulate(k);
  o = order(t2[,"Fat tissue"]>.1, -t2[,"Tumor"])
  k = match(k, o); t2 = t2[o,]; rownames(t2) = NULL;
  
  Np = tapply(idC[,1], k, function(i) length(unique(i)))
  
  bb = factor(cli[idC[,'id'], "barPB"]);
  tp = do.call(rbind, tapply(seq_along(k), k, function(i)
  { tapply(i, bb[i], function(j) length(unique(idC[j,"id"])));
  }));
  
  desc = apply(t2, 1, function(i) 
  { if (i["Lymphocyte"]>.1) { b=paste0(" Lym(", round(i["Lymphocyte"]*100),"%)")} else { b=""; }
    if (i[1]>.1) { return(paste0("Fat(",round(i[1]*100), "%)", b)); }
    if (i["Tumor"]>.25) { return(paste0("Tum(", round(i["Tumor"]*100), "%)",b)) }
    return(paste0("Stro(", round(i["Stroma"]*100), "%)",b))
  } )

  return(list(k=k, tp=tp, Np=Np, desc=desc));
}

calcPres = function(k, idC)
{ pres = do.call(rbind, tapply(k, idC[,"id"], function(i) tabulate(i, nbins=max(k))))#>0;
  colnames(pres) = paste0("k", 1:ncol(pres))
  pres = pres[rownames(cli),]>0;
  pres;
}

## Display
############
ff = function(x, q) { x = (x-q[1])/(q[2]-q[1]); x[x<0]=0; x[x>1]=1; x; }

spotXY = function(im, spots, diam=50)
{ xy = sapply(spots[,c("pixel_x","pixel_y")], as.integer); di2 = floor(diam);
  base = cbind(-di2:di2, rep(-di2:di2, each=2*di2+1));
  base = base[base[,1]^2+base[,2]^2 <= diam^2,]
  xys = base[rep(1:nrow(base), nrow(xy)),] + xy[rep(1:nrow(xy), each=nrow(base)),]
  ids = rep(1:nrow(xy), each=nrow(base))
  xs = xys[,1]+xys[,2]*nrow(im);
  w = xs>0 & xs<=prod(dim(im)[1:2]); xs = xs[w]; ids = ids[w];
  return(list(xy=xs, id=ids))
}

makeImSpot = function(im, spots)
{ xy = spots[,5:6];
  imSpot = Image(0, dim(imB)[1:2])
  for (i in 1:nrow(xy))
  { imSpot = .Call(EBImage:::C_drawCircle, imSpot, c(as.integer(xy[i,]), 0L, 50L), c(i,0,0), fill=TRUE) }
  spl = split(seq_along(imSpot), as.vector(imSpot))[-1];
  xy = unlist(spl); id = as.integer(rep(names(spl), sapply(spl, length))); 
  names(xy) = NULL;
  return(list(imSpot=imSpot, xy=xy, id=id))
}

categorizeSpots = function(imR, imSpot, cls)
{ spotCateg = do.call(rbind, tapply(imR[imSpot$xy], imSpot$id, function(i) tabulate(i, length(cls)) )); 
  colnames(spotCateg) = cls;
  return(spotCateg);
}

plt = function(cols="palette")
{ if (cols=="palette") { return(c(palette(), 'white', palette()[-1])); }
 #require(RColorBrewer)
 #if (cols %in% rownames(brewer.pal.info)) { return(c('black', brewer.pal(12, cols))) }
 #c('black', brewer.pal(12, "Set3"))
 library(Polychrome)
 c("black", glasbey.colors()[-5])
}

colLabls = function(x, cols=plt())
{ cols = t(col2rgb(cols)/255);
  cols = rbind(cols, 1, cols[-1,]);
  
  y = cols[as.vector(x)+1,]
  y = array(y, dim=c(dim(x), 3))
  im = Image(y, colormode="Color")
  return(im);
}

makeCol = function(n, rescale=TRUE, asRGB=FALSE)
{ #if (any(is.na(n)))
  #{ m = rep(NA, length(n));
  #  if (is.vector(n) || is.factor(n) || any(dim(n)==1)) { w = !is.na(n); n=n[w]; }
  #  else { w = rowAlls(!is.na(n)); n = n[w,]; }
  #  m[w] = 1:sum(w);
  #  ww = imSpot$id %in% which(w); imSpot$id=m[imSpot$id[ww]]; imSpot$xy=imSpot$xy[ww];
  #}
  if (is.logical(n)) { n = factor(n); }
  if (is.factor(n))
  { if (!is.null(cols) && all(levels(n) %in% names(cols))) { cols = cols[levels(n)]; }
    n = unclass(n);
  }
  if (is.integer(n))
  { if (min(n)<1) { n = n+1-min(n); }
    if (is.null(cols))
    { if (max(n, na.rm=TRUE)>12) { stop("Max 12 colors ATM"); }
      cols = palette();
    }
    cols = col2rgb(cols)/255;
    n1 = cols[,n,drop=FALSE];
    n1R = n1[1,]; n1G = n1[2,]; n1B = n1[3,];
  }
  else
  { if (is.vector(n) || any(dim(n)==1))
    { if (rescale) { q = quantile(n, c(.05, .95), na.rm=TRUE); n = ff(n, q); } else { q = 0:1; }
      n1R = 1-n*2; n1R[n1R<0] = 0;
      n1G = n*2-1; n1G[n1G<0] = 0;
      n1B = rep(0, length(n1G))
    }
    else
    { if (rescale) { n = apply(n,2, function(i) { q = quantile(i, c(.05, .95), na.rm=TRUE); ff(i,q); }) } else { q=0:1; }
      n1R = n[,1]; n1G = n[,2];
      if (ncol(n)==3) { n1B = n[,3]; } else { n1B = rep(0, length(n1G)); }
    }
  }
  i = seq(from=0,to=1, len=100);
  if (asRGB) { rgb=cbind(R=n1R, G=n1G, B=n1B); } else
  { w = !is.na(n1R); rgb=rep(NA, length(n1R)); rgb[w]=rgb(n1R[w], n1G[w], n1B[w]);
  }
  return(list(rgb=rgb, scale=list(lims=q, valCols=i, cols=rgb(pmax(0, 1-2*i), pmax(2*i-1,0), 0))));
}

colSpots = function(n, im, imSpot, rescale=TRUE, cols=NULL, retScale=FALSE)
{ if (is.list(im)) { imSpot = im$imSpot; im = im$im; }
  if (length(dim(im))==2) { stop("Only color images (not greyscale)")}
  
  co = makeCol(n, rescale=rescale, asRGB=FALSE); n1 = co$rgb;
  delt = dim(im)[1]*dim(im)[2]
  im@.Data[imSpot$xy] = n1[imSpot$id, "R"];
  im@.Data[imSpot$xy+delt] = n1[imSpot$id, "G"];
  im@.Data[imSpot$xy+2*delt] = n1[imSpot$id, "B"];
  if (retScale)
  { return(list(im=im, scale=co$scale))#lims=q, valCols=i, cols=rgb(pmax(0, 1-2*i), pmax(2*i-1,0), 0))) 
  }
  return(im);
}

colSpotsOff = function(nm, n, ...)
{ im = readImage(paste0("~/Data/Spatial/TNBC/spots/", sub("_", "/", nm), "/small.jpg"))
  spots = readRDS(paste0("~/Data/Spatial/TNBC/spots/", sub("_", "/", nm), "/selectionSpots.RDS"))
  spots[,c("pixel_x", "pixel_y")] = spots[,c("pixel_x", "pixel_y")]/8;
  if (is.vector(n)) { spots = spots[names(n),]; }
  if (is.matrix(n)) { spots = spots[rownames(n), ]}
  imSpot = spotXY(im, spots, diam=6)
  colSpots(n=n, im=im, imSpot=imSpot, ...)
}

dispAll = function(dta, toDisp=NULL, f=NULL, rescale=TRUE)
{ par(mfrow=c(1,length(dta)))
  for (i in seq_along(dta))
  { if (is.null(toDisp)) { x = dta[[i]]$cnts; } else { x = toDisp[[i]]; }
    if (!is.null(f)) { if (any(colnames(x)%in%f)) { x = x[,f] } else { x=x[f,]; if (is.matrix(x)) { x=t(x) } } }
    display(colSpots(x, dta[[i]]$im, dta[[i]]$imSpot, rescale=rescale), 'r')
  }
}

dispAllOld = function(imAll, x, nr=NULL, method='r', ...)
{ if (is.character(x)) # single gene then
  { x = lapply(imAll$ge, function(i) as.numeric(i[x,]));
  }
  r = lapply(seq_along(x), function(i) colSpots(x[[i]], imAll$im[[i]], imAll$spotPos[[i]]));
  di = sapply(r, dim)[1:2,];
  if (any(apply(di, 1, function(i) length(unique(i))>1))) { message("Not all same dim"); return(r); }
  r = do.call(combine, r);
  if (method == "n") { return(r); }
  display(r, method=method, ...)
  return(invisible(r));
}

dispX = function(n, prec)
{ n=floor(nrow(cols)*(n/(1+1e-10)))+1;
  x = as.array(prec$x);
  x[,,1][prec$w] = x[,,1][prec$w]*.5 + .5*cols[n[prec$idS],1]
  x[,,2][prec$w] = x[,,2][prec$w]*.5 + .5*cols[n[prec$idS],2]
  #browser();
  x = Image(x, dim=dim(prec$x), colormode="Color")
  #display(x, method='browser') # Somehow does not work
  return(invisible(x))
}

dispX2 = function(i, f, css, imSpots, method='raster', resize=NULL)
{ j = 2*i - 1;
  if (any(colnames(css[[j]])==f)) { n1 = css[[j]][,f]; n2 = css[[j+1]][,f]; }
  else { n1 = css[[j]][f,]; n2 = css[[j+1]][f,]; }
  nn = c(n1,n2); nn = nn[nn > -2.5];
  q = quantile(nn, c(.05, .95));
  
  n1 = ff(n1, q); n2 = ff(n2, q);
  x1 = dispX(n1, imSpots[[j]]);
  x2 = dispX(n2, imSpots[[j+1]]);
  x = combine(x1, x2)
  if (!is.null(resize)) { x = resize(x, h=nrow(x)*resize); }
  display(x, method=method, all=TRUE);
  return(invisible(q))
}

dispX2b = function(id, f, css, imSpots, method='raster', resize=NULL)
{ if (id %% 2 == 0) { j = id-1; } else { j = id+1; }

  if (any(colnames(css[[id]])==f)) { n1 = css[[id]][,f]; n2 = css[[j]][,f]; }
  else { n1 = css[[id]][f,]; n2 = css[[j]][f,]; }
  nn = c(n1,n2); nn = nn[nn > -2.5];
  q = quantile(nn, c(.05, .95));
  ff = function(x, q) { x = (x-q[1])/(q[2]-q[1]); x[x<0]=0; x[x>1]=1; x; }
  n1 = ff(n1, q); n2 = ff(n2, q);
  x = dispX(c(n1,n2), imSpots[[id]]);
  if (!is.null(resize)) { x = resize(x, h=nrow(x)*resize); }
  display(x, method=method, all=TRUE);
  return(invisible(x))
}

dispCL = function(im1, im2, ..., resize=NULL)
{ set.seed(10); j = colorLabels(im1);
  set.seed(10); j2 = colorLabels(im2);
  x = combine(j, j2);
  if (!is.null(resize)) { x = resize(x, h=nrow(x)*resize); }
  display(x, all=TRUE, ...)
}

drawScale = function(a, cols, log=FALSE, ylab)
{ par(mar=c(.5,3,.5,.2), mgp=c(1.5,.5,0))
  if (log) { log='y'; } else { log=''; }
  plot(1, type='n', xlim=c(0,1), ylim=a, log=log, xaxt='n', xlab='', ylab=ylab, xaxs='i', yaxs='i')
  if (log=='y') { s = exp(seq(log(a[1]), log(a[2]), len=nrow(cols)+1)) }
  else { s = seq(a[1], a[2], len=nrow(cols)+1); }
  coo = rgb(cols[,1], cols[,2], 0);
  for (i in 1:nrow(cols))
  { rect(0, s[i], 1, s[i+1], lty=0, col=coo[i])
  }
}

plotPie = function(x, y, radius, values, colors, Nsegments=100)
{ xx = cbind(radius*sin(2*base::pi*(0:Nsegments)/Nsegments), 
    radius*cos(2*base::pi*(0:Nsegments)/Nsegments));
  
  pos = lapply(seq_along(x), function(ii)
  { va=values[ii,];
    v = c(0, round(cumsum(Nsegments*va/sum(va)))); 
    w = which(diff(v)>0); v = c(0, v[w+1])+1; co=colors[w];
  
    xif = matrix(NA, nrow=nrow(xx) + 3*length(w), ncol=2); 
    for (i in seq_along(w))
    { xif[(v[i]:v[i+1]) + 3*(i-1), ] = xx[v[i]:v[i+1], ]
      xif[v[i+1] + 3*(i-1)+1, ] = 0;
    }
    xif = xif+rep(c(x[ii], y[ii]), each=nrow(xif));
    return(list(xif, co));
  })
  
  x2 = do.call(rbind, lapply(pos, function(i) i[[1]]));
  co = unlist(lapply(pos, function(i) i[[2]]));
  
  #browser();
 
  polygon(x2[,1], x2[,2], col=co, border=NA);
}

pointsDta = function(x, ...)
{ fa = 1/(max(dim(x$im)[1:2])/238);
  #message(fa);
  points(x$spots[,"pixel_x"], x$spots[,"pixel_y"], ..., cex=fa, pch=16)
}

#########################
## Shiny...
registration = function(im1, im2)
{ shinyApp(
    ui = fluidPage(
      verticalLayout(
        flowLayout(numericInput("theta", "Theta", 0, -180, 180),
          numericInput("dx", "Translation x", 0),
          numericInput("dy", "Translation y", 0),
          sliderInput("f1", "Opacity fixed", min=0, max=100, value=50, step=5),
          sliderInput("f2", "Opacity moving", min=0, max=100, value=50, step=5),
          checkboxInput("inv", "Flip!")),
          plotOutput("im", height="800px")#,
          #sliderInput("zoom", "Zoom", min=10, max=400, value=100, step=10)
    )), 
    server = function(input, output) {
      dispIt = function ()
      { par(mfrow=c(1,3));
        if (input$inv) { x = flip(im2); } else { x = im2; }
        if (input$theta != 0 ) { x = rotate(x, input$theta, filter='none', output.dim=dim(im2)[1:2], bg.col='white') }
        x = translate(x, c(input$dx, -input$dy), bg.col='white')
        display(im1, 'r');
        i = round(seq(0, nrow(im1), len=10));
        segments(0, i, nrow(im1), i); segments(i, 0, i, nrow(im1))
        display(x, 'r');
        segments(0, i, nrow(im1), i); segments(i, 0, i, nrow(im1))
        display((im1*input$f1+x*input$f2)/100, 'r')
      }
      output$im = renderPlot(dispIt());#
    }
  )
}

##########################
## RNA seq stuff
normRNAseq = function(d, lim=20, small=1e-3)
{ if (any(d<0, na.rm=TRUE)) { stop("Negative value(s) present, probably already in log"); }
  d = d[rowSums(d, na.rm=TRUE) > lim,];
  d = d/rep(colMeans(d, na.rm=TRUE), each=nrow(d));
  d = log10(d + small);
  return(d);
}


#####################################
## Meta...
findClusts = function(x, W, Niter=10)
{ x = x[rowSums(x)>100,];
  g = intersect(rownames(x), rownames(W)); x=x[g,]; W=W[g,];
  a = x/rep(colSums(x), each=nrow(x));
  b = rowMeans(a) / rowMeans(W);
  bs = list(); Hs=list();
  for (iter in 1:Niter)
  { bs[[iter]] = b;
    message("Iter:", iter);
    W2 = W*b;
    #H = do.call(cbind, mclapply(1:ncol(a), function(i)
    #  nnlm(W2,x[,i,drop=FALSE],check.x=FALSE, loss='mkl')$coefficients))
    H = nnlm(W2, x, n.threads=options("mc.cores"), loss='mkl',check.x=FALSE)$coefficients
    Hs[[iter]] = H;
    rec = W2 %*% H;
    #o = 1/unlist(mclapply(1:nrow(rec), function(i)
    #  optimize(toOpt, c(.01, 10), X=a[i,], Y=rec[i,])$minimum))
  
    b = b * rowSums(x) / rowSums(rec)
  }
  bs[[Niter+1]] = b;
  H = H/rep(colSums(H), each=nrow(H));
  return(list(H=H, b=b, bs=bs, Hs=Hs));
}

cPar = function(p)
{ if (length(p)==0) { return(c(0, 0, rep(NA, 12))); }
  Np = sum(p);
  p = p/sum(p);
  o = sort(p, decreasing = TRUE); names(o) = NULL;
  co = cumsum(o);
  itv = findInterval(pp<-c(.1, .25, .5, .75, .9, .95), co)+1;
  names(itv) = paste0("N", pp);
  c(N=length(p), Np=Np, evenness = -sum(p * log(p), na.rm = TRUE)/log(sum(!is.na(p))),
    entropy = -sum(p * log(p), na.rm = TRUE),
    gini = gini(p[!is.na(p)]),
    `gini-simpson` = 1 - sum(p^2, na.rm = TRUE),
    top=o[1], `2nd`=o[2], itv);
}

calcP = function(cc, cmps)
{ ok = which(!is.na(cmps));
  ps = do.call(cbind, tapply(seq_along(cmps), cmps, function(w)
  { apply(cc,2,function(i) wilcox.test(i[w], i[setdiff(ok, w)])$p.value)
  }))
  p = ps; p[] = p.adjust(ps, method='fdr');
  return(list(p=ps, fdr=p));
}

########################################
## Visium
cleanDta = function(dta)
{ dta$spots[,"pixel_x"] = dta$spots[,"imagerow"] * dta$scale$lowres;
  dta$spots[,"pixel_y"] = dta$spots[,"imagecol"] * dta$scale$lowres;
  dta$im = Image(dta$im, colormode="Color")
  if (ncol(dta$im)<600) { dta$im = abind(dta$im, Image(0, colormode='Color', dim=c(600, 600-ncol(dta$im), 3)), along=2); }
  if (nrow(dta$im)<600) { dta$im = abind(dta$im, Image(0, colormode='Color', dim=c(600-nrow(dta$im), 600, 3)), along=1); }
  dta$cnts = as.matrix(dta$cnts); dta$cnts = dta$cnts[rowSums(dta$cnts)>0,]
  dta$imSpot = spotXY(dta$im, dta$spots, diam=2);
  dta;
}

#####################################
## heatmap.3 - like heatmap, but allow for multiple annotations
##############################
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){

    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || anyNA(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || anyNA(Colv))
        Colv <- FALSE
    else if (isTRUE(all.equal(Colv,"Rowv")) && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors) && !is.null(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
        }

        if (!missing(RowSideColors)&& !is.null(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)&& !is.null(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(colnames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!missing(ColSideColors) && !is.null(ColSideColors)) {

        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }

        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}