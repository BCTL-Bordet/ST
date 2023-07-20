source("~/prog/ngs/R scripts/STscripts.R")

projDir = "~/Documents/BCTL/ST/TNBC/projections/";

cli = readRDS(paste0(dataDir, "Clinical/Clinical.RDS"))

# Project regressed annotations
##################################
ist = readRDS(paste0(dataDir, "classification/classifAll.RDS"))
ist = lapply(ist, function(ii)
{ ii = ii$pr;
  ii = ii[,c('Fat tissue','in situ','Lactiferous duct','Lymphoid nodule','Necrosis',
    "Tumor", "Stroma", "Lymphocyte", 'Vessels')];
})
ist = ist[rownames(cli)] # Classification per spot

co = colAnn2;
names(co)[names(co)=="Stroma cell"] = "Stroma";
co = co[colnames(ist[[1]])];

dir.create(paste0(projDir, "annotRegressed"));
dir.create(paste0(projDir, "tumorStroma"));
ers = mclapply(rownames(cli), function(nm)
{ img = readRDS(paste0(dataDir, "Robjects/imagesSmall/TNBC", nm, ".RDS"));
  i2 = ist[[nm]]; i2 = i2/rowSums(i2);
  id = rep(seq_along(img), sapply(img, function(i) nrow(i$spots)));
  if (length(id) != nrow(i2)) { stop("Not same number of spots for ", nm); }
  pdf(paste0(projDir, "annotRegressed/TNBC", nm, ".pdf"), height=3.5, width=3.5*length(img)+1.5)
  par(mfrow=c(1, length(img)), cex=.66, mex=.66, omi=c(0,0,0,1.5), mar=c(0,0,0,0))
  for (i in seq_along(img))
  { display(img[[i]]$im, 'r')
    fa = 1/(max(dim(img[[i]]$im)[1:2])/238);
    plotPie(img[[i]]$spots[,"pixel_x"], img[[i]]$spots[,"pixel_y"], values=i2[id==i,], colors=co, radius=fa*1.6, Nsegments=100)
  }
  w = colMaxs(i2)>.05;
  legend(par('usr')[2]+diff(par('usr')[1:2])*.01, par('usr')[4]+10, colnames(i2)[w], fill=co[w], xpd=NA)
  dev.off();
  
  i2 = cbind(Tumor=i2[,"Tumor"], Stroma=rowSums(i2[,!(colnames(i2) %in% c("Tumor", "Necrosis","in situ"))]),
      Other=rowSums(i2[,c("Necrosis","in situ")]));
  if (length(id) != nrow(i2)) { stop("Not same number of spots for ", nm); }
  pdf(paste0(projDir, "tumorStroma/TNBC", nm, ".pdf"), height=3.5, width=3.5*length(img)+1)
  par(mfrow=c(1, length(img)), cex=.66, omi=c(0,0,0,1))
  for (i in seq_along(img))
  { display(img[[i]]$im, 'r')
    fa = 1/(max(dim(img[[i]]$im)[1:2])/238);
    plotPie(img[[i]]$spots[,"pixel_x"], img[[i]]$spots[,"pixel_y"], values=i2[id==i,], colors=c('red', 'green', 'black'),
      radius=fa*1.6)
  }
  legend('topright', c("Tumor", "Stroma", "Other"), fill=c('red', 'green', 'black'), inset=c(-.25, 0), xpd=NA, cex=.8)
  dev.off();
}, mc.cores=4)

# Signatures and so on by spots
################################
plot3 = function(dta, x, col, id)
{ cols = color(x, col=col);
  par(mfrow=c(1, length(dta)), mar=c(0,0,0,0), cex=.66)
  for (j in seq_along(dta))
  { display(dta[[j]]$im, 'r')
    pointsDta(dta[[j]], col=cols[id==j])
  }
}
color = function(x, minx=min(x), maxx=max(x), col=colorRampPalette(c('blue', 'yellow'))(100))
{ x = x-minx; x = x/maxx; x[x>1]=1; x[x<0]=0;
  x = 1+99*x;
  col[round(x)];
}

# Get the mean and standard deviations as on a full tumor
geneMap = readRDS(paste0(dataDir, "misc/geneMap.RDS"))
PBraw = readRDS(paste0(dataDir, "PB_count.RDS"))
PBraw = rowsum(PBraw, geneMap$PB[rownames(PBraw)]);
w = rowSums(PBraw)>200;
pb = normIt(PBraw[w,]);
m = rowMeans(pb); s = rowSds(pb); names(s) = rownames(pb);
rm(PBraw, pb);

sigTLS = readRDS(paste0(dataDir, "misc/signature ST.RDS"))[["TLS ST"]]
col = colorRampPalette(c('black', 'green'))(100); mi=0; ma=0.75;

dir.create(paste0(projDir, "/subtype"));
dir.create(paste0(projDir, "/sigTLS"));

ers = mclapply(rownames(cli), function(nm)
{ message(nm);
  img = readRDS(paste0(dataDir, "Robjects/imagesSmall/TNBC", nm, ".RDS"));
  cnt = readRDS(paste0(dataDir, "Robjects/counts/TNBC", nm, ".RDS"));
  id = rep(seq_along(img), sapply(img, function(i) nrow(i$spots)));
  
  x = t(cnt$cnts[,colnames(cnt$cnts)%in%names(m)]);
  x = rowsum(x, geneMap$PB[rownames(x)]);
  x = normIt(x);
  x = (x-m[rownames(x)])/s[rownames(x)];
  bar = TNBCclassif(x, version="bareche", shortName=TRUE, rescale=FALSE)
  
  pdf(paste0(projDir, "/subtype/TNBC", nm, ".pdf"), height=3.5, width=3.5*length(img)+1)
  par(mfrow=c(1, length(img)), cex=.66, omi=c(0,0,0,1))
  for (i in seq_along(img))
  { display(img[[i]]$im, 'r')
    #fa = 1/(max(dim(dta[[i]]$im)[1:2])/238);
    pointsDta(img[[i]], col=colBar[bar[id==i]])
  }
  legend('topright', names(colBar), fill=colBar, inset=c(-.25, .05), xpd=NA, cex=.8)
  dev.off();
  
  x = normIt(t(cnt$cnts));
  sig = calcSig(x, sigTLS);
  pdf(paste0(projDir, "/sigTLS/TNBC", nm, ".pdf"), height=3.5, width=3.5*length(img)+.5);
  par(omi=c(0,0,0,.5))
  plot3(img, sig, col, id);
  plotScale(col, v=c(ma, mi), atV=c(0,1), 5+par('usr')[2], 10, width=1, height=.8)
  dev.off();
}, mc.cores=2)

# Clustering and megaclustering
#################################
# For this one you need the data from "ST TNBC start.R" (source...)

pb = normIt(PBraw[w,]);
m = rowMeans(pb); s = rowSds(pb); names(s) = rownames(pb);
rm(pb);

#m = do.call(rbind, ms);
#scaleFact = nnlm(x=m,y=matrix(1, nrow=nrow(m)))$coef


dir.create(bd<-paste0(projDir, "fullClusts"));
for (nm in rownames(cli))
{ message(nm);
  load(paste0(dataDir, "clustering/intraPatientClust/TNBC", nm, ".RData")) # clusts km qv N scg
  img = readRDS(paste0(dataDir, "Robjects/imagesSmall/TNBC", nm, ".RDS"));
  id = rep(seq_along(img), sapply(img, function(i) nrow(i$spots)));
  mix = readRDS(paste0(dataDir, "clustering/clustPrototypes/TNBC", nm , ".RDS"))$fit$mix;
  bst = which.max(qv-1/N)
  cl = as.integer(km[,bst]); lvl=1:max(cl);
  cols = plt("a")[seq_along(lvl)]
  x = clusts[[bst]]; x = rowsum(x, geneMap$PB[rownames(x)]);
  x = normIt(x[intersect(names(m), rownames(x)),])
  x = (x-m[rownames(x)])/s[rownames(x)];
  barOrig = TNBCclassif(x, version='bar', shortName=TRUE, rescale=FALSE)

  m = m2s[[nm]];
  m2 = m; #m2[m2<1e-2]=0;
  # m2 = m2*rep(scaleFact, each=nrow(m2));
  #m2 = m2/rowSums(m2);
  sM = colSums(m2); wCo=which(sM>1);
  
  pdf(paste0(bd, "/TNBC", nm, ".pdf"), width=3.5*length(img)+2, height=14)
  par(mfcol=c(4, length(img)), xpd=NA, omi=c(0,0,0,2), cex=.65, mar=c(0,0,0,0))
  for (j in seq_along(img))
  { display(img[[j]]$im, 'r')
    pointsDta(img[[j]], col=cols[cl[id==j]])
    display(img[[j]]$im, 'r')
    plotPie(img[[j]]$spots[,"pixel_x"], img[[j]]$spots[,"pixel_y"],
      values=mix[id==j,], radius=2, colors=cols, Nsegments=50)
    display(img[[j]]$im, 'r')
    pointsDta(img[[j]], col=MC.colors[K[klink[paste(nm, lvl)]][cl[id==j]]])
    display(img[[j]]$im, 'r')
    plotPie(img[[j]]$spots[,"pixel_x"], img[[j]]$spots[,"pixel_y"], values=m2[id==j,], radius=2,
      colors=MC.colors, Nsegments=100)

  }
  legend(par('usr')[2]+diff(par('usr')[1:2])*.01, grconvertY(.99, 'ndc', 'user'), fill=cols, xpd=NA,
    legend=paste0("C", lvl, " KM:", barOrig, " Dec:", idC[klink[paste(nm, lvl)],"bar"])) #tnP[paste(nm, lvl)]
  legend(par('usr')[2]+diff(par('usr')[1:2])*.01, grconvertY(.25, 'ndc', 'user'),
    legend=paste0('MC', wCo, " (N spots=", round(sM[wCo]),")"),
    fill=MC.colors[wCo], inset=c(-.45, .01), xpd=NA)
  mtext(paste(c("PB:", "Tumor:", "Stroma:"), sapply(cli[nm,c("barPB","barT.an", "barS.an")], as.character)),
      line=-c(3,5,7), side=1, at=par('usr')[2]*1.01, adj=0, cex=.7)
  dev.off();
}