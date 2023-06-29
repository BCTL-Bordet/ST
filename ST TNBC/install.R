source("~/prog/ngs/R scripts/STscripts.R")

library(xlsx)
library(EnsDb.Hsapiens.v86)
library(EBImage)
library(parallel);options(mc.cores=detectCores()-1);

dataDir = "~/Data/Spatial/TNBC/datadist/";
ids = readRDS(paste0(dataDir, "/Clinical/ids.RDS")) # Array to patients data

## Install counts per spots
#############################
r = mclapply(1:nrow(ids), function(i)
{ message(rownames(ids)[i]);
  bd = paste0(dataDir, "spots/", ids[i,"slide"], "/", ids[i,"subslide"], "/");
  if (file.exists(paste0(bd, "/all.RData"))) { return("Already there"); }
  spotsA = lapply(c(selection="selection.*tsv", all="data-all.*tsv"), function(i)
    read.table(dir(bd, full.names=TRUE, pattern=i), header=TRUE, sep="\t"));
  aa = as.matrix(read.table(paste0(dataDir, "rawCountsMatrices/", ids[i, "count"]), sep="\t", header=TRUE, row.names=1))
  for (nm in names(spotsA))
  { spots = spotsA[[nm]]; cnts = aa;
    rn = paste0(spots$x, "x", spots$y)
    dup = unique(rn[duplicated(rn)])
    if (length(dup)>0)
    { spots = spots[!((rn %in% dup) & spots[,"selected"]==0), ]
      rn = paste0(spots$x, "x", spots$y)
    }
    rownames(spots) = rn;
    spots = spots[rownames(spots) %in% rownames(cnts), ]
    cnts = cnts[rownames(spots), ]
    cnts = cnts[, colSums(cnts)>0]
    save(cnts, spots, file=paste0(bd, nm, ".RData"))
  }
  return("OK")
} , mc.cores=7, mc.preschedule=FALSE)

## Import annotations
######################

## Install data by patient
#############################
dir.create(paste0(dataDir, "/res/dta")); # Will contain all data by patient
dir.create(paste0(dataDir, "/res/cntt")); # Wil contain the count data by patient
dir.create(paste0(dataDir, "/res/imgsOnly")); # Wil contain the images by patient
dir.create(paste0(dataDir, "/res/imgsOnlySmall")); # Wil contain the downscaled images by patient

# transf contains the rotations / translations necessary to superimpose (register) the arrays
transf = read.xlsx(paste0(dataDir, "/misc/registration.xlsx"), 1); rownames(transf) = transf$moving;
transf = transf[!is.na(transf$theta), ]
geneMap = ensembldb::select(EnsDb.Hsapiens.v86, keys=keys(EnsDb.Hsapiens.v86, keytype = "GENEID"), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
g = geneMap[,1]; names(g) = geneMap[,2]; geneMap = g;

autoCut = readRDS(paste0(dataDir, "misc/autocut.RData"))

ers = lapply(unique(ids$id), function(wId)
{ message(wId);
  #if (file.exists(paste0("~/Data/Spatial/TNBC/res/dta/TNBC", wId, ".RData"))) { return("Already there"); }
  
  j = which(ids$id==wId);
  if (any(transf$pts==wId))
  { jj = which(rownames(ids) == transf$fixed[transf$pts==wId][1]); j = c(jj, setdiff(j, jj))
  }

  dta = lapply(j, function(i)
  { message(rownames(ids)[i]);
    r = loadData(paste0(dataDir, "/spots/", ids$slide[i], "/", ids$subslide[i], "/"), geneMap=geneMap, autoCut=autoCut);
    if (nrow(r$cnts) < 20) { message("Problem with",j,"!"); return(NULL); }
    return(r)
  } ); names(dta) = rownames(ids)[j]
  if (is.null(dta[[1]])) { stop("What?"); }
  dta = dta[!sapply(dta, is.null)]

  transfo = function(p, x, ctr)
  { ctr + rep(p[2:3], each=nrow(x)) + ((x-ctr) %*% matrix(c(cos(p[1]), -sin(p[1]), sin(p[1]), cos(p[1])), nrow=2))
  }
  for (i in 1:length(dta))
  { if (!any(rownames(transf)==rownames(ids)[j[i]]))
    { dta[[i]]$im = rotate(dta[[i]]$im, 0, filter='none', output.dim=dim(dta[[i]]$im)[1:2]+600, bg.col='white');
      next;
    }
    tr = transf[rownames(ids)[j[i]],]
    spots = dta[[i]]$spots;
    if (!is.na(tr$inv)) { x = flip(dta[[i]]$im); spots[,"pixel_y"]=dim(x)[1]-spots[,"pixel_y"]; }
    else { x = dta[[i]]$im; }
    spots[,c("pixel_x", "pixel_y")] = transfo(c(tr$theta*base::pi/180, tr$dx, -tr$dy),
      as.matrix(spots[,c("pixel_x", "pixel_y")]), ctr=dim(dta[[i]]$im)[1]/2)
    w = rowAlls(spots[,c('pixel_x', 'pixel_y')]>0);
    x = rotate(x, tr$theta, filter='none', output.dim=dim(dta[[i]]$im)[1:2]+600, bg.col='white')
    x = translate(x, c(tr$dx, -tr$dy), bg.col='white')
    dta[[i]]$spots = spots; dta[[i]]$im = x;
  }
  for (i in seq_along(dta)) { dta[[i]]$imSpot=NULL; }

  # Recenter
  sp = colRanges(do.call(rbind, lapply(dta, function(i) as.matrix(i$spots[,c("pixel_x", "pixel_y")]))));
  sp[,1] = sp[,1]-50; sp[,2]=sp[,2]+50;
  for (i in seq_along(dta))
  { x = translate(dta[[i]]$im, -sp[,1], bg.col='white');
    dta[[i]]$im = x[300+(0:diff(sp[1,])), 300+(1:diff(sp[2,])), ];
    dta[[i]]$spots[,c("pixel_x", "pixel_y")] = dta[[i]]$spots[,c("pixel_x", "pixel_y")] - rep(sp[,1], each=nrow(dta[[i]]$spots))
  }

  # Full cnts
  g = unique(unlist(lapply(dta, function(i) rownames(i$cnts))))
  cntt = do.call(rbind, lapply(dta, function(i) t(fullMat(i$cnts, g))));
 
  i = rep(sapply(dta, function(i) dim(i$cnts)[2]))
  id = rep(seq_along(dta), i);
  
  sp = do.call(rbind, lapply(dta, function(i) i$spots))
  sp$slide = names(dta)[id];
  
  saveRDS(list(cnts=cntt, spots=sp), file=paste0(dataDir, "/countsNonCorrected/TNBC", wId, ".RDS")) 
  
  im = lapply(dta, function(i) list(imgs=i$im, spots=i$spots));
  saveRDS(im, file=paste0(dataDir, "/images/TNBC", wId, ".RDS"))
  
  for (j in seq_along(im))
  { im[[j]]$spots[,c("pixel_x", "pixel_y")] = im[[j]]$spots[,c("pixel_x", "pixel_y")]/5;
    im[[j]]$imgs = EBImage::resize(dta[[j]]$imgs, w = dim(im[[j]]$imgs)[1]/5, antialias=TRUE);
  }
  saveRDS(im, file=paste0(dataDir, "/imagesSmall/TNBC", wId, ".RDS"))
});

## Batch effect correction
############################
library(glmGamPoi)
library(scran);

dir.create(paste0(dataDir, "misc/BatchCorrection")); # Batch correction by patient
dir.create(paste0(dataDir, "counts")); # Batch corrected data (if needed)

# Genes...
PBraw = readRDS(paste0(dataDir, "/PB_count.RDS"))
PB = normRNAseq(PBraw, lim=5e3)
avg = rowSums(PBraw)
g = rownames(PB); rm(PBraw, PB);

NBfits = lapply(dir(paste0(dataDir, "countsNonCorrected/")), function(f)
{ message(f);
  cnt = readRDS(paste0(dataDir, "countsNonCorrected/", f))
  id = factor(cnt$spots$slide);
  x = t(cnt$cnts[, intersect(colnames(cnt$cnts), g)]);
  rm(cntt);
  fit = glm_gp(x, design=~id, size_factors="deconvolution")
  
  b = fit$Beta[,-1,drop=FALSE];
  saveRDS(b, file=paste0(dataDir, "misc/BatchCorrection/",f))
  bl = colMedians(abs(b));
  if (all(bl>10) || !any(bl >.2)) # No batch effect / bogus batch effect, do not correct
  { #file.remove(o<-paste0(dataDir, "/res/cnttCorMine/", f));
    file.symlink(paste0(dataDir, "countsNonCorrected/", f), paste0(dataDir, "counts/", f))
    return(b);
  }
  
  pbs = t(rowsum(t(x), id))
  cc = cor(avg[rownames(pbs)], pbs, method='s')
  cc2 = cor(avg[rownames(pbs)], rowSums(pbs), method='s')
  
  b = cbind(0, b);
  b[rowAnys(abs(b)>10),] = 0;
  if (cc2>max(cc)) { b = b-rowMeans(b); } else { b = b-b[,which.max(cc)]; }
  
  cntt = t(x/exp(b[,id]))
  saveRDS(list(cnts=cntt, spots=cnt$spots), file=paste0(dataDir, "counts/", f))
  return(b);
})

# Make table with correction info
tmp = do.call(rbind,  mclapply(dir(paste0(dataDir, "misc/BatchCorrection/")), function(f)
{ b = readRDS(paste0(dataDir, "misc/BatchCorrection/",f));
  b = fit$Beta[,-1,drop=FALSE];
  bl = colMedians(abs(b));
  
  cnt = readRDS(paste0(dataDir, "countsNonCorrected/", f))
  id = factor(cnt$spots$slide);
  x = t(cnt$cnts[, intersect(colnames(cnt$cnts), g)]);
  
  msg = "Corrected"
  if (all(bl>10) || !any(bl >.2)) { msg = "Not corrected"; }
  
  pbs = t(rowsum(t(x), id))
  cc = cor(avg[rownames(pbs)], pbs, method='s')
  cc2 = cor(avg[rownames(pbs)], rowSums(pbs), method='s')
  
  if (cc2>max(cc)) { base='Mean' } else { base=which.max(cc); }
  
  if (length(bl)==1) { bl = c(bl, ""); cc=c(cc, ""); }
  return(c(msg, base, bl, cc, cc2))
}, mc.preschedule=FALSE))
tbl = data.frame(tmp); for (i in 3:ncol(tbl)) { tbl[,i] = as.numeric(tbl[,i]); }
colnames(tbl) = c("Correction", "Base", "b1", "b2", "cc1", "cc2", "cc3", "cc Mean")
rownames(tbl) = sub("(TNBC[0-9]+).RDS","\\1", dir(paste0(dataDir, "misc/BatchCorrection/")))
#write.xlsx(tbl, file=paste0(dataDir, "misc/batchCor.xlsx"), 'orig') # Commented to protect it

# Info for batch correction
bc = read.xlsx(paste0(dataDir, "misc/batchCor.xlsx"), "modif", row.names=1)
bc$OK.for.actual.version. = sub(" +$", "", bc$OK.for.actual.version.)
bc$final = gsub("[^1-3]+", "", bc$Final.base)


ers = lapply(rownames(bc)[bc$OK.for.actual.version. == "no"], function(f)
{ message("Doing ", f);
  o = paste0(dataDir, "counts/TNBC", f, ".RDS");
  file.remove(o);
  if (bc[f,"Do.additional.correction"] == "batch corr not necessary")
  { file.symlink(paste0(dataDir, "countsNonCorrected/", f), o)
    return("No corr");
  }
  
  b = readRDS(paste0(dataDir, "misc/BatchCorrection/TNBC",f,".RDS"));
  cnt = readRDS(paste0(dataDir, "countsNonCorrected/TNBC", f, ".RDS"))
  id = factor(cnt$spots$slide);
  x = t(cnt$cnts[, intersect(colnames(cnt$cnts), g)]);
  rm(cntt);
  
  base = as.integer(strsplit(bc[f,"final"], "")[[1]]);
  if (length(base)==1)
  { b = fit$Beta[,-1,drop=FALSE];
    b = cbind(0, b);
    b[rowAnys(abs(b)>10),] = 0;
    b = b-b[,base];
    cntt = t(x/exp(b[,id]))
    saveRDS(cntt, file=o);
    return("Cor base 1 slide");
  }
  
  id2 = !(id %in% base);
  fit = glm_gp(x, design=~id2, size_factors="deconvolution")

  b = fit$Beta[,"id2TRUE"]; b[abs(b)>10] = 0;
  cntt = x; cntt[,id2] = cntt[,id2]/exp(b); cntt=t(cntt);
  saveRDS(cntt, file=o);
  
  return("Corrected base 2 slides");
})
