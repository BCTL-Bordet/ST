computeMC = function(x, maxIter=3, rescale=TRUE, mccores=options("mc.cores"), quiet=FALSE)
{ data(MC);
  clustersInOtherDS(x, prot=proto, sigma=sigma, maxIter=maxIter, rescale=rescale,
    mccores=mccores, quiet=quiet)
}

computeET = function(MC)
{ data(ET);
  h = do.call(cbind, lapply(1:ncol(MC), function(i) quantile(faN[,i],
    (rank(MC[,i])-1)/(nrow(MC)-1), type=7)));
  di = knnx.index(log10(faN), log10(h), k=1);
  ecot[di];
}


rescaleForClassif = function(fm, cnts)
{ cnts = cnts[intersect(rownames(cnts), fm$g), ]
  x = 1e4*cnts/rep(colSums(cnts), each=nrow(cnts));
  x2 = matrix(0, nrow=length(fm$center), ncol=ncol(x)); colnames(x2) = colnames(x); rownames(x2) = names(fm$center);
  i = intersect(rownames(x), names(fm$center)); x2[i,] = x[i,]; x = x2; rm(x2);
  x = log(x+1);
  x = (x-fm$center)/fm$scale;
 
  xx = t(qr.solve(fm$rot, x));
}


classifySpots = function(x, type="TNBC")
{ if (type!="TNBC") { stop("Only TNBC at this stage"); }

  a = .libPaths();
  ok = sapply(a, function(i) file.exists(paste0(i, "/STstuff")));
  if (!any(ok)) { stop("For some reason cannot locate the package directory, sorry"); }
  
  a = paste0(a[which(ok)[1]], "/STstuff/extdata/");
  
  base = readRDS(paste0(a, "baseClassif.RDS"));
  
  z = rescaleForClassif(base, x);
  
  d = dir(a, pattern=".xgb");
  fmGlob = lapply(d, function(what) xgb.load(path.expand(paste0(a, what))));
  names(fmGlob) = sub("Reg.xgb", "", d);

  r = do.call(cbind, lapply(fmGlob, function(fm) predict(fm, newdata=z[,1:250])));
}