## Helper functions
#########################
symb2Entrez = function(d, tryHard=FALSE)
{ if (length(dim(d)) == 2)
  { rownames(d) = symb2Entrez(rownames(d));
    return(d);
  }
  
  if (!require(org.Hs.eg.db)) { stop("Need library org.Hs.eg.db"); }
  d = factor(d);
  dd = try(suppressMessages(biomaRt::select(org.Hs.eg.db, levels(d), columns="ENTREZID", keytype="SYMBOL")), silent = TRUE);
  if (is(dd, 'try-error'))
  { dd = data.frame(SYMBOL=levels(d), ENTREZID=NA, stringsAsFactors=FALSE);
  }
  dd = dd[!duplicated(dd[,1]),];
  if (tryHard)
  { w = which(is.na(dd[,2]));
    if (length(w)>0)
    { t2 = try(suppressMessages(biomaRt::select(org.Hs.eg.db, dd[w,1], columns="ENTREZID", keytype="ALIAS")), silent = TRUE);
      if (!is(t2, 'try-error'))
      { t2 = t2[!duplicated(t2[,1]),];
        dd[w,2] = t2[,2];
      }
    }
  }
  dd[is.na(dd[,2]),2] = "";
  levels(d) = dd[,2];
  return(as.character(d));
}

entrez2Symb = function(d)
{ if (length(dim(d)) == 2)
  { rownames(d) = entrez2Symb(rownames(d));
    return(d);
  }
  
  if (!require(org.Hs.eg.db)) { stop("Need org.Hs.eg.db"); }
  f = factor(d);
  lvl = suppressMessages(biomaRt::select(org.Hs.eg.db, levels(f), columns="SYMBOL"))[,2];
  lvl[is.na(lvl)] = "";
  levels(f) = lvl;
  if (is.factor(d)) { return(f); }
  else { return(as.character(f)); }
}


## And now about sigs themselves
##################################

calcSig = function(d, sig, dropEmpty=TRUE, balanced=FALSE, isCount=FALSE)
{ if (is(d, "Matrix")) { colMeans=Matrix::colMeans; colSums=Matrix::colSums; }
  if (is.list(sig) && !is.data.frame(sig))
  { ok = rep(TRUE, length(sig)); names(ok) = names(sig);
    w = grep(".norm$", names(sig))
    r1 = list();
    for (i in w)
    { n = sub(".norm$", "", names(sig)[i])
      if (is.null(sig[[n]])) { next; }
      x = calcSig(d, sig[[i]], balanced=balanced, isCount=isCount)
      d2 = d-rep(x, each=nrow(d));
      r1[n] = list(calcSig(d2, sig[[n]], balanced=balanced, isCount=isCount))
      ok[n] = ok[i] = FALSE;
    }
    r1 = do.call(cbind, r1);
    ret = sapply(sig[ok], function(i) calcSig(d, i, balanced=balanced, isCount=isCount));
    if (dropEmpty)
    { ret = ret[,!colSums(!is.na(ret))==0];
    }
    if (!is.null(r1)) { ret=cbind(ret, r1); }
    return(ret);
  }

  s = as.character(sig[,1]);
  s[s==""] = NA;
  x = intersect(rownames(d), s);
  x = x[!is.na(x)];
  if (length(x) == 0)
  { if (colnames(sig)[1] == "name") { other = "entrez"; } else { other = "name"; }
    s = as.character(sig[,other]);
    s[s==""] = NA;
    x = intersect(rownames(d), s);
    x = x[!is.na(x)];
    if (length(x) == 0)
    { warning("No gene in common");
      return(rep(NA, ncol(d)));
    }
  }
  
  sig = sig[match(x, s),];
  
  if (isCount) 
  { fCalc = function(d, x, coef)
    { r = colSums(d[x,,drop=FALSE]*coef, na.rm=TRUE)/colSums(d, na.rm=TRUE);
      sign(r)*log10(1+abs(r)*100);
    }
  } else { fCalc = function(d, x, coef) { colMeans(d[x,,drop=FALSE]*coef); } }
    
  if (balanced && any(!(w <- (sig[,"coefficient"] > 0))))
  { val1 = fCalc(d, x[w], sig[w,"coefficient"]);
    val2 = fCalc(d, x[!w], sig[!w,"coefficient"]);
    val = (val1+val2)/2;
  }
  else
  { val = fCalc(d, x, sig[,"coefficient"]);
  }
  return(val);
}

###################################################################
## Lehmann
TNBCclassif = function(x, version=c("lehmann", "bareche", "burstein", "burstein2", "jiang"), shortName=FALSE,
  coef=FALSE, sig=NULL, rescale=TRUE)
{ version = match.arg(version);
  
  data(signaturesTNBC)
  if (version%in%c("burstein", "burstein2"))
  { if (is.null(sig)) { sig = signaturesTNBC[[version]]; }
    g = intersect(rownames(x), rownames(sig))
    if (length(g)<10) { stop("Less than 10 genes OK, must be in gene names"); }
    if (rescale)
    { y = (rowRanks(x[g,])-1)/(ncol(x)-1);
      #y = t((colRanks(x[g,])-1)/(length(g)-1)); colnames(y) = colnames(x);
    } else
    { if (min(x)<0 || max(x)>1) { stop("x must be between 0 and 1 for burstein"); }
      y = x[g,];
    }
    
    #
    if (version=="burstein")
    { scor = matrix( colMeans( (y[,rep(1:ncol(y), ncol(sig))] -
        sig[g, rep(1:ncol(sig), each=ncol(y))] )^2, na.rm=TRUE ), nrow=ncol(y))
      colnames(scor) = colnames(sig); rownames(scor) = colnames(x);
      scor = sqrt(scor);
      blis_subtype = colnames(scor)[apply(scor,1,function(i)
        {w=which.min(i); if (length(w)==0) {NA} else {w}})]
    } else
    { scor = cor(y, sig[g,], method="spearman")
      blis_subtype = colnames(scor)[apply(scor,1,function(i)
        {w=which.max(i); if (length(w)==0) {NA} else {w}})]
    }
    if (coef) { return(scor); }
    
    return(blis_subtype);
  }
  
  if (rescale) { s = rowSds(x, na.rm=TRUE); y = (x-rowMeans(x, na.rm=TRUE))/s; y=y[which(s>0),]; }
  else { y = x; }
  
  if (version=="jiang")
  { library(pamr);
    fm = signaturesTNBC[["jiang"]]
    ok = !is.na(y[1,]);
    cl = pamr.predict(fm, y[rownames(fm$centroids),ok], threshold=2, type=c("class", "posterior")[coef+1]);
    if (coef) { r = matrix(NA, nrow=ncol(x), ncol=ncol(cl), dimnames=list(colnames(x), colnames(cl))); r[ok,] = cl;}
    else { r = rep(NA, ncol(y)); names(r) = colnames(y); r[ok] = as.character(cl); }
    return(r);
  }
  
  if (is.null(sig)) { sig = signaturesTNBC[["lehmann"]]; }
  if (version=="bareche") { sig$basal_like_2 = NULL; }
  if (shortName) { names(sig) = c(basal_like_1=ifelse(version=="lehmann", "BL1", "BL"), basal_like_2="BL2",
    immunomodulatory="IM",luminal_ar="LAR",
    mesenchymal="M",mesenchymal_stem_like="MSL")[names(sig)] }
  leh = calcSig(y, sig, balanced=TRUE)
  if (coef) { return(leh); }
  cl = colnames(leh)[apply(leh,1,function(i) { r=which.max(i); if (length(r)==1) {r} else { NA; } }) ];
  names(cl) = colnames(y); 
  return(cl);
}

calcTIME= function(d, p1=.6, p2=.454)
{ require(survcomp);
  cs = calcSig(d, getSig("TIME"))

  x = cbind(immune=cs[,"CDSig1"], fib=-cs[,"CDSig3"]);
  x = apply(x,2,function(i) pnorm((i-median(i, na.rm=TRUE))/mad(i, na.rm=TRUE), lower.tail=FALSE))
  x = cbind(x, meta1=apply(x, 1, combine.test, na.rm=TRUE))
  x = cbind(x, meta2=apply(1-x, 1, combine.test, na.rm=TRUE))

  ImmHfibL = x[,"meta1"]<quantile(x[,"meta1"], p1, na.rm=TRUE)
  ImmLfibH = x[,"meta2"]<quantile(x[,"meta2"], 1-p1, na.rm=TRUE)

  ImmHfibL[ImmLfibH] = FALSE;

  x2 = cbind(interferon=cs[,"EDSig2"], fibrosis=-cs[,"EDSig5"])
  x2[ImmLfibH,] = NA;
  x2 = apply(x2,2,function(i) pnorm((i-median(i, na.rm=TRUE))/mad(i, na.rm=TRUE), lower.tail=FALSE))
  x2 = cbind(x2, meta1=apply(x2, 1, combine.test, na.rm=TRUE));
  x2 = cbind(x2, meta2=apply(1-x2, 1, combine.test, na.rm=TRUE));

  IntHCholL = x2[,"meta1"]<quantile(x2[,"meta1"], p2, na.rm=TRUE)
  IntLCholH = x2[,"meta2"]<quantile(x2[,"meta2"], 1-p2, na.rm=TRUE)

  IntHCholL[IntLCholH] = FALSE;

  type = rep(NA, ncol(d));
  type[ImmLfibH] = "Margin restricted";
  type[which(ImmHfibL & IntLCholH)] = "Stroma Restricted"
  type[which(ImmHfibL & IntHCholL)] = "Fully Inflamed"
  return(type);
}