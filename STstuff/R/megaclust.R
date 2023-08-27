## Extreme prototypes and NB NNMF
###################################

chooseApply = function(mccores)
{ if (is.list(mccores) || mccores==0) { return(function(...) lapply(...)); }
  if (isTRUE(mccores))
  { if (is.list(options("mc.cores"))) { stop("mccores is TRUE but options('mccores') is not defined"); }
    return(function(...) mclapply(..., mc.cores=options("mc.cores")))
  }
  return(function(...) mclapply(..., mc.cores=mccores))
}

## NB nnmf
NNMFfromKM = function(X, km, Niter=3, killList=c(), mccores=0, quiet=FALSE)
{ appl = chooseApply(mccores);

  # Protos
  fitNB = function(x, W, init)
  { toOpt = function(p, x, W)
    { p = exp(p); #print(p[2])
      #-sum(pmax(-20, dnbinom(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=p[1]+.2, log=TRUE)));
      #-sum(dnbinom(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=p[1]+.2, log=TRUE));
      -sum(dnbinom2(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=pmin(p[1]+.2, 1e4)));
    }
    exp(optim(c(0, init), toOpt, x=x, W=W, method="BFGS")$par);
  }

  # Mix
  fitNB2 = function(x, H, sigma)
  { toOpt = function(p, x, H)
    { p = exp(p); #print(p[2])
      #-sum(pmax(-20, dnbinom(x, mu=rowSums(H*rep(p, each=nrow(H))), size=sigma, log=TRUE)));
      #-sum(dnbinom(x, mu=rowSums(H*rep(p, each=nrow(H))), size=sigma, log=TRUE));
      -sum(dnbinom2S(x, mu=rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4)));
    }
    exp(optim(rep(0, ncol(H)), toOpt, x=x, H=H, method="BFGS")$par);
  }


  Ncl = max(km); cl = 1:Ncl; if (length(killList)>0) { cl = cl[-killList]; }
  W = matrix(0, nrow=nrow(X), ncol=length(cl)); colnames(W) = cl;
  w = km %in% cl;
  W[cbind(which(w), match(km[w], cl))] = 1;
  
  #W[!w,] = 1/ncol(W);
  H = matrix(colMeans(X), nrow=ncol(X), ncol=ncol(W));

  Hs = list(); prevMix=list();
  for (iter in 1:Niter)
  { if (!quiet) { message("     Iter:", iter); }
    tmp = do.call(rbind, appl(1:ncol(X), function(i)
    { fitNB(X[w,i], W[w,], log(pmax(.1, H[i,])))
      #fitNB(X[,i], W, log(H[i,]))
    }))
    w = rep(TRUE, nrow(X));
    H = tmp[,-1]; sigma = tmp[,1]+.2;
    Hs[[iter]] = H;

    W = do.call(rbind, appl(1:nrow(X), function(i)
    { fitNB2(X[i,], H, sigma)
    }));
    prevMix[[iter]] = W;
  }
  return(list(proto=H, prevProto=Hs, mix=W, prevMix=prevMix, sigma=sigma, cl=cl));
}

#expressionInAnnotations = function(X, mix, mccores=0)
#{ appl = chooseApply(mccores);
#
#  fitNB = function(x, W)
#  { toOpt = function(p, x, W)
#    { p = exp(p);
#      -sum(dnbinom2(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=pmin(p[1]+.2, 1e4)));
#    }
#    exp(optim(c(0, rep(log(mean(x)), ncol(W))), toOpt, x=x, W=W, method="BFGS")$par);
#  }
#  X = X[,colSums(X>0)>5];
#  tmp = do.call(rbind, appl(1:ncol(X), function(i)
 # { fitNB(X[,i], mix)
#  }))
#  rownames(tmp) = colnames(X);
#  proto = tmp[,-1]; sigma = tmp[,1]+.2;
#  return(list(proto=proto, sigma=sigma));
#}

allLikelihood = function(X, W, H, sigma)
{ mu = W %*% t(H);
  p = dnbinom(X, mu=mu, size=rep(sigma, each=nrow(X)), log=TRUE);
}


## Do the deconvolution
deconvoluteClusters = function(X, k, clean=TRUE, rmCutOff=0, quiet=FALSE, killK = c(), Niter=10,
  mccores=options("mc.cores"), anchor=1e3)
{ if (clean)
  { X = X[, colSums(X)>0]
    c2 = X[,colMeans(X>2)>.05];
  } else { c2 = X; }
  
  c3 = scale(log(1+c2/rowSums(c2)*1e4))
  proto = rowsum(c3, k)/as.vector(table(k))

  if (ncol(X) > 500)
  { withins = colSums( (c3-proto[k,])^2 ) / nrow(c3)
    o = order(withins)[1:500];
    X = c2[,o];
  }
  
  x = rowsum(c2, k);
  killK = c(killK, which(rowSums(x)<1e4))

  prevProt = list();
  while (TRUE)
  { if (!quiet) { message("  Kill: ", paste(killK, collapse=",")) };
    aa = NNMFfromKM2(X, k, Niter=Niter, killList=killK, mccores=mccores, quiet=quiet, anchor=anchor)
    prevProt = c(prevProt, list(aa));
    mix = aa$mix; mix = mix/rowSums(mix);
    ma = colQuantiles(mix, probs=.98);
    if (any(ma<rmCutOff-.2)) { killK = c(killK, aa$cl[rmCutOff-.2]); next; }
    if (any(ma<rmCutOff)) { killK = c(killK, aa$cl[which.min(ma)]); next; }
    break;
  }
  if (!quiet) { message("  Full prototypes"); }
  proto = NNMFproto(c2, aa$mix, mccores=mccores)
  colnames(proto[[1]]) = aa$cl
  return(list(proto=proto, prevProt=prevProt, fit=aa, kOrig=k, Nreads=rowSums(x)));
}

NNMFproto = function(X, mix, mccores=0)
{ appl = chooseApply(mccores);

  fitNB = function(x, W)
  { toOpt = function(p, x, W)
    { p = exp(p); #print(p[2])
      #-sum(pmax(-20, dnbinom(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=p[1]+.1, log=TRUE)));
      -sum(dnbinom(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=p[1]+.2, log=TRUE));
    }
    exp(optim(c(0, rep(log(mean(x)), ncol(W))), toOpt, x=x, W=W, method="BFGS")$par);
  }
  X = X[,colSums(X>0)>5];
  tmp = do.call(rbind, appl(1:ncol(X), function(i)
  { fitNB(X[,i], mix)
  }))
  rownames(tmp) = colnames(X);
  proto = tmp[,-1]; sigma = tmp[,1];
  return(list(proto=proto, sigma=sigma));
}


## NB nnmf
NNMFfromKM2 = function(X, km, Niter=3, killList=c(), mccores=0, quiet=FALSE, anchor=1)
{ appl = chooseApply(mccores);

  # Protos
  fitNB = function(x, W, init)
  { toOpt = function(p, x, W)
    { p = exp(p); #print(p[2])
      #-sum(pmax(-20, dnbinom(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=p[1]+.2, log=TRUE)));
      #-sum(dnbinom(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=p[1]+.2, log=TRUE));
      -sum(dnbinom2(x, mu=rowSums(W*rep(p[-1], each=nrow(W))), size=pmin(p[1]+.2, 1e4)));
    }
    exp(optim(c(0, init), toOpt, x=x, W=W, method="BFGS")$par);
  }

  # Mix
  fitNB2 = function(x, H, sigma, mixO, anchor)
  { toOpt = function(p, x, H)
    { p = exp(p); #print(p[2])
      cc = suppressWarnings(cor(mixO, p)); if (!is.finite(cc)) { cc=0; }
      #-sum(pmax(-20, dnbinom(x, mu=rowSums(H*rep(p, each=nrow(H))), size=sigma, log=TRUE)));
      #-sum(dnbinom(x, mu=rowSums(H*rep(p, each=nrow(H))), size=sigma, log=TRUE));
      -sum(dnbinom2S(x, mu=rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4))) +
        anchor*length(mixO)*(1-cc)^2; #anchor*sum((mixO-p)^2);
    }
    exp(optim(rep(-.5, length(mixO)), toOpt, x=x, H=H, method="BFGS")$par);
  }


  Ncl = max(km); cl = 1:Ncl; if (length(killList)>0) { cl = cl[-killList]; }
  W = matrix(0, nrow=nrow(X), ncol=length(cl)); colnames(W) = cl;
  w = km %in% cl;
  W[cbind(which(w), match(km[w], cl))] = 1;
  mixO = W; # Mix based on kmeans
  
  #W[!w,] = 1/ncol(W);
  H = matrix(colMeans(X), nrow=ncol(X), ncol=ncol(W));

  Hs = list(); prevMix=list();
  for (iter in 1:Niter)
  { if (!quiet) { message("     Iter:", iter); }
    tmp = do.call(rbind, appl(1:ncol(X), function(i)
    { fitNB(X[w,i], W[w,], log(pmax(.1, H[i,])))
      #fitNB(X[,i], W, log(H[i,]))
    }))
    w = rep(TRUE, nrow(X));
    H = tmp[,-1]; sigma = tmp[,1]+.2;
    Hs[[iter]] = H;

    W = do.call(rbind, appl(1:nrow(X), function(i)
    { fitNB2(X[i,], H, sigma, mixO[i,], anchor)
    }));
    prevMix[[iter]] = W;
  }
  return(list(proto=H, prevProto=Hs, mix=W, prevMix=prevMix, sigma=sigma, cl=cl));
}

####################################################

dnbinom2 = function(x, mu, size)
{ #lgamma(size+x)-lfactorial(x)-lgamma(size) + size*log(size/(size+mu)) + x*log(mu/(mu+size))
  #-lbeta(size, x) - log(x) + size*log(size/(size+mu)) + x*log(mu/(mu+size))
  -lbeta(size, x+1) - log(x+size) + size*log(size/(size+mu)) + x*log(mu/(mu+size))
}

dnbinom2S = function(x, mu, size) # Version if keeping size fixed
{ #lgamma(size+x)-lfactorial(x)-lgamma(size) + size*log(size/(size+mu)) + x*log(mu/(mu+size))
  #-lbeta(size, x) - log(x) + size*log(size/(size+mu)) + x*log(mu/(mu+size))
  size*log(size/(size+mu)) + x*log(mu/(mu+size))
  #log( ((size/(size+mu))^size) * ((mu/(mu+size))^x) )
}

# NB version
clustersInOtherDS = function(x, prot, sigma, maxIter=3, rescale=TRUE, mccores=options("mc.cores"), quiet=FALSE)
{ appl = chooseApply(mccores);

  fitNB2 = function(x, H, sigma, sc, init=rep(1, ncol(H)))
  { toOpt = function(p, x, H, sc, sigma)
    { p = exp(p);
      #-sum(pmax(-20, dnbinom2(x, mu=sc*rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4))));
      -sum(dnbinom2(x, mu=sc*rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4)));
    }
    opt = optim(log(init), toOpt, x=x, H=H, sc=sc, sigma=sigma, method="BFGS");
    c(opt$value, exp(opt$par));
  }
  fitNB2S = function(x, H, sigma, sc, init=rep(1, ncol(H)))
  { toOpt = function(p, x, H, sc, sigma)
    { p = exp(p);
      #-sum(pmax(-20, dnbinom2(x, mu=sc*rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4))));
      -sum(dnbinom2S(x, mu=sc*rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4)));
    }
    opt = optim(log(init), toOpt, x=x, H=H, sc=sc, sigma=sigma, method="BFGS");
    c(opt$value, exp(opt$par));
  }

  fitSigma2 = function(x, prot, W, s)
  { toOpt = function(i, x, prot, W)
    { i = exp(i);
      #-sum(pmax(-20, dnbinom2(x, mu=i[2]*rowSums(W*rep(prot, each=nrow(W))), size=min(i[1],1e4))));
      -sum(dnbinom2(x, mu=i[2]*rowSums(W*rep(prot, each=nrow(W))), size=min(i[1],1e4)));
    }
    opt=optim(log(s), toOpt, x=x, prot=prot, W=W, method="BFGS")
    c(opt$value, exp(opt$par));
  }
  
  fitSigma = function(x, prot, W, s)
  { toOpt = function(i, x, prot, W)
    { i = exp(i);
      #-sum(pmax(-20, dnbinom2(x, mu=rowSums(W*rep(prot, each=nrow(W))), size=min(i[1],1e4))));
      -sum(dnbinom2(x, mu=rowSums(W*rep(prot, each=nrow(W))), size=min(i[1],1e4)));
    }
    opt=optim(log(s), toOpt, x=x, prot=prot, W=W, method="BFGS")
    c(opt$value, exp(opt$par));
  }

  # Init
  g = intersect(rownames(x), rownames(prot)); x=x[g,]; prot=prot[g,];
  s = sigma[rownames(prot)]; # Sigma per gene
  s[s>100] = 100; s[s<.5] = .5;
  
  W = matrix(colSums(x), nrow=ncol(x), ncol=ncol(prot))/sum(prot); # Init W
  #W = t(nnlm(prot, x, loss='mkl', n.threads = options("mc.cores"))$coefficients);
  #W[W<1e-2] = 1e-2;
  
  if (rescale) { sc = rowMeans(x)/ rowMeans(prot %*% t(W)); sc[sc<.1]=.1; sc[sc>10]=10; } # Scaling factor 
  else { sc = rep(1, nrow(prot)); }
  
  Wb = 0;

  iter = 1; qb=+Inf; Wp = list();
  while((m<-max(abs(W-Wb)/(W+Wb+1)))>1e-3 && iter<maxIter)
  { if (!quiet) { message("Iter: ", iter, " rel diff: ", signif(m, 4)); }
    Wb=W; 
    # Update W
    
    #W1 = do.call(rbind, mclapply(1:ncol(x), function(i)
    #{ fitNB2(x[,i], prot, s, sc, init=W[i,]) }));
    W = do.call(rbind, appl(1:ncol(x), function(i)
    { fitNB2S(x[,i], prot, s, sc, init=W[i,]) }));
    W[,1] = W[,1] -colSums(-lbeta(s, x+1) - log(x+s))
    qual<-sum(W[,1])
    if (!quiet) { message("  W qual:   ", qual) }
    Wp[[iter]] = W = W[,-1];
    #if (qual>qb) { browser(); }
    qb=qual;
    
    iter=iter+1; if (iter==maxIter) { break; } # No point of redoing sigma for last iter
    # Update sigma and scaling
    if (rescale)
    { tmp = do.call(rbind, appl(1:nrow(x), function(i) fitSigma2(x[i,], prot[i,], W, c(s[i], sc[i]))))
      sc = tmp[,3]
    }
    else { tmp = do.call(rbind, appl(1:nrow(x), function(i) fitSigma(x[i,], prot[i,], W, s[i]))) }
    s=tmp[,2]; s[s>1e4-1] = 1e4-1;
    qual<-sum(tmp[,1])
    if (!quiet) { message("  Sisc qual:", qual) }
    if (qual>qb) { browser(); }
    qb=qual;
  }
  names(s) = names(sc) = g;
  return(list(W=W, Wp=Wp, sigma=s, sc=sc, qual=qb))
}


expressionInAnnotations = function(x, f, sigma, maxIter=5, mccores=TRUE, quiet=FALSE)
# f: matrix of predictions, rows: spots; cols: components
# x: matrix of counts, rows: spots; cols: genes
{ fitNB2 = function(x, f, sc, sigma, init=rep(1, ncol(f)))
  { toOpt = function(p, x, f, sc)
    { p = exp(p); sigma=p[1]+.2; p=p[-1];
      #-sum(pmax(-20, dnbinom2(x, mu=sc*rowSums(H*rep(p, each=nrow(H))), size=pmin(sigma, 1e4))));
      -sum(dnbinom2(x, mu=rowSums(sc*f*rep(p, each=nrow(f))), size=pmin(sigma, 1e3)));
    }
    opt = optim(log(c(pmax(sigma-.2, .01), init)), fn=toOpt, x=x, f=f, sc=sc, method="BFGS");
    r = exp(opt$par); r[1] = r[1]+.2;
    c(opt$value, r);
  }
  
  fitSc = function(x, f, W, sc, sigma)
  { toOpt = function(p, x, f, sigma, W)
    { sc = exp(p);
      -sum(dnbinom2S(x, mu=sc*colSums(f*W), size=sigma));
    }
    opt = optim(log(sc), fn=toOpt, x=x, f=f, sigma=sigma, W=W, method="BFGS");
    c(opt$value, exp(opt$par));
  }
  
  fitSigma = function(x, f, W, sc, s)
  { toOpt = function(i, x, f, W, sc)
    { i = exp(i);
      #-sum(pmax(-20, dnbinom2(x, mu=rowSums(W*rep(prot, each=nrow(W))), size=min(i[1],1e4))));
      -sum(dnbinom2(x, mu=sc*rowSums(f*rep(W, each=nrow(f))), size=min(i[1],1e4)));
    }
    opt=optim(log(s), fn=toOpt, x=x, f=f, W=W, sc=sc, method="BFGS")
    c(opt$value, exp(opt$par));
  }

  appl = chooseApply(mccores);

  # Init
  if (nrow(x) != nrow(f)) { stop("Not the same set of spots in x and f"); }
  #g = intersect(rownames(x), rownames(f)); x=x[g,]; f=f[g,];
  if (missing(sigma)) { s = rep(1, ncol(x)) } else { s = sigma; } 
  s[s>100] = 100; s[s<.5] = .5;
  
  W = matrix(colSums(x), ncol=ncol(x), nrow=ncol(f), byrow=TRUE)/sum(f); # Init W (components, col=genes)
  sc = rowMeans(x)/ rowMeans(f %*% W); sc[sc<.1]=.1; sc[sc>10]=10; # Scaling factor 
  
  Wb = 0; Winit=W;

  iter = 1; qb=+Inf; Wp = sp = scp = list();
  while(iter<=maxIter) #(m<-max(abs(W-Wb)/(W+Wb+1)))>1e-3 && )
  { if (!quiet) { message("Iter: ", iter) }#, " rel diff: ", signif(m, 4)); }
    Wb=W; 
    # Update W
    
    #W1 = do.call(rbind, mclapply(1:ncol(x), function(i)
    #{ fitNB2(x[,i], prot, s, sc, init=W[i,]) }));
    W = do.call(cbind, appl(1:ncol(x), function(i)
    { fitNB2(x[,i], f=f, sigma=s[i], sc=sc, init=Winit[,i]) }));
    #W[1,] = W[1,] -colSums(-lbeta(s, x+1) - log(x+s))
    qual<-sum(W[1,])
    s = W[2,]; s[s>1e3-10] = 1e3-10;
    sp[[iter]]=s;
    if (!quiet) { message("  W qual:    ", qual) }
    Wp[[iter]] = W = W[-(1:2),];
    #if (qual>qb) { browser(); }
    qb=qual;
    
    iter=iter+1; if (iter==maxIter) { break; } # No point of redoing sigma for last iter
    
    # Update sigma
    #tmp = do.call(rbind, appl(1:ncol(x), function(i) fitSigma(x[,i], f, W[,i], sc, s[i])))
    #s=tmp[,2]; s[s>1e4-1] = 1e4-1;
    #qual<-sum(tmp[,1])
    #if (!quiet) { message("  Sigma qual:", qual) }
    #if (qual>qb) { browser(); }
    #qb=qual;
    
    # Update scaling
    tmp = do.call(rbind, appl(1:nrow(x), function(i) fitSc(x[i,], f[i,], W, sc[i], s)))
    sc=scp[[iter]] = tmp[,2];
    qual<-sum(tmp[,1])
    if (!quiet) { message("  SC qual:   ", qual) }
  }
  W = t(W); colnames(W) = colnames(f); rownames(W) = colnames(x);
  names(s) = colnames(x); names(sc) = rownames(x);
  return(list(W=W, Wp=Wp, sigma=s, sc=sc, qual=qb, sp=sp, scp=scp))
}