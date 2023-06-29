plotTable = function(x, colHeader='grey80', space=.01, colCell=c('grey96', 'grey90'),
  colSize=NULL, rowSize=NULL, tableLabels=NULL, showZeros=TRUE, parseCells=FALSE)
{ plot.new();
  rn = !is.null(rownames(x)); cn = !is.null(colnames(x));
  
  if (is.null(colSize)) { co = rep(1, ncol(x)+rn); } else
  { if (length(colSize) != ncol(x)+rn) { stop("Expected length(colSize)==", ncol(x)+rn); }
    co = colSize;
  }
  if (is.null(rowSize)) { ro = rep(1, nrow(x)+cn); } else
  { if (length(rowSize) != nrow(x)+cn) { stop("Expected length(rowSize)==", nrow(x)+cn); }
    ro = rowSize;
  }
  a = c(0, cumsum(co)); coL = a[-length(a)]+space; coR = a[-1]-space; coC = (coL+coR)/2; 
  a = c(0, cumsum(ro)); roL = a[-length(a)]+space; roR = a[-1]-space; roC = (roL+roR)/2;
  
  xx = (1:ncol(x))+rn; yy = (1:nrow(x))+cn;
  
  plot.window(xlim=c(coL[1],coR[length(coR)]), ylim=c(roR[length(roR)], roL[1]))
  
  if (!is.null(colHeader))
  { if (cn) { rect(coL[xx], roL[1], coR[xx], roR[1], col=colHeader, border=NA); }
    if (rn) { rect(coL[1], roL[yy], coR[1], roR[yy], col=colHeader, border=NA); }
  }
  if (!is.null(colCell))
  { a = rep(colCell, ceiling(nrow(x)/2))[1:nrow(x)];
    rect(rep(coL[xx], each=nrow(x)), roL[yy], rep(coR[xx], each=nrow(x)),
      roR[yy], col=a, border=NA);
  }
  if (cn) { text(coC[xx], roC[1], labels=colnames(x), font=2) }
  if (rn) { text(coC[1], roC[yy], labels=rownames(x), font=2) }

  if (!showZeros) { w = x != 0; } else { w = rep(TRUE, length(x)); }
  if (parseCells) { cls = parse(text=x[w]); } else { cls=x[w]; }
  text(x=rep(coC[xx], each=nrow(x))[w], rep(roC[yy], ncol(x))[w], cls);
  
  if (!is.null(tableLabels))
  { if (identical(tableLabels, TRUE) && is.table(x)) { tableLabels=names(dimnames(x)); }
    mtext(tableLabels, side=2:3, line=0, cex=par('cex'), font=2, at=c((roL[yy[1]]+roR[yy[length(yy)]])/2,
      (coL[xx[1]]+coR[xx[length(xx)]])/2 ))
  } 
}

basicForest = function(x, a, adj=NULL, xlim=NULL, xlog=FALSE, xlab="", col='black', cex.axis=.7,
  titles=names(x), annotDir=NULL, lineHeight=1.5, cexAnnot=0.7, colWidth=NULL, pch=16)
{ if (is.data.frame(x)) { nr = nrow(x); } else { nr = length(x[[1]]); }
  if (nr != nrow(a)) { stop("Not same number of rows in x and a"); }
  if (is.null(adj)) { adj = c("left", rep("right", length(x)-1)) }
  if (length(col)==1) { col = rep(col, nr); }
  
  plot.new(); plot.window(xlim=c(0,1), ylim=c(1, 0), xaxs='i', yaxs='i') 
  
  if (is.null(colWidth))  
  { cw = colMaxs(sapply(x, function(i) sapply(i, strwidth)));
    cw = pmax(cw, strwidth(names(x)));
    cw = cw + strwidth("M");
  }
  else { cw = colWidth[-length(colWidth)]/sum(colWidth); } 
   
  ch = abs(strheight("!"));
  if (is.numeric(lineHeight))
  { if (length(lineHeight)==1) { li = (1:(1+nr))*ch*lineHeight; }
    else { li = lineHeight*ch; }
  }
  else
  { if (lineHeight=="full")
    { li = seq(ch, 1-7*ch, len=nr+1); }
  }
  
  if (max(li)>1) { warning("Too many lines, does not fit in window"); }
  st = cumsum(c(0,cw));
  left = st[-length(st)]; right=st[-1]; center=(left+right)/2;
  cpos = left; cpos[adj=="right"]=right[adj=="right"]; cpos[adj=="center"]=center[adj=="center"];
  
  li2 = li-strheight("M")*.4;
  ad = c(left=0, right=1, center=.5)[adj]
  for (i in seq_along(x))
  { for (j in seq_along(x[[i]]))
    { text(labels=x[[i]][[j]], x=cpos[i], y=li2[j+1], adj=c(ad[i], 0)); }
      #rep(cpos[i], each=nr), y=li[-1], adj=ad[i]) }
    text(labels=titles[[i]], x=cpos[i], y=li2[1], font=2, adj=c(ad[i], 0))
  }
  
  if (is.null(xlim)) { xlim=range(a, na.rm=TRUE) };
  xx = grconvertX(c(st[length(st)]+1.5*strwidth("M"), 1), from='user', to='nfc')
  yy = grconvertY(c(li[length(li)],li[2]), from='user', to='nfc')
  old.par <- par( c('plt', 'usr', 'mgp', "cex.axis") )
  on.exit(par(old.par))
  
  par(new=TRUE, plt=c(xx,yy), mgp=c(1.5,.5,0), cex.axis=cex.axis);
  xl = range(a, na.rm=TRUE);
  plot.new(); plot.window(xlim=xlim, ylim=c(1,0), yaxs='i', xaxs='i', log=ifelse(xlog, "x", ""))
  axis(1, line=1); title(xlab=xlab, line=2.5, cex.lab=cex.axis);
  lines(c(xlog, xlog), c(1-strheight("M")*1.5, 0+strheight("M")), col='lightgrey', xpd=TRUE);
  
  wp = which(!is.na(a[,2]));
  ty=rep(0, nrow(a));
  ty[which(a[,1]<xlim[1])]=ty[which(a[,1]<xlim[1])]+1;
  ty[which(a[,3]>xlim[2])]=ty[which(a[,3]>xlim[2])]+2;
  
  le=.05;
  w = which(a[,2]>=xlim[1] & a[,2]<=xlim[2]);
  a[which(a[,2]<xlim[1]),2] = (xlim[1]+a[which(a[,2]<xlim[1]),3])/2
  a[which(a[,2]>xlim[2]),2] = (xlim[2]+a[which(a[,2]>xlim[2]),1])/2
  
  li2 = li[-1]-li[2]; li2=li2/li2[length(li2)];
  
  for (i in wp)
  { suppressWarnings(arrows(a[i,2], li2[i], min(xlim[2], a[i,3]), li2[i], angle=c(30,90)[(a[i,3]<=xlim[2])+1],
      length=le, col=col[i], xpd=NA, lwd=2*par('cex')))
    suppressWarnings(arrows(a[i,2], li2[i], max(xlim[1], a[i,1]), li2[i], angle=c(30,90)[(a[i,1]>=xlim[1])+1],
      length=le, col=col[i], xpd=NA, lwd=2*par('cex')))
  }
  points(a[w,2], li2[w], pch=pch, col=col[w], xpd=NA, cex=2, family="Arial Unicode MS")
  
  if (!is.null(annotDir))
  { h = 1-strheight("M")*5;
    v = strwidth(xlab, cex=cex.axis)*2;
    med = mean(log(xlim));
    arrows(exp(med-v), h, xlim[1], h, xpd=TRUE, length=le)
    arrows(exp(med+v), h, xlim[2], h, xpd=TRUE, length=le)
    text(annotDir, x= sqrt(c(exp(med-v), exp(med+v))*xlim), y=h-strheight("M"), xpd=NA, cex=cexAnnot);
  }
}

colsAnnot = function(what=c("pam50", "IHC", "IntSub", "SCM", "leh", "all"))
{ what=match.arg(what);
  colPam50 = c(Basal="#E53935", Her2="#F48FB1", LumA="#3949AB", LumB="#03A9F4", Normal="#4CAF50")#, Bof="black"); 
  colIHC = colPam50[c("Basal", "Her2", "LumB")]; names(colIHC) = c("TNBC", "HER2+", "HR+/HER2-")
  colIntSub = colIHC; names(colIntSub) = c("ER-/HER2-", "HER2+", "ER+/HER2-")
  colSCM = colPam50[c("Basal", "Her2", "LumA", "LumB")]; names(colSCM) = c("ER-/HER2-", "HER2+", "ER+/HER2- Low Prolif", "ER+/HER2- High Prolif")
  colLeh = c(BL1="#1976D2", BL2="#4fc3f7", IM="#66BB6A", LAR="#EF5350", M="#FFEE58", MSL="#FFA726", UNS="#8D6E63")

  r = list(pam50=colPam50, IHC=colIHC, IntSub=colIntSub, SCM=colSCM, leh=colLeh);
  if (what=="all") { return(r); }
  return(r[[what]])
}

addAnnot = function(name, cols, line, side=1, at=seq_along(cols), heights=1, pch=NA, point.col=NA, point.cex=1)
{ if (side %in% c(1,3))
  { lineW = strheight("M\nM")-strheight("M");
    if (side==1) { y0 = par('usr')[3]-(line)*lineW; }
    else { y0 = par('usr')[4]+(line)*lineW; };
  } else 
  { lineW = strwidth("M"); # Not very satisfying...
    if (side==2) { y0 = par('usr')[1]-(line)*lineW; }
    else { y0 = par('usr')[2]+(line)*lineW; };
  }
  
  a = c(at[1]-(at[2]-at[1]), at, at[length(at)]+(at[length(at)]-at[length(at)-1]))
  x = (a[-1] + a[-length(a)])/2;
  if (side %in% c(1,3)) { rect(x[-length(x)], y0, x[-1], y0-heights*lineW, xpd=NA, border=NA, col=cols) }
  else { rect(y0, x[-length(x)], y0-heights*lineW, x[-1],  xpd=NA, border=NA, col=cols) }
  if (any(!is.na(pch)))
  { if (side %in% c(1,3)) { points(at, rep(y0-heights*lineW/2, length(at)), pch=pch, col=point.col, xpd=NA, cex=point.cex);  }
    else { points(rep(y0-heights*lineW/2, length(at)), at, pch=pch, col=point.col, xpd=NA, cex=point.cex);  }
  }
  mtext(name, side=side, line=line, at=a[1], adj=1)
}

formatNdig = function(x, n=2)
{ m = floor(log10(x))+1;
  if (m>=n) { return(format(x, digits=1)); }
  #m2 = m; if (m2<0) { m2 = 0; }
  format(x, nsmall=n-m, digits=n)
}

formatNice = function(x, parse=TRUE, lim=1e-20, nDigits=2)
{ if (!is.finite(x))
  { if (is.na(x)) { return("NA"); }
    if (x==Inf) { return("Inf"); }
    if (x==-Inf) { return("-Inf"); }
    return("NA");
  }
  if (x < lim) { return("0"); }
  if (abs(x)>=1e-3 && abs(x)<=1e3) { return(formatNdig(x, n=nDigits)); }
  ex = ceiling(-log10(x));
  ma = formatNdig(x*10^ex, n=nDigits);
  #ret = paste(ma,"%*%10^{",-ex,"}");
  if (ma=="1") { ret=paste0("10^{",-ex,"}")}
  else { ret = paste0("paste('", ma, "x', 10^{",-ex,"})"); }
  if (parse) { return(parse(text=ret)); } 
  return(ret);
}

formatNiceP = function(x, descr=NULL, parse=TRUE, lim=1e-20, pName="p", sepDescr=" - ", italic=TRUE, boldDescr=FALSE)
{ p = formatNice(x, parse=FALSE, lim=lim);
  if (italic) { pName = paste0("italic(", pName, ")"); }
  if (!grepl("\\^",p)) { p = paste0("'", p, "'")}
  if (is.null(descr)) { ret = paste0("paste(", pName, ",' = ', ", p, ")"); }
  else
  { if (boldDescr) { descr = paste0("bold('", descr, "')"); } else { descr = paste0("'", descr, "'")}
    ret = paste0("paste(", descr, ",'", sepDescr, "', ", pName, ",' = ', ", p, ")");
  }
  if (parse) { return(parse(text=ret)); } 
  return(ret);
}

fpDrawCI.old = function (lower_limit, estimate, upper_limit, size, y.offset = 0.5, 
    clr.line, clr.marker, lwd, lty = 1, vertices, vertices.height = 0.1, 
    ...) 
{
    forestplot:::prFpDrawLine(lower_limit = lower_limit, upper_limit = upper_limit, 
        clr.line = clr.line, lwd = lwd, lty = lty, y.offset = y.offset, 
        vertices = vertices, vertices.height = vertices.height)
    box <- convertX(unit(estimate, "native"), "npc", valueOnly = TRUE)
    skipbox <- box < 0 || box > 1
    if (!skipbox) {
        if (!is.unit(size)) {
            size <- unit(size, "snpc")
        }
        grid.rect(x = unit(estimate, "native"), y = y.offset, 
            width = size, height = size, gp = gpar(fill = clr.marker, 
                col = clr.line));
    }
}

fpDrawCI = function (lower_limit, estimate, upper_limit, size, y.offset = 0.5, 
                     clr.line, clr.marker, lwd = 2, lty = 1, vertices, vertices.height = 0.1, shapes_gp = fpShapesGp(), 
                     shape_coordinates = structure(c(1, 1), max.coords = c(1, 1)),
                     ...) 
{
  forestplot:::prFpDrawLine(lower_limit = lower_limit, upper_limit = upper_limit, 
                            clr.line = clr.line, lwd = lwd, lty = lty, y.offset = y.offset, 
                            vertices = vertices, vertices.height = vertices.height)
  box <- convertX(unit(estimate, "native"), "npc", valueOnly = TRUE)
  skipbox <- box < 0 || box > 1
  if (!skipbox) {
    if (!is.unit(size)) {
      size <- unit(size, "cm")
    }
    grid.circle(x = unit(estimate, "native"), y = y.offset, r = size,                                              # this will plot circles instead of squares
                gp = prGetShapeGp(shapes_gp, shape_coordinates, "box", default = gpar(fill = clr.marker,  
                                                                                      col = clr.marker)))
  }
}

allForest = function(x, y, subset=NULL, control=~1, sig=.05, sigOnFDR=FALSE, fdr=TRUE, new_page=FALSE, clip=c(.2,5),
  ticks=NULL, dispN=FALSE, useWilcox=FALSE, genesItal=FALSE, parseNames=FALSE, rmNAp=TRUE, normReal=TRUE,
  capHR=TRUE, weights=NULL, useNested=TRUE, columns=c("OR", "CI", "P", "FDR"), boxsize=.2, 
  version=3, FirthIfNeeded=TRUE, graphwidth = unit(4, "cm"), plot=TRUE, ...)
{ if (!require(forestplot)) { stop("Need forestplot"); }
  f = formatNice;
  withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
        myWarnings <<- c(myWarnings, list(w))
        invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
  } 


  if (version==1)
  { plBlue = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...)
    { fpDrawCI.old(lower_limit, estimate, upper_limit, size, y.offset=.5,
        clr.line="darkblue", clr.marker="royalblue1");
    }
    plYellow = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...)
    { fpDrawCI.old(lower_limit, estimate, upper_limit, size, y.offset=.5,
        clr.line="yellow4", clr.marker="orange");
    }
  }
  if (version==2)
  { plBlue = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...)
    { fpDrawCI(lower_limit, estimate, upper_limit, size, y.offset=.5, 
               # clr.line="darkblue", 
               clr.line="#00AFBB",       # "#00AFBB" , "royalblue1                ################# In this three functions (plBlue, plYellow, plRed) you can change the colors.
               clr.marker="#00AFBB"
               ,vertices = TRUE);
    }
    plYellow = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...) # not yellow anymore, but same name of the function
    { fpDrawCI(lower_limit, estimate, upper_limit, size, y.offset=.5, 
               clr.marker = "slategray4", # "#6C7B8BFF", #alpha("slategray4", 1), # "orange" , "darkgrey"
               clr.line="slategray4"
               # clr.line="yellow4" 
               ,vertices = TRUE)
    }
    plRed = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...)
    { fpDrawCI(lower_limit, estimate, upper_limit, size, y.offset=.5, 
               # clr.line="darkblue", 
               clr.line="#FC4E07",        # "#FC4E07", "red4"
               clr.marker="#FC4E07",
               vertices = TRUE)
    }
  }

  
  if (!is.null(weights) && useWilcox) { stop("Must have useWilcox as FALSE when using weights"); }
  if (!is.null(weights)) { useNested=FALSE; seCol="robust se"; FirthIfNeeded=FALSE; } else { seCol="se(coef)"; }
  if (class(y) == "Surv") { xl = "HR"; } else { xl = "OR"; }

  fbO = update(y~1, control);
  
  nms = colnames(x);
  if (is.data.frame(x) & normReal)
  { x = lapply(x, function(i) if (is.numeric(i) & !is.integer(i))
      { (i-mean(i, na.rm=TRUE))/sd(i, na.rm=TRUE); } else { i; } );
    x = do.call(data.frame, x);
  }
  else { x = data.frame(x); }
  names(nms) = names(x);
  
  if (is.null(subset)) { w = rep(TRUE, nrow(x)); } else { w = subset; }
  if (!is.logical(w)) { w2 = rep(FALSE, nrow(x)); w2[w]=TRUE; w=w2; } 
  w = w & !is.na(y);
  if (length(all.vars(control))>0)
  { w = w & !rowAnys(is.na(x[,all.vars(control),drop=FALSE]));
  }
  x = x[w,];
  y = y[w];
  weights=weights[w];
  
  vOk = setdiff(colnames(x), all.vars(control));
  pch = rep(16, length(vOk));
  
  tt = list();
  tt[[1]] = as.list(c("", nms[vOk]));
  if (genesItal)
  { if (!require(org.Hs.eg.db)) { stop("Need library org.Hs.eg.db"); }
    w = tt[[1]] %in% keys(org.Hs.eg.db, "ALIAS");
    for (i in which(w)) { tt[[1]][[i]] = parse(text=paste0("italic('",tt[[1]][[i]],"')")) }
  }
  if (parseNames)
  { tt[[1]] = lapply(tt[[1]], function(i) parse(text=i))
  }
  vls = list("", parse(text=paste0("bold(",xl,")")), expression(bold(CI)),
      expression(bolditalic(p)), expression(bold(FDR)));
  if (!fdr) { vls=vls[-length(vls)]; }
  for (i in 2:length(vls))
  { tt[[i]] = list(vls[[i]]);
  }
  
  tn = list();
  col = list();#rep("black", length(vOk));
  pp = rep(NA, length(vOk));
  for (i in seq_along(vOk))
  { fb = fbO;
    if (!identical(control, ~1))
    { if (any(is.na(x[,vOk[i]]))) # Look if some parameters got collapsed
      { av = all.vars(fb)[-1];
        Noks = sapply(av, function(id) length(unique(x[[id]][!is.na(x[,vOk[i]])])))<2;
        for (id in names(which(Noks))) { fb = update(fb, formula(paste("~ . -", id))); }
      }
    }
    if (useNested)
    { if (class(y) == "Surv")
      { a0 = coxph(fb, data=x, subset=!is.na(x[,vOk[i]]), weights=weights); }
      else { a0 = glm(fb, data=x, family=binomial, subset=!is.na(x[,vOk[i]]), weights=weights); }
    }
    
    fn = update(fb, formula(paste("~", vOk[i],"+.")));
    if (class(y) == "Surv")
    { tmp = withWarnings(coxph(fn, data=x, weights=weights)); a = tmp$value; warn = tmp$warnings;
      if (identical(control, ~1))
      { if (useWilcox && (is.factor(x[,vOk[i]]) || is.logical(x[,vOk[i]])))
        { fm = survdiff(fn, data=x)
          p = pchisq(fm$chisq, length(fm$n) - 1, lower.tail = FALSE)
        }
        else { p = summary(a)$logtest["pvalue"]; }
      }
      else
      { if (useNested) { an = anova(a0, a); p = an[["Pr(>|Chi|)"]][2]; }
        else { p=coef(summary(a)); p=p[2,"Pr(>|z|)"]; }
      }
    }
    else
    { a = glm(fn, data=x, family=binomial, weights=weights);
      if (is.na(coef(a)[2])) {  tn[[i]] = c(NA, NA, NA); pp[i] =NA; tt[[2]][[i+1]]=tt[[3]][[i+1]]=tt[[4]][[i+1]]=""; next; } 
      if (identical(control, ~1) & useWilcox)
      { if (is.numeric(x[,vOk[i]])) { p = wilcox.test(x[unclass(factor(y))==1,vOk[i]], x[unclass(factor(y))==2,vOk[i]])$p.value; }
        else { p = fisher.test(x[,vOk[i]], y)$p.value; }
      }
      else
      { if (useNested) {an = anova(a0, a); p = pchisq(an$Deviance[2], an$Df[2], lower.tail=FALSE); }
        else { p = coef(summary(a)); p=p[2,"Pr(>|z|)"]; }
      }
    }
    
    if (class(y) == "Surv")
    { if (!FirthIfNeeded || is.null(warn))
      { co = coef(summary(a)); 
        me = co[1,"coef"]; se = co[1,seCol]; cint=exp(me + se*c(-1.96, 1.96));
      }
      else
      { if (length(warn)>1 || !grepl("Loglik converged before variable", warn[[1]])) { warning("Strange warning...", unlist(warn)); }
        if (!require(coxphf, quietly=TRUE)) { stop("Need the coxphf library to be able to use Firth"); }
        x2 = x; x2$y=y; x2 = x2[!is.na(x2[,vOk[i]]),];
        co = coxphf(fn, data=x2);
        me = co$coefficients[1];
        cint = c(co$ci.lower[1], co$ci.upper[1]);
      }
    }
    else
    { co = coef(summary(a));
      me = co[2,"Estimate"];
      if (FirthIfNeeded && (exp(me)<.01*clip[1] || exp(me)>100*clip[2]))
      { if (!require(logistf, quietly=TRUE)) { stop("Need the logistf library to be able to use Firth"); }
        a = logistf(fn, data=x);
        me = a$coefficients[2];
        cint = exp(c(a$ci.lower[2], a$ci.upper[2]))
        if (!useNested) { p = a$prob[2]; }
      }
      else { se = co[2,"Std. Error"]; cint=exp(me + se*c(-1.96, 1.96)); }
    }
    if (useWilcox && !class(y) == "Surv" && !is.numeric(x[,vOk[i]]) && length(setdiff(unique(x[,vOk[i]]), NA))==2)
    { fm = fisher.test(x[,vOk[i]], y);
      me = log(fm$estimate);
      cint = fm$conf.int
    }
    if (!rmNAp || !is.na(me))
    { if (capHR && (me<log(clip[1]) || me>log(clip[2])))
      { if (exp(me)<clip[1]*1e-2) { me=-Inf; }; if (exp(me)>clip[2]*1e2) { me=+Inf;  } 
        tt[[2]][[i+1]] = f(exp(me));
        if (me<0) { pch[i]=-9668; me = log(clip[1]); }
        else { pch[i]=-9658; me = log(clip[2]); }
      } 
      else { tt[[2]][[i+1]] = f(exp(me)); }
      tt[[3]][[i+1]] = parse(text=paste("paste(",f(cint[1]), ",' to ',", f(cint[2]), ")"));
      tn[[i]] = c(me, log(cint));
      tt[[4]][[i+1]] = f(p);
      pp[i] = p;
    } else { tn[[i]] = c(NA, NA, NA); pp[i] =NA; tt[[2]][[i+1]]=tt[[3]][[i+1]]=tt[[4]][[i+1]]=""; }
  }
  tn = do.call(rbind, tn);
  colnames(tn) = c("mean", "lower", "upper");
  tn = exp(tn);
  
  tn = rbind(NA, tn);
  
  if (fdr)
  { fdr = p.adjust(pp, method='fdr');
    tt[[5]][2:(length(fdr)+1)] = lapply(fdr,function(i) if (rmNAp & is.na(i)) { "" } else f(i));
    names(tt)[5]="FDR";
  }
  
  if (!is.null(sig))
  { if (sigOnFDR) { psig = fdr; } else { psig = pp; }
    if (version<3) { col = lapply(psig, function(i) if (!is.na(i) & i < sig) { plBlue } else { plYellow }); }
    else { col = sapply(psig, function(i) if (!is.na(i) & i < sig) { "#00AFBB" } else { "slategray4" }); }
  }
  
  if (dispN)
  { N = colSums(!is.na(x[,vOk]));
    tt = c(tt[1], list(c(list(expression(bold(N))), as.list(N))), tt[2:length(tt)]);
  }
  
  is = c(FALSE, rep(FALSE, ncol(x)));
  col = c(col[[1]], col);

  wn = which(tn[,3]==tn[,2]);
  tn[wn,] = NA;
  
  names(tt)[1:4] = c("name", "OR", "CI", "P");
  tt = tt[names(tt) %in% c("name", columns) ];
  
  if (version<3)
  { if (is.null(ticks))
    { ticks = c(.1, .5, 1, 2, 10)
      if (max(exp(abs(tn)), na.rm=TRUE) < 5) { ticks = c(.2, .5, 1, 2, 5); }
      if (max(exp(abs(tn)), na.rm=TRUE) < 2) { ticks = c(.5, .75, 1, 1.5, 2); }
    } 
    attr(ticks, "labels")=as.character(ticks); ticks=log(ticks); tn=log(tn);

    grob<-forestplot(tt, tn, new_page = new_page, is.summary=is, xticks=ticks,
      fn.ci_norm=col, xlab=xl, clip=log(clip), boxsize=boxsize,
      txt_gp=fpTxtGp(xlab=gpar(cex=.7), ticks=gpar(cex=.7)), graphwidth=graphwidth, ...)
    if (plot) { print(grob); }
    return(invisible(list(tt=tt, tn=tn, p=pp, grob=grob)));
  }
  else
  { nm = lapply(tt, function(i) i[[1]]); tt2 = lapply(tt, function(i) i[-1]);
    if (plot) { basicForest(tt2, tn[-1,c(2,1,3)], xlog=TRUE, titles=nm, xlim=clip, xlab=xl, col=col[-1], pch=pch, ...) }
    return(invisible(list(tt=tt, tn=tn, p=pp)))
  }
}

######################################################################
## Scale...
plotScale = function(cols, v, atV, posx, posy, width=1, height=10, title=NULL, title.font=2,
  horizontal=FALSE)
{ fx = function(x0, dx)
  { if (par('xlog')) { exp(log(x0)+dx*strwidth("M")); } else { x0+dx*strwidth("M"); } }
  fy = function(y0, dy)
  { if (par('ylog')) { exp(log(y0)+dy*strheight("M")); } else { y0+dy*strheight("M"); } }
  
  if (horizontal) { legend_image = as.raster(matrix(cols, nrow=1)) }
  else { legend_image = as.raster(matrix(rev(cols), ncol=1)) }
  rasterImage(legend_image, posx, posy, fx(posx, width), fy(posy, height), xpd=NA)
  rect(posx, posy, fx(posx, width), fy(posy, height), border='black', xpd=NA);
  
  if (!is.null(title)) { text(posx, fy(posy, height+1.5-horizontal), title, adj=c(0,0), font=title.font, xpd=NA); }
  
  if (horizontal)
  { lines(as.vector(rbind(fx(posx, atV*width), fx(posx, atV*width), NA)),
      rep(c(fy(posy, c(0, -.5)), NA), length(atV)), xpd=NA)
    text(fx(posx, atV*width), fy(posy, -1.5), labels=v, xpd=NA, adj=c(.5,0)); 
  }
  else
  { lines(rep(c(fx(posx, width+c(0, .5)), NA), length(atV)),
      as.vector(rbind(fy(posy, atV*height), fy(posy, atV*height), NA)), xpd=NA)
    text(fx(posx, width+1), fy(posy, atV*height), labels=v, xpd=NA, adj=c(0,.5)); 
  }
}

logAxis = function(from, to, side=1, atLog=TRUE, granularity=0, inverted=FALSE, lim=1e-20)
{ # Find range
  r = from:to;
  if (granularity==1) { r = as.vector(rbind(r+log10(.5), r)); r=r[-1]; }
  if (granularity==2) { r = as.vector(rbind(r+log10(.2), r+log10(.5), r)); r=r[-1]; }
  
  #a = axTicks(side);
  #if (length(a)<length(r)) { r = r[seq(from=1, to=length(r), length.out=length(a))]; }
  
  if (inverted) { leg = sapply(10^-r, formatNice, nDigits=1, lim=lim) }
  else { leg = sapply(10^r, formatNice, nDigits=1, lim=lim) }
  
  if (!atLog) { r = 10^r; }
  axis(side, at=r, labels=leg);
  return(invisible(list(at=r, labels=leg)));
}

addAnnotArrows = function(side=1, xlab="OR", annotDir=c("Better", "Worse"), le=.05, cex=1, inset=.5)
{ if (side!=1) { stop("Does not work on y-axis (yet)"); }
  h = par('usr')[3]-par("cex.axis")*par("cex")*par("cxy")[2]*(inset+par("mgp")[1]);
  v = strwidth(paste0(xlab, "   "), cex=par("cex.axis")*par("cex"));
  arrows(m1<-mean(par('usr')[1:2])-v, h, par('usr')[1], h, xpd=TRUE, length=le)
  arrows(m2<-mean(par('usr')[1:2])+v, h, par('usr')[2], h, xpd=TRUE, length=le)
  mtext(annotDir, side=1, line=par("mgp")[1]+inset, at=(par('usr')[1:2]+c(m1,m2))/2, cex=cex);
}

#######################################################################
## Dot plot
dotPlot = function(x, cl, efi=NULL, lbls=NULL, toDisp=NULL, maxRange=4, oma=c(0,0,0,5), inMa=.5,
  col.lbl=NULL, maxP=.1, horizontal=FALSE, pch=c("full", "left", "right"), add=FALSE,
  family="Arial Unicode MS", blackBorder=TRUE, at=NULL, axPos=1, srt=45, pt.lwd=.5, cex.pch=1, ...)
{ pch= match.arg(pch); 
  if (pch == "full") { blackBorder=FALSE; }
  pch = c(full=21, left=-0x25D6L, right=-0x25D7L)[pch];
  if (is.null(lbls)) { lbls = colnames(x); }
  
  calcPetc = function(x, ws)
  { if (length(unique(ws))==1) { return(list(p=rep(NA, ncol(x)), efi=rep(NA, ncol(x)))) }
    p = apply(x, 2, function(i)
    { if (is.numeric(i)) { wilcox.test(i ~ ws)$p.value }
      else
      { if (any(rowSums(table(i,ws))==0)) { return(NA); }
        fisher.test(i, r)$p.value
      }
    });
    efi =  apply(x, 2, function(i)
    { if (is.numeric(i)) { i = scale(i); }
      r = as.numeric(try(coef(glm(ws ~ i, family=binomial))["i"], silent=TRUE))
      
      if (abs(r) > log(2*maxRange))
      { if (!require(logistf, quietly=TRUE)) { stop("Need the logistf library to be able to use Firth"); }
        x = data.frame(ws, i);
        a = logistf(ws~i, data=x);
        r = a$coefficients[2];
      }
      return(r);
    }); efi[is.na(efi)]=NA;
    #if (!is.matrix(efi)) { efi = matrix(efi, ncol=1); }
    return(list(p=p, efi=efi));  
  }
  
  if (is.character(cl)) { cl=factor(cl); }
  if (is.factor(cl))# && length(levels(cl))>2)
  { cl2 = lapply(levels(cl), function(i) list(which(cl==i), which(cl!=i)));
    names(cl2) = levels(cl); cl=cl2;
  }

  if (!is.null(oma)) { par(oma=oma); }
  if (is.null(efi))
  { if (is.list(cl))
    { pif = lapply(cl, function(ws)
      { y = rep(NA, nrow(x)); y[ws[[1]]]=TRUE; y[ws[[2]]]=FALSE; y=factor(y);
        calcPetc(x, y)
      }); names(pif) = names(cl);
    
      pi = do.call(cbind, lapply(pif, function(i) i$p))
      efi = do.call(cbind, lapply(pif, function(i) i$efi))
    }
    else
    { if (!is.factor(cl)) { cl = factor(cl); }
      tmp = calcPetc(x, cl);
      pi = tmp$p; efi = tmp$efi;
    }
    p = pi; p[] = p.adjust(pi, method='fdr');
    if (!is.null(toDisp)) { efi[!toDisp] = NA; }
    else { efi[p>maxP]=NA }
  } else { pi=NULL; }
  efi[efi>log(maxRange)]=log(maxRange); efi[efi< -log(maxRange)] = -log(maxRange);
  efi = efi/log(maxRange); # Map to -1/1
  if (!add) { plot.new(); }
  #par(mar=c(3,7,.5,.5), mgp=c(1.5,.5,0));
  if (pch==21) { fill=ifelse(efi>0, "green3", "red"); col='black'; } else { fill=NULL; col=ifelse(efi>0, "green3", "red"); }
  
  if (is.null(at)) { at = 1:nrow(efi); }
  
  if (!horizontal)
  { if (!add)
    { plot.window(xlim=c(1-inMa, ncol(efi)+inMa), ylim=range(at)+c(-inMa, inMa), xaxs="i", yaxs="i");
      if (srt==0) { axis(1, labels=colnames(efi), at=1:ncol(efi), gap.axis=0) }
      else
      { axis(1, labels=FALSE, at=1:ncol(efi), gap.axis=0);
        text(x=1:ncol(efi), y=par('usr')[1]-strheight("M"), labels=colnames(efi), srt=srt, adj=1, xpd=TRUE);
      }
      if (is.null(col.lbl)) { axis(2, labels=lbls, at=at, las=2) }
      else
      { Map(axis, side=2, at=at, col.axis=col.lbl, labels=lbls, lwd=0, las=2)
        axis(2,at=at,labels=FALSE)
      }
      if (!is.null(oma))
      { legend(x=par('usr')[2], y=par('usr')[4], bty='n', xpd=NA, legend=c("≤.25", ".5", ".66",
        "1.5", "2", "≥4"), col='black', pt.bg=rep(c("red", "green3"), each=3), title="Effect size",
        pt.cex=4*abs(log(c(.25, .5, 2/3, 3/2, 2, 4)))/log(4), pch=21, pt.lwd=.25, y.intersp=1.5);
      }
    }
    if (blackBorder)
    { points(col(efi), at[row(efi)], col='black', bg=fill, pch=pch, cex=abs(efi)*4*cex.pch, family=family, lwd=.2);
      points(col(efi), at[row(efi)], col=col, bg=fill, pch=pch,
        cex=cex.pch*(abs(efi)*4-pt.lwd/2), family=family, lwd=.2);
      
      z = strheight("M", cex=par('cex')*cex.pch);
      xs = rbind(as.vector(col(efi)), as.vector(col(efi)), NA);
      ys = rbind(as.vector(at[row(efi)]-abs(efi)*4*z*.56), as.vector(at[row(efi)]+abs(efi)*4*z*.56), NA);
      lines(as.vector(xs), as.vector(ys), col='black', lwd=pt.lwd)
    }
    else { points(col(efi), at[row(efi)], col=col, bg=fill, pch=pch, cex=cex.pch*abs(efi)*4, family=family, lwd=pt.lwd); }
  }
  else
  { if (!add)
    { plot.window(ylim=c(1-inMa, ncol(efi)+inMa), xlim=range(at)+c(-inMa, inMa), xaxs="i", yaxs="i");
      axis(2, labels=colnames(efi), at=1:ncol(efi), las=2, gap.axis=0)
     # if (is.null(col.lbl)) { axis(1, labels=rownames(efi), at=1:nrow(efi), las=2) }
      if (axPos==1) { axis(1,at=at,labels=FALSE); cpy = par()$usr[3]-par()$cxy[2]*.5; adj=1; }
      else { axis(3,at=at,labels=FALSE); cpy = par()$usr[4]+par()$cxy[2]*.5; adj=0; }
      text(x=at, y=cpy, labels=lbls, srt=srt, adj=adj, xpd=NA, col=col.lbl);
      if (!is.null(oma))
      { legend(x=par('usr')[2], y=par('usr')[4], bty='n', xpd=NA, legend=c("≤.25", ".5", ".66",
        "1.5", "2", "≥4"), col='black', pt.bg=rep(c("red", "green3"), each=3), title="Effect size",
        pt.cex=4*abs(log(c(.25, .5, 2/3, 3/2, 2, 4)))/log(4), pch=21, pt.lwd=pt.lwd, y.intersp=1.5);
      }
    }
    #points(at[row(efi)], col(efi), col=col, bg=fill, pch=pch, cex=abs(efi)*4, family=family, lwd=.2);
    if (blackBorder)
    { points(at[row(efi)], col(efi), col='black', bg=fill, pch=pch, cex=abs(efi)*4, family=family, lwd=pt.lwd);
      points(at[row(efi)], col(efi), col=col, bg=fill, pch=pch, cex=abs(efi)*4-pt.lwd, family=family, lwd=pt.lwd);
    }
    else { points(at[row(efi)], col(efi), col=col, bg=fill, pch=pch, cex=abs(efi)*4, family=family, lwd=pt.lwd); }
  }    
  return(invisible(list(pi=pi, efi=efi)));
}

legendDotPlot = function(x, y, double=FALSE, doubleAnnot, interline=1,
  family="Arial Unicode MS", pt.lwd=.5, horizontal=FALSE, cex.pch=1)
{ li = strheight("M", cex=par('cex'))*1.5*interline;
  em = strwidth("M", cex=par('cex'));
  
  par(xpd=NA);
  
  if (horizontal) {text(x+7.5*em,y,"Effect size", adj=c(.5,1), font=2); }
  else { text(x,y,"Effect size", adj=c(.5,1)); }
  y2 = y;
  
  if (horizontal)
  { xs = x + 3*(0:5)*em; dxs=0;
    ys = rep(y2-2*li, 6); dys=-2*li; adj = c(.5,0);
  }
  else { xs = rep(x-em, 6); dxs = 2*em; ys = y2 - (0:5)*li*1.5-2.5*li; dys = 0; adj=c(0,.5) }
  
  v = 4*abs(log(c(.25, .5, 2/3, 3/2, 2, 4)))/log(4)
  points(xs, ys, pch=21, lwd=pt.lwd, cex=v*cex.pch, xpd=NA,
     bg = rep(c("red", "green3"), col='black', each=3))
  text(xs + dxs, ys+dys, c("≤.25", ".5", ".66", "1.5", "2", "≥4"), adj=adj, xpd=NA)
  if (double) # OK bug if horizontal
  { z = strheight("M", cex=par('cex'));
    for (i in 1:6)
    { lines(c(xs[i], xs[i]), ys[i]+c(-z, z)*v[i]*.6, col='black', lwd=pt.lwd, xpd=NA) }
    y3 = ys[6]-1.5*li;
    text(x-em+c(-em, em)/2, rep(y3, 2), doubleAnnot, pos=c(2,4), xpd=NA)
    text(x-em, y3, "|", xpd=NA)
    #arrows(x-em+3*c(em, -em), rep(y3-3*li, 2), x-em+.3*c(em, -em), c(y3, y3), lwd=1, length=.1)
    #text(x-em-3*c(em, -em), rep(y3-3*li, 2), doubleAnnot, adj=c(1,.5), srt=90)
  }
}

#######################################################################
## Box plot functions
myBoxplot = function(x, f, subset=NULL, subset2=NULL, baseLab=NULL, index=NULL, ylim, subLine=NA, subSide=1, subAt=NA,
  colPoints='blue', Plog=NULL, alternative="two.sided", parseLeg=FALSE, leg45=FALSE, pch=1, normGroup=NULL, compToAll=FALSE,
  compToConstant=NULL, isComp=FALSE, cexSignif=1, alpha=1, starsOnFDR=FALSE,
  xlab="", ylab="", pvalLim=1e-20, xaxt, sepBaselab=" - ", boldBaselab=FALSE, paramTest=FALSE, pairsToShow=NULL,
  log=ifelse(is.null(Plog) || is.na(Plog), "", "y"), addN=FALSE, main=NULL, colXlbl='black', ...)
# x: matrix
# f: grouping factor; if missing compare columns of x or items in x (if list)
# subset: if only using a subset
# subset2: remove from calc, but plot
# baseLab: annotation
# subLine: where the pvalue is shown
# Plog: if !missing, adds a value to everything and use log axis with the right rescaling
{ paired = FALSE;
  if (!is.null(subset) && !is.null(subset2)) { subset2 = subset2[subset]; }
  if (missing(f))
  { if (length(dim(x))==2)
    { if (is.null(colnames(x))) { colnames(x) = 1:ncol(x); }
      if (length(colPoints)==1) { colPoints = rep(colPoints, nrow(x)); }
      if (length(colPoints) == nrow(x)) { colPoints=rep(colPoints, ncol(x)); }
      if (!is.null(subset)) { x = x[subset,]; colPoints=colPoints[subset]; }
      xb = x; if (missing(isComp)) { isComp = TRUE; }
      f = factor(rep(colnames(x), each=nrow(x)), levels=colnames(x));
      x = as.vector(x);
      if (ncol(xb)==2 & isComp) { paired=TRUE; }
      if (length(subset2) == length(f)/2) { subset2 = rep(subset2, 2); }
    }
    if (is.list(x))
    { if (is.null(names(x))) { names(x) = seq_along(x); }
      f = factor(rep(names(x), sapply(x, length)), levels=names(x));
      x = xb = unlist(x);
    }
    if (length(colPoints)==1) { colPoints = rep(colPoints, length(f)); }
  }
  else
  { if (length(colPoints)==1) { colPoints = rep(colPoints, length(f)); }
  if (!is.null(subset)) { x=x[subset]; f=factor(f[subset]); subset2=subset2[subset]; colPoints=colPoints[subset]; }
    xb = x;
  }
  if (is.null(subset2)) { subset2 = rep(TRUE, length(x)); }
  if (alpha != 1) { require(scales); colPoints = alpha(colPoints, alpha); }

  if (!is.factor(f)) { f = factor(f); }
  if (addN) { t=table(f[!is.na(x)]); levels(f) = paste0(levels(f), " (N=", t[levels(f)], ")"); } 
  
  if (length(levels(f))==2)
  { if (paired)
    { if (paramTest) { p = try(t.test(xb[,2]-xb[,1], alternative=alternative, subset=subset2)$p.value, silent=TRUE); }
      else { p = try(wilcox.test(xb[,2]-xb[,1], alternative=alternative, subset=subset2)$p.value, silent=TRUE); }
    }
    else
    { if (paramTest) { p = try(t.test(x~f, alternative=alternative, subset=subset2, var.equal=TRUE)$p.value, silent=TRUE); }
      else { p = try(wilcox.test(x~f, alternative=alternative, subset=subset2)$p.value, silent=TRUE); }
    }
  }
  else
  { if (is.ordered(f)) { p = cor.test(x[subset2],unclass(f[subset2]),method='k',
        alternative=alternative, na.rm=TRUE)$p.value; }
    else
    { if (paramTest) { p = summary(aov(x ~ f, subset=subset2))[[1]][1,"Pr(>F)"] }
      else { p = try(kruskal.test(x[subset2], f[subset2])$p.value, silent=TRUE); }
    }
  }
  if (is(p, "try-error")) { p=NA; }

  w = !is.na(x) & !is.na(f);
  x = x[w]; f=f[w]; colPoints = colPoints[w]; subset2 = subset2[w];
  
  if (!(is.null(Plog) || is.na(Plog))) { x = x+Plog; xb = xb+Plog; }
  if (missing(ylim)) { ylim=range(x); }
  
  if (!is.null(pairsToShow))
  { pp = apply(pairsToShow, 1, function(i) wilcox.test(x[f==i[1]], x[f==i[2]])$p.value)
    pairsToShow = pairsToShow[pp<.05,,drop=FALSE]; pp = pp[pp<.05]
    if (nrow(pairsToShow)!=0) 
    { # Find how to pile those up
      ob = matrix(match(pairsToShow, levels(f)), ncol=2);
      ob = cbind(pmin(ob[,1], ob[,2]), pmax(ob[,1], ob[,2]))
      o = order(ob[,1], ob[,2]-ob[,1]);
      pairsToShow = pairsToShow[o,,drop=FALSE]; ob = ob[o,,drop=FALSE]; pp=pp[o];
      lvl = rep(NA, nrow(pairsToShow)); mFill=rep(-Inf, nrow(pairsToShow))
      for (i in sort(unique(ob[,1])))
      { fst = 1;
        for (j in which(ob[,1]==i))
        { while(TRUE) { if (mFill[fst]<=ob[j,1]) { mFill[fst] = ob[j,2]+((ob[j,2]-ob[j,1])==1); lvl[j] = fst; break; }; fst=fst+1; }
        }
      }
      sz = strheight("I", unit="figure");
      if (log=="") { maxy = ylim[2]; ylim[2] = ylim[2] + (ylim[2]-ylim[1]) * max(lvl)*sz*6 }
    }
    else { pairsToShow=NULL; }
  }
  
  plot.new(); plot.window(xlim=c(.5, length(levels(f))+.5), ylim=ylim, log=log);
  title(main=main);
  
  fl = unclass(f);
  r = list();
  for (i in unique(fl))
  { points(r[[i]]<-runif(sum(fl==i))*.4+i-.2, x[fl==i], col=colPoints[fl==i], pch=pch);
  }
  if (!is.null(index)) { addFigIndex(LETTERS[index]); }
  
  if (isComp)
  { cols = adjustcolor(c("red", "grey", "blue"), alpha.f=.5);
    for (i in 1:(ncol(xb)-1))
    { w1 = which(!is.na(xb[,i])); w2 = which(!is.na(xb[,i+1]));
      id1 = id2 = rep(NA, nrow(xb)); id1[w1] = seq_along(w1); id2[w2] = seq_along(w2);
      for (j in 1:nrow(xb))
      { if (any(is.na(xb[j,i:(i+1)]))) { next; }
        lines(c(r[[i]][id1[j]], r[[i+1]][id2[j]]), xb[j,i:(i+1)], col=cols[2+sign(xb[j,i]-xb[j,i+1])])
      }
    }
  }
  
  if (!is.null(normGroup))
  { nor = which(levels(f)==normGroup);
    for (i in setdiff(seq_along(levels(f)), nor))
    { p = wilcox.test(x[unclass(f)==i], x[unclass(f)==nor])$p.value;
      if (p<.1)
      { Signif = symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 
          0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "));
        mtext(paste0("(", Signif, ")"), side = 3, line = -par("mgp")[1], outer = FALSE, at = i, cex=cexSignif)
      }
    }
  }
  
  if (compToAll)
  { for (i in seq_along(levels(f)))
    { p2 = wilcox.test(x[unclass(f)==i], x[unclass(f)!=i])$p.value;
      if (p2<.1)
      { Signif = symnum(p2, corr = FALSE, na = FALSE, cutpoints = c(0, 
          0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "));
        mtext(paste0("(", Signif, ")"), side = 3, line = -par("mgp")[1], outer = FALSE, at = i, cex=cexSignif)
      }
    }
  }
  
  if (!is.null(compToConstant))
  { if (log=='y') { z = log(x); compToConstant=log(compToConstant) } else { z = x; }
    p2 = sapply(levels(f), function(i) { wilcox.test(z[f==i], mu=compToConstant)$p.value; } )
    if (starsOnFDR) { p2 = p.adjust(p2, method='fdr'); }
    for (i in seq_along(levels(f)))
    { if (p2[i]<.1)
      { Signif = symnum(p2[i], corr = FALSE, na = FALSE, cutpoints = c(0, 
          0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "));
        mtext(paste0("(", Signif, ")"), side = 3, line = -par("mgp")[1], outer = FALSE, at = i, cex=cexSignif)
      }
    }
  }
  
  if (!is.null(pairsToShow))
  { sz = strheight("I"); lvl=lvl-1;
    for (i in 1:nrow(pairsToShow))
    { #browser();
      lines(c(ob[i,1], ob[i,1], ob[i,2], ob[i,2]), maxy + 2.5*sz*lvl[i] + (c(0, sz/2, sz/2, 0)))
      text((ob[i,1]+ob[i,2])/2, maxy+sz*(lvl[i]*2.5-.25), formatNiceP(pp[i]), pos=3)
    }
  }
  
    
  if (missing(xaxt)) { xaxt = ifelse(leg45, "n", "s"); }
  if (parseLeg) { nms = parse(text=levels(f)) } else { nms = levels(f); }
  if (is.null(Plog) || is.na(Plog))
  { a = boxplot(x~f, outline=FALSE, ylim=ylim, xaxt=xaxt, subset=subset2, col=NA, log=log, add=TRUE, names=nms, ...); }
  else
  { a = boxplot(x~f, outline=FALSE, ylim=ylim, log=log, yaxt="n", subset=subset2, xaxt=xaxt, col=NA, add=TRUE, names=nms, ...);
    axisP(Plog);
  }
  title(xlab=xlab, ylab=ylab);
  if (leg45)
  { axis(1, at = seq_along(a$name), labels=FALSE);
    if (par("ylog")) { cpy = 10^(sum(par("usr")[3:4] * c(1.05, -.05))) }
    else { cpy = par()$usr[3]-par()$cxy[2]*.5; }
    text(x=seq_along(a$names), y=cpy, labels=nms, srt=45, adj=1, xpd=TRUE, col=colXlbl);
  }
  
  if (!is.null(subSide))
  { mtext(formatNiceP(p, baseLab, lim=pvalLim, sepDescr=sepBaselab, boldDescr=boldBaselab),
      line=subLine,  side=subSide, cex=par("cex"), at=par("usr")[1]+subAt*diff(par("usr")[1:2]));
  }
  return(p);
}

combineFact = function(d, useName=TRUE)
{ if (is.null(dim(d)) || length(dim(d))==1)
  { if (is.factor(d)) { return(d) } else { return(factor(d)); } }
  d = data.frame(d);
  if (useName) { d2 = do.call(paste, lapply(names(d), function(i) paste0(i, ":", d[[i]]))) }
  else { d2 = do.call(paste, d); }
  d2[!complete.cases(d)] = NA;
  return(factor(d2))
}


plotSurv = function(su, x, pval=TRUE, col=NULL, legend.title=NULL, subset=NULL, addMedian=FALSE, medianUnit="",
  addN=FALSE, hr=FALSE, ppos=.5, pposx=par("usr")[2]*.95, legSep=" ", useName=TRUE,
  legendPos='topright', addTable=FALSE, marTable=0,
    leftTable=3, cexTable=par("cex"), test="survdiff", useCox=FALSE, parseCateg=FALSE, mark.time=TRUE,
    nperm=1e4, args.legend=NULL, moveTimeOne=0, ...)
{ y = combineFact(x, useName=useName);
  if (is.null(col)) { col = 1:length(levels(y)); }
  if (!is.null(subset)) { su=su[subset]; y = factor(y[subset, drop=FALSE]); }
  if (addTable) { marBak = par('mar'); bl = par('mar')[1]; par(mar=par("mar")+c(length(levels(y))+1+marTable, 0,0,0)) }
  plot(fmSurv<-survfit(su ~ y), col=col, mark.time=mark.time, ...)
  leg = levels(y);
  N = ""; if (addN) { N = paste0("N=", table(y), " "); }
  if (addMedian) { N = paste0(N, "med=", round(quantile(survfit(su~y), .5)$quantile, digits=1),
    medianUnit) }
  if (addN || addMedian)
  { if (parseCateg) { leg = paste0("paste(",leg, ",' (", N,  ")')")}
    else { leg = paste0(leg, legSep, "(", N, ")") }
  }
  if (parseCateg) { leg = parse(text=leg); }
  if (!is.null(legendPos))
  { do.call(legend, c(list(x=legendPos, fill=col, legend=leg,
      title=legend.title, y.intersp=ifelse(legSep=="\n", 1.6, 1)), args.legend))
  }
  if (!is.null(pval))
  { if (is.numeric(pval)) { p=pval; }
    else
    { if (is.ordered(y)) { y = unclass(y); p = getCoxP(coxph(su~y)); }
      if (is.numeric(y)) { p = getCoxP(coxph(su~y)) }
      else
      { if (missing(test)) { if (useCox) { test="cox"; } else { test="survdiff"; } }
        if (test=="cox") { p = getCoxP(coxph(su~y)); }
        if (test=="survdiff") { fm = survdiff(su~y); p=pchisq(fm$chisq, length(fm$n)-1, lower.tail=FALSE)}
        if (test=="perm") { library(clinfun); data = data.frame(su=su, y=y); p = permlogrank(su~y, data=data, nperm=nperm)$perm.p; }
      }
    }
    text(pposx, ppos, formatNiceP(p), pos=2)
    if (hr)
    { fm = coef(summary(coxph(su ~ y)));
      if (nrow(fm)!=1) { stop("No hr if more than one categ"); }
      fm = fm[1,];
      #hr = paste0("HR=", signif(exp(fm["coef"]), 2), " CI=",
      #  signif(exp(fm["coef"]-fm["se(coef)"]*1.96), 2), " to ", signif(exp(fm["coef"]+fm["se(coef)"]*1.96), 2))
      if (abs(log10(exp(fm["coef"])))>5) { hr="HR not defined"; } else
      { hr = paste0("HR=", formatNdig(exp(fm["coef"]), 2), " CI=",
        formatNdig(exp(fm["coef"]-fm["se(coef)"]*1.96), 2), " to ", formatNdig(exp(fm["coef"]+fm["se(coef)"]*1.96), 2))
      }
      text(pposx, ppos-1.5*strheight("I"),hr, pos=2)
    }
  }
  if (addTable)
  { i = axTicks(1); i[1]=i[1]+moveTimeOne;
    tbl = summary(fmSurv, times=i)
    tbl = data.frame(time=tbl$time, nrisk=tbl$n.risk, strata=tbl$strata)
    levels(tbl$strata) = sub("y=", "", levels(tbl$strata))
    mtext("Number at risk", font=2, line=bl, at = 0, side=1, adj=0, cex=cexTable)
    leg=levels(tbl$strata); if (parseCateg) { leg = parse(text=paste("bold(", leg, ")")); }
    mtext(leg, side=1, font=2, col=col, at=0-strwidth("I")*leftTable, adj=1,
      line = bl+seq_along(levels(tbl$strata)), cex=cexTable)
    for (j in seq_along(levels(tbl$strata)))
    { w = tbl$strata==levels(tbl$strata)[j];
      mtext(tbl$nrisk[w], line=bl+j, at=tbl$time[w], side=1, cex=cexTable)
    }
    par(mar=marBak)
  }
}

###############################################################
## barplot with confidence intervals
barplotConfint = function(a, b, length=.1, addN=FALSE, addP=FALSE, lineP=3, sideP=1, p, ...)
{ if (missing(b))
  { if (ncol(a)==2) { a = t(a); }
    if (nrow(a)==2) { b = colSums(a); a=a[1,]; } else { stop("not a matrix and not 2 params")}
  }
  w = barplot(100*a/b, ...);
  ci = 100*sapply(seq_along(a), function(i) { if (b[i]>0) { binom.test(a[i], b[i])$conf.int; } else { c(NA,NA); } })
  arrows(x0=w, y0=ci[1,], y1=ci[2,], angle=90, code=3, length=length)
  if (addN) { mtext(paste0("N=",b), line=2, at=w, side=1, cex=par('cex')); }
  if (addP)
  { if (missing(p)) { p=fisher.test(rbind(a,b-a))$p.value; }
    mtext(formatNiceP(p), line=lineP, side=sideP, cex=par('cex'));
    w = list(at=w, p=p);
  }
  return(invisible(w))
}

#########################################################
# Other functions
axisP = function(Plog, side=2, nice=FALSE, maxVal=NULL)
{ ticks = lbl = axTicks(side);
  if (!is.na(Plog))
  { if (ticks[1]>Plog*.7) { ticks = lbl = c(0, ticks); }
    ticks[-1] = ticks[-1]+Plog;
    lbl[1] = 0; ticks[1] = Plog;
  }
  if (nice)
  { lbl = sapply(lbl, formatNice, nDigits=1);
    sw = abs(sapply(lbl, strwidth)) # In log...
    btw = log10(ticks[1]/ticks[2]) # Space between 2 ticks
    b = sw[1]/2+btw/5+btw;
    w = rep(TRUE, length(lbl))
    for (i in 2:length(lbl))
    { if (i*btw - sw[i]/2 < b) { w[i] = FALSE; } else { b = i*btw + sw[i]/2+btw/5; }
    }
    ticks = ticks[w]; lbl = lbl[w]
  }
  if (!is.null(maxVal)) { w = ticks-Plog<maxVal; ticks=ticks[w]; lbl=lbl[w]; }
  #browser();
  axis(side, at=ticks, labels=lbl);
}

## Add fig indice
.FigIndex = new.env(); .FigIndex$index = "A";

initFigIndex = function(m="A", line, cex, delta, clear=FALSE)
{ if (clear) { .FigIndex <<- new.env(); }
  .FigIndex$index = m;
  if (!missing(line)) { .FigIndex$line=line}
  if (!missing(delta)) { .FigIndex$delta=delta}
  if (!missing(cex)) { .FigIndex$cex=cex}
}

addFigIndex = function(m, line=1, cex=1.5, delta=0.1)
{ if (missing(m)) { m = .FigIndex$index; } else { .FigIndex$index = m;  }
  if (missing(line) && !is.null(.FigIndex$line)) { line = .FigIndex$line; }
  if (missing(cex) && !is.null(.FigIndex$cex)) { cex = .FigIndex$cex; }
  if (missing(delta) && !is.null(.FigIndex$delta)) { delta = .FigIndex$delta; }
  .FigIndex$index = rawToChar(as.raw(as.numeric(charToRaw(m))+1));
  u = par("usr"); v = par("plt")
  sc = (u[2]-u[1])/(v[2]-v[1]);
  u[1] = u[1] - v[1]*sc; u[2] = u[2] + (1-v[2])*sc
  at=u[1]+(u[2]-u[1])*delta;
  if (par("xlog")) { at = 10^at; }
  mtext(eval(parse(text=paste0("expression(bold(", m, "))"))), side=3, line=line,
    at=at, cex=cex, adj=1);
}

####################################################################
## Make all subgroups
makeAllSubgroups = function(x, groups)
{ if (is.matrix(x) || is.data.frame(x))
  { ret = lapply(1:ncol(x), function(i) makeAllSubgroups(x[,i], groups));
    ret = do.call(cbind, lapply(seq_along(ret), function(i) 
    { r = ret[[i]]; colnames(r)=paste(colnames(x)[i], colnames(r));
      return(r);
    }));
    return(ret);
  }
  
  if (is.vector(groups)) { groups=data.frame(groups)}
  for (i in 1:ncol(groups)) { if (!is.factor(groups[,i])) { groups[,i] = factor(groups[,i]); } }
  
  idxM = sapply(groups, function(i) length(levels(factor(i))));
  x = matrix(x, ncol=prod(idxM+1), nrow=length(x));
  idx = rep(0, ncol(groups));
  
  k=1;
  cn = character(ncol(x));
  while (1)
  { for (i in seq_along(idx))
    { if (idx[i]>0)
      { x[unclass(groups[,i])!=idx[i] | is.na(groups[,i]), k] = NA; } 
    }
    cn[k] = paste(sapply(which(idx>0), function(i) levels(groups[[i]])[idx[i]]), collapse=" ");
    sapply(groups[idx>0], function(i) levels(i))
    k = k+1;
    for (i in seq_along(idx))
    { if (idx[i]<idxM[i]) { idx[i] = idx[i]+1; break; }
      idx[i] = 0;
    }
    if (all(idx==0)) { break; }
  }
  cn[1] = "All";
  colnames(x) = cn;
  x = x[,colSums(!is.na(x))>0];
  return(x);
}

#####################################
myAddTextLabels <- function(xCoords, yCoords, labels, cex.label=1, col.label="red", col.line="black", col.background=NULL,
                          lty=1, lwd=1, border=NA, avoidPoints=TRUE, keepLabelsInside=TRUE, cex.pt=1, lblBorder=0,
                          figBorder=0, heightPad=0.5, widthPad=0.04, axisLimits=graphics::par("usr")){

  #######################################################
  # Check that the input data are in the correct format #
  #######################################################
  
  # Are each of coordinate vectors the same length?
  if(length(xCoords) != length(yCoords)){
    stop("addTextLabels() The vectors containing the X and Y coodinates must be the same length.")
  }
  if(length(xCoords) != length(labels)){
    stop("addTextLabels() The vector of labels must be the same length as the coordinate vectors.")
  }
  
  #######################
  # Get the axis limits #
  #######################

  # Get the axis limits
  xDelt = figBorder*(axisLimits[2]-axisLimits[1]); yDelt = figBorder*(axisLimits[4]-axisLimits[3]);
  axisLimits = axisLimits + c(xDelt, -xDelt, yDelt, -yDelt); 

  ############################
  # Check for NA coordinates #
  ############################
  
  # Check if any NA coordinates present
  indicesOfNAs <- which(is.na(xCoords) | is.na(yCoords))
  if(length(indicesOfNAs) > 0){
    
    # Send warning
    warning("NA values present in coordinates provided. These are ignored.")
    
    # Check for each of the parameters that can have multiple parameters
    if(length(col.line) == length(xCoords)){
      col.line = col.line[-indicesOfNAs]
    }
    if(length(col.background) == length(xCoords)){
      col.background = col.background[-indicesOfNAs]
    }
    if(length(lty) == length(xCoords)){
      lty = lty[-indicesOfNAs]
    }
    if(length(lwd) == length(xCoords)){
      lwd = lwd[-indicesOfNAs]
    }
    if(length(border) == length(xCoords)){
      border = border[-indicesOfNAs]
    }

    # Remove the NA coordinates
    xCoords <- xCoords[-indicesOfNAs]
    yCoords <- yCoords[-indicesOfNAs]
    
    # Remove the respective labels
    labels <- labels[-indicesOfNAs]
  }
  
  ############################
  # Check if axes are logged #
  ############################

  # Check X axis
  xAxisLogged <- FALSE
  if(graphics::par("xlog")){

    # Note that X axis was logged
    xAxisLogged <- TRUE

    # Log the X coordinates
    xCoords <- log10(xCoords)

    # Reset the X axis logged flag - fools points and polygon commands below
    graphics::par(xlog=FALSE)
  }

  # Check Y axis
  yAxisLogged <- FALSE
  if(graphics::par("ylog")){

    # Note that Y axis was logged
    yAxisLogged <- TRUE

    # Log the Y coordinates
    yCoords <- log10(yCoords)

    # Reset the Y axis logged flag - fools points and polygon commands below
    graphics::par(ylog=FALSE)
  }

  ###############################
  # Store the point information #
  ###############################

  # Store the input coordinates and labels
  pointInfo <- list("X"=xCoords, "Y"=yCoords, "Labels"=labels, "N"=length(xCoords), "cex"=cex.pt)

  # Set the amount to pad onto height and width
  if(!is.null(col.background)){
    heightPad <- lblBorder
    widthPad <- lblBorder
  }

  # Calculate the label heights and widths
  wLbl = which(!is.na(labels))
  pointInfo <- basicPlotteR:::calculateLabelHeightsAndWidths(pointInfo=pointInfo, cex=cex.label,
                                              heightPad=heightPad, widthPad=widthPad)
                                            

  ###########################################
  # Produce a list of alternative locations #
  ###########################################

  # Generate the alternative locations
  alternativeLocations <- basicPlotteR:::generateAlternativeLocations(axisLimits)

  # Calculate the distance between the actual and alternative points - rescale X axis remove axis range bias
  distances <- basicPlotteR:::euclideanDistancesWithRescaledXAxis(pointInfo, alternativeLocations, axisLimits)

  ###############################################################
  # Create a list to store the information about plotted points #
  ###############################################################

  # Initialise the list to store the information about plotted labels
  plottedLabelInfo <- list("X"=c(), "Y"=c(), "Height"=c(), "Width"=c(), "N"=0)

  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################

  # Plot the point label
  for(i in wLbl){

    # Set the colours for plotting the label - allows multiple colours and cycling through colours
    labelColour <- basicPlotteR:::setOption(options=col.label, index=i)
    backgroundColour <- basicPlotteR:::setOption(options=col.background, index=i)
    borderColour <- basicPlotteR:::setOption(options=border, index=i)

    # Set the line characteristics
    lineColour <- basicPlotteR:::setOption(options=col.line, index=i)
    lineType <- basicPlotteR:::setOption(options=lty, index=i)
    lineWidth <- basicPlotteR:::setOption(options=lwd, index=i)

    # Get the information for the current point
    x <- pointInfo$X[i]
    y <- pointInfo$Y[i]
    label <- pointInfo$Labels[i]
    height <- pointInfo$Heights[i]
    width <- pointInfo$Widths[i]

    # Get a new location
    newLocationIndex <- basicPlotteR:::chooseNewLocation(pointInfo, i, alternativeLocations, distances, plottedLabelInfo, axisLimits, keepLabelsInside)

    #if (i==39) { browser(); }
    # Is the current point too close to others?
    if(alternativeLocations$N != 0 && newLocationIndex != -1 && 
       (avoidPoints == TRUE || basicPlotteR:::tooClose(x, y, height, width, plottedLabelInfo) || basicPlotteR:::outsidePlot(x, y, height, width, axisLimits))){

      # Get the coordinates for the chosen alternate location
      altX <- alternativeLocations$X[newLocationIndex]
      altY <- alternativeLocations$Y[newLocationIndex]

      # Add line back to previous location
      basicPlotteR:::addLineBackToOriginalLocation(altX=altX, altY=altY, x=x, y=y, label=label,
                                    cex=cex.label, col=lineColour, lty=lineType, lwd=lineWidth, heightPad=heightPad,
                                    widthPad=widthPad)

      # Add label
      basicPlotteR:::addLabel(x=altX, y=altY, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour, heightPad=heightPad, widthPad=widthPad)

      # Append the plotted label information
      plottedLabelInfo <- basicPlotteR:::addPlottedLabel(x=altX, y=altY, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)

      # Remove the alternative plotting location used
      alternativeLocations$X <- alternativeLocations$X[-newLocationIndex]
      alternativeLocations$Y <- alternativeLocations$Y[-newLocationIndex]
      alternativeLocations$N <- alternativeLocations$N - 1
      distances <- distances[, -newLocationIndex]

    }else{

      # Add label
      basicPlotteR:::addLabel(x=x, y=y, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour,
               heightPad=heightPad, widthPad=widthPad)

      # Append the plotted label information
      plottedLabelInfo <- basicPlotteR:::addPlottedLabel(x=x, y=y, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)
    }
  }

  #####################################################################################
  # Return axes logged flags to original state - for if person makes any future plots #
  #####################################################################################

  graphics::par(xlog=xAxisLogged)
  graphics::par(ylog=yAxisLogged)

}