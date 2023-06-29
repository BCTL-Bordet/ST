figDir = "~/Documents/BCTL/ST/TNBC/figsTest/";

# Figure files
#
# Fig 1d - morphoDistUp.pdf
# Fig 2b - alluvialDec.2.pdf
# Fig 2c - TvsS.bar.dec.Tumor.Sigs.horiz.pdf
# Fig 2d - TvsS.bar.dec.Stroma.Sigs.horiz.pdf
# Fig 2f - TvsS.Mbp.hm.pdf
# Fig 3b - tlsViolin.pdf
# Fig 3c - tlsGene.pdf
# Fig 3e/f/g - tlsAllCliSmall.pdf
# Fig 4c - barByClust.pdf
# Fig 4d - tblsMC2.pdf
# Fig 4e - forestMC.PB.pdf
# Fig 4f - forestMC.All3.pdf
# Fig 5a - heatmapMC.pdf
# Fig 5b - barplotMCvsET.pdf
# Fig 5c - barEcot.pdf
# Fig 5d - EcoSingle.horiz.pdf
# Fig 5e - forestET.All3.pdf
# Fig 5f - survEcotypeNew2.pdf
#
# Supplementary figures
# Fig S2a - aucs.pdf
# Fig S3a - TvsS.bar.dec.Tumor.*.horiz.pdf
# Fig S3b - TvsS.bar.dec.Stroma.*.horiz.pdf
# Fig S4a/b/c - survTSdec.final.2.pdf
# Fig S4d - TvsS.bar.dec.Mpb.Cli.horiz.pdf, TvsS.bar.dec.Mpb.Annot.horiz.pdf
# Fig S5a - tlsXcl.pdf
# Fig S5b - sigCharoentongInTLS.FP.pdf
# Fig S5c - sigCharoentongInTLS.pdf
# Fig S5d - tlsGO.VSlympho.pdf
# Fig S5e - tlsGO.VSrest.pdf
# Fig S6a/b - aucTLS.pdf
# Fig S6c - tlsAllCliRest.pdf
# Fig S6d - tlsIspy.pdf
# Fig S6e - fpTLSimmuno.pdf
# Fig S6f/g/h - tlsAllCliRest.pdf
# Fig S9a - meanAUC.pdf
# Fig S9b/c - MCcomp.Sigs.pdf
# Fig S10 - MCcomp.Genes.pdf
# Fig S11 - MCcomp.xCell.pdf
# Fig S12 - forestMC.PB.pdf
# Fig S13 - forestMC.Bulk.pdf
# Fig S14 - forestMC.MB.pdf
# Fig S15 - forestMC.SCANB.pdf
# Fig S16 - forestMC.All3.pdf
# Fig S17a - Ncluster.pdf
# Fig S17b - barEcot.pdf
# Fig S17c - EcoComp.Cli.pdf
# Fig S17d - EcoComp.Annot.pdf
# Fig S17e/f - ET.HM.xc.pdf
# Fig S18a/b - ET.HM.sig.pdf
# Fig S19a/b - ET.HM.genes.pdf
# Fig S20a - confMatEcot.pdf
# Fig S20b - EcoComp.Sigs.Bulk.pdf
# Fig S20c - EcoComp.Sigs.MB.pdf
# Fig S21 - EcoComp.Sigs.ScanB.pdf
# Fig S22a - EcoComp.xCell.Bulk.pdf
# Fig S22b - EcoComp.xCell.MB.pdf
# Fig S22c - EcoComp.Sigs.ScanB.pdf
# Fig S23a - EcoComp.Genes.Bulk.pdf
# Fig S23b - EcoComp.Genes.MB.pdf
# Fig S24 - EcoComp.Genes.ScanB.pdf
# Fig S25 - survEco1vAll.DRFS.pdf / survEco1vAll.Bulk.DRFS.pdf / survEco1vAll.MB.DRFS.pdf / survEco1vAll.SCANB.DRFS.pdf
# Fig S26 - survEco1vAll.iBCFS.pdf / survEco1vAll.Bulk.iBCFS.pdf / survEco1vAll.MB.EFS.pdf / survEco1vAll.SCANB.iBCF.pdf
# Fig S27 - survEco1vAll.OS.pdf / survEco1vAll.Bulk.OS.pdf / survEco1vAll.MB.OS.pdf / survEco1vAll.SCANB.OS.pdf
# Fig S28 - survEco1vAll.RFS.pdf / survEco1vAll.Bulk.RFS.pdf
# Fig S29a/b - forestET.PB.pdf
# Fig S30a/b - forestET.Bulk.pdf
# Fig S31a/b - forestET.MB.pdf
# Fig S32a/b - forestET.SCANB.pdf
# FIg S33a/b - forestET.All3.pdf

# Replace 4d, S6e, S6a/b, S2a

# Morpho etc
##############
n2 = cli$annotClean[,c("Tumor", "Stroma", "Acellular stroma", "Lymphocyte", 
  "Fat tissue", "Necrosis", "in situ", "Lymphoid nodule", "Vessels", "Lactiferous duct",
  "Nerve", "Heterologous elements")];
co = colAnn2[colnames(n2)]; co[2] = colAnn2["Stroma cell"];
o = do.call(order, data.frame(cli$barPB, n2));

n3 = do.call(rbind, tapply(1:nrow(n2), cli$barPB, function(i) colMeans(n2[i,])))
n3 = n3[5:1,];
pdf(paste0(figDir, "morphoDistUp.pdf"), height=4, width=16)
par(mfcol=c(1,2), mar=c(3,3,1,.5), mgp=c(1.5,.5,0), oma=c(0,0,0,12));
for (j in 1:2)
{ if (j==1) { yl = c(0,100); } else { yl = c(95, 100); }
  h = barplot(t(n3), col=co, ylab="% pixels", border=NA, ylim=yl, xpd=FALSE)
  par(xpd=TRUE)
  if (j==1)
  { lines(c(h[length(h)]+.5, par('usr')[2]+2.5*strwidth("M")), c(100, 100), xpd=NA, lwd=2)
    lines(c(h[length(h)]+.5, par('usr')[2]+2.5*strwidth("M")), c(95, 0), xpd=NA, lwd=2)
  }
 
  if (j==2) { legend('topright', inset=c(-max(strwidth(colnames(n2), units="figure"))*1.3,0),
    legend=rev(colnames(n2)), fill=rev(co), xpd=NA, bty='n')} 
}
dev.off();

#####################
## Tumor / stroma PB
#####################
load(paste0(dataDir, 'classification/loosReg.RData'));
a = c("Tumor", "Stroma", "Necrosis", "Fat tissue", "Vessels", "in situ", "Lymphoid nodule",
  "Lymphocyte", "Lactiferous duct");
pdf(paste0(figDir, "aucs.pdf"), height=2, width=length(a)*2)
par(mfrow=c(1,length(a)), mar=c(3,3,1.5,.5), mgp=c(1.5,.5,0))
for (i in a)
{ roc(n[,i]>.25, loos[[i]], direction="<", plot=TRUE, mar=par("mar"), mgp=c(1.5,.5,0),
    print.auc=TRUE, print.auc.cex=1.3, smooth=FALSE, print.auc.y=.2, print.auc.x=.6)
  title(main=i);
}
dev.off();

pdf(paste0(figDir, "aucs2.pdf"), height=4, width=length(a)*2)
par(mfcol=c(2,length(a)), mar=c(3,3,1.5,.5), mgp=c(1.5,.5,0))
for (i in a)
{ roc(n[,i]>.25, loos[[i]], direction="<", plot=TRUE, mar=par("mar"), mgp=c(1.5,.5,0),
    print.auc=TRUE, print.auc.cex=1.3, smooth=FALSE, print.auc.y=.2, print.auc.x=.6)
  title(main=i);
  plotWcc(n[,i]*100, loos[[i]]*100, subSide=3, subLine=0, xlab=paste0("%", i, " annotated"),
    ylab=paste0("%", i, " predicted %"), pch='.', method='p')
}
dev.off();

f = function(i) factor(i, levels=c("IM", "BL", "M", "MSL", "LAR"));
pdf(paste0(figDir, "alluvialDec.2.pdf"));
par(mfrow=c(1,1), mar=c(3,3,0.5,.5), mgp=c(1.5,.5,0), xpd=NA)
#t = plyr::count(data.frame(Tum=f(tT), Tumor=f(cli$barTS.T), Full=f(cli$barPB), Stroma=f(cli$barTS.S),
#  Str=f(tS)));
t = plyr::count(data.frame(Tumor=f(cli$barT.an), PB=f(cli$barPB), Stroma=f(cli$barS.an),
  TIME=factor(cli$Immunophenotype.pathologist, levels=rev(c('ID', 'MR', 'SR', 'FI')))));
t = t[complete.cases(t),];
alluvial(t[,1:4], freq=t$freq, col=colBar[as.character(t$PB)], blockCol=c(colBar,  colTIME),
  draw_ticks=FALSE)
dev.off();

nfo = lapply(1:3, function(i) read.xlsx("~/Data/Spatial/TNBC/infos/List of signatures, genes of interest, cell type enrichment analysis computed in the study.xlsx", sheetIndex=i))
for (i in 1:3) { nfo[[i]][,2] = sub("^ *([^ ].+[^ ]) *$", "\\1", nfo[[i]][,2]); }
nfo[[1]][nfo[[1]][,2]=="VCpred_TN",2]="VCpred TN";
nfo[[1]][,2] = sub("^TNFA signaling.+", "TNFA signaling via NF-kB", nfo[[1]][,2]);
nfo[[2]][nfo[[2]][,2]=="mTOR",2]="MTOR";
for (i in 1:3) { rownames(nfo[[i]]) = nfo[[i]][,2]; }
rownames(nfo[[2]]) = rownames(geneList)[match(rownames(nfo[[2]]),
  geneList$List.of.interesting.genes)];
names(nfo) = c("Sigs", "Genes", "xCell");

f = function(x) factor(x, levels=names(colBar));
bl = paste(cli$barPB, cli$Immunophenotype.pathologist);
bl[!(bl %in% c("BL SR", "IM FI", "IM SR"))] = NA; bl[cli$barT.an!="BL"] = NA;
bl = factor(bl, levels=c("BL SR", "IM SR", "IM FI"))
me = paste(cli$barPB, cli$barS.an);
me[!(me %in% c("MSL MSL", "M MSL", "M M"))] = NA; me[cli$barT.an!="M"] = NA;
me2 = paste("Stroma", cli$barS.an); me2[cli$barPB!="M" | !(cli$barS.an %in% c("M", "MSL"))] = NA;
rownames(geneList) = geneList$symbol;
for (what in c("Sigs", "xCell", "Genes", "Annot", "Cli"))
{ for (version in c("Tumor", "Stroma", "BL", "M", "Mpb"))
  { v2 = c(Tumor2="Tumor", Tumor="Tumor", Stroma="Stroma", Stroma2="Stroma",
      BL="Stroma", M.Tumor="Tumor", M="Stroma", Mpb="Stroma")[version]
    cc = list(Tumor=list(Sigs=csAnT, SigsH=csAnT, xCell=xcAnT, Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol),]), Annot=cli$annotClean, Cli=ccn),
      Stroma=list(Sigs=csAnS, SigsH=csAnS, xCell=xcAnS, Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol),]), Annot=cli$annotClean, Cli=ccn))[[v2]][[what]]
    cmps = list(Tumor2=f(cli$barT.an), Tumor=f(cli$barT.an), Stroma=f(cli$barS.an), Stroma2=f(cli$barS.an), BL=bl,
      M=me, Mpb=me2)[[version]]
    if (version=="Tumor") { cmps[cmps=="IM"]=NA; }
    if (version=="Stroma") { cmps[cmps=="BL"]=NA; }
    if (what=="Genes" && v2=="Stroma")
    { cc = cc[,!(geneList[colnames(cc), "Role.of.genes"] %in% c("Methylation", "HRD", "Cycling", "MMR"))]
    }
    
    bubullePlot(cc, cmps, what=what,
      fileName=paste0(figDir, "TvsS.bar.dec.",version, ".pdf"))
  }
}

for (what in c("Sigs", "xCell", "Genes"))
{ for (version in c("Tumor", "Stroma", "BL", "M", "Mpb"))
  { v2 = c(Tumor2="Tumor", Tumor="Tumor", Stroma="Stroma", Stroma2="Stroma",
      BL="Stroma", M.Tumor="Tumor", M="Stroma", Mpb="Stroma")[version]
    cc = list(Tumor=list(Sigs=csAnT, SigsH=csAnT, xCell=xcAnT, Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol),])),
      Stroma=list(Sigs=csAnS, SigsH=csAnS, xCell=xcAnS, Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol),])))[[v2]][[what]]
    cc = cc[, grepl(c(Tumor="T", Stroma="S")[v2], nfo[[what]][colnames(cc), 3]),drop=FALSE]
    cmps = list(Tumor2=f(cli$barT.an), Tumor=f(cli$barT.an), Stroma=f(cli$barS.an), Stroma2=f(cli$barS.an), BL=bl,
      M=me, Mpb=me2)[[version]]
    if (version=="Tumor") { cmps[cmps=="IM"]=NA; }
    if (version=="Stroma") { cmps[cmps=="BL"]=NA; }
    #if (what=="Genes" && v2=="Stroma")
    #{ cc = cc[,!(geneList[colnames(cc), "Role.of.genes"] %in% c("Methylation", "HRD", "Cycling", "MMR"))]
    #}
    
    bubullePlot(cc, cmps, what=what,
      fileName=paste0(figDir, "TvsS.spec.",version, ".pdf"))
  }
}

# Tailored version...
aa = lapply(c("Tumor compartment", "Stroma compartment"), function(i)
{ a = read.xlsx("~/Data/Spatial/TNBC/infos/Molecular features selection.xlsx", i);
  a = a[!is.na(a$List_selected),]; a$List_selected = sub("^ *([^ ].+[^ ]) *$", "\\1", a$List_selected)
  a = a[order(factor(a$Main_classes, levels=c("Immune", "Stress response", "Stroma", "Pathway",
    "Metabolism", "Proliferation"))),]
  return(a);
}); names(aa) = c("T", "S");
col = c(colDisp, Proliferation=colDisp[["Oncogenic"]], Pathway=MC.colors[[1]]);

for (what in c("T", "S"))
{ a = aa[[what]];
  if (what=="T")
  { x = list(Sigs=csAnT, xCell=xcAnT, Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol), ]))
    cmps = f(cli$barT.an); cmps[cmps=="IM"]=NA; 
  }
  else
  { x = list(Sigs=csAnS, xCell=xcAnS, Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol), ]))
    cmps = f(cli$barS.an); cmps[cmps=="BL"]=NA;
  }
  x = lapply(names(x), function(i) x[[i]][,grepl(what, nfo[[i]][colnames(x[[i]]), 3]),drop=FALSE]) # Tumor/stroma features only
  colnames(x[[3]]) = geneList[colnames(x[[3]]), "List.of.interesting.genes"]
  cc = do.call(cbind, x);

  cmps=factor(cmps);

  p1 = calcP(cc, cmps); 

  ww = a$List_selected
  lbls=paste0("'", ww, "'");
  w = a$Feature_type=="Single genes"; lbls[w] = paste0('italic(', lbls[w], ')');
  lbls = parse(text=lbls);
  colLbls = c('grey45', 'black', 'darkblue'); names(colLbls) = c("Cell types", "Signatures", "Single genes");

  z = p.adjust(p1$p[ww,], method='fdr');
  cutOff = max(z[p1$fdr[ww,]<=.05])+1e-10;

  leftMar = 1+max(nchar(ww))/2; 
  sl = 0.06*max(nchar(colnames(p1$p))+8); sr = 0.06*(max(nchar(ww))+3)*.7; sb=0.06*leftMar*.7+2;
  cairo_pdf(paste0(figDir, c(T="tumor", S="stroma")[what], "Compartment.pdf"),
        height=(ncol(p1$p)*.45+sb)/1.5-.5, width=(2*(length(ww)+7)/7+sl+.5)/1.5)
  par(mai=c(.2, sl/2, sb/3, sr/2+.5), mgp=c(1.5,.5,0), cex=.6);
  dotPlot(cc[,ww,drop=FALSE], cmps, maxP=cutOff, oma=NULL, cex.pch=1, srt=30, horizontal=TRUE, axPos = 3, lbls=lbls,
    col.lbl=colLbls[a$Feature_type]);
  wi = which(!duplicated(a$Main_classes))[-1]
  rect(wi-.95, rep(par('usr')[4]-.2,length(wi)), wi-.05, rep(par('usr')[4]+.02, length(wi)), col='white', xpd=NA, border=NA);
  wi = c(1, wi, nrow(a)+1);
  mtext(txt<-a$Main_classes[wi[-length(wi)]+1], side=1, line=.3, at=(wi[-1]+wi[-length(wi)])/2-.5, cex=.7)
  for (i in 1:(length(wi)-1))
  { lines(wi[(0:1)+i]+c(-.3,-.7), c(.5,.5), xpd=NA, col=col[txt[i]], lwd=2)
  }
  if (what == "T")
  { legendDotPlot(length(ww)+4,par('usr')[4]+1.5, horizontal=FALSE, interline=1.6)
    #i = names(colLbls) %in% a$Feature_type;
    #legend(5,-1*strheight("M", cex=.6), legend=names(colLbls[i]), fill=colLbls[i], bty='n', ncol=3);
    legend(length(ww)+4,par('usr')[4]-2, legend=names(colLbls), fill=colLbls, bty='n', xjust=.5, border=NA);
  }
  dev.off();
}

# Survival
set.seed(123);
pdf(paste0(figDir, "survTSdec.2.pdf"), height=2.5, width=5);
par(mfrow=c(1,2), mar=c(3,3,1.5,.5), mgp=c(1.5,.5,0), cex=2/3)
plotSurv(cli$DRFS, bl, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="BL in tumor", legend.title="PB/TIME",
  legendPos='bottomleft', test="perm");
plotSurv(cli$DRFS, me, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="M in tumor", legend.title="PB/Stroma",
  legendPos='bottomleft', test="perm");
dev.off();

set.seed(123);
pdf(paste0(figDir, "survTSdec.pdf"), height=3, width=10);
par(mfrow=c(1,4), mar=c(3,3,1.5,.5), mgp=c(1.5,.5,0))
a = cli$barS.an; a[!(a %in% c("MSL", "M"))] = NA; a[cli$barT.an!="M"]=NA; a = factor(as.character(a));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="M in tumor", legend.title="Stroma",
  legendPos='bottomleft', test="perm", ppos=.1);
a = cli$barS.an; a[!(a %in% c("MSL", "M"))] = NA; a[cli$barPB!="M"]=NA;a = factor(as.character(a));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="M in PB", legend.title="Stroma",
  legendPos='bottomleft', test="perm", ppos=.1, nperm=1e5);
a = cli$barS.an; a[!(a %in% c("BL", "IM"))] = NA; a[cli$barT.an!="BL"]=NA;a = factor(as.character(a));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="BL in tumor", legend.title="Stroma",
  legendPos='bottomleft', test="perm", ppos=.1, nperm=1e5);
a = cli$barS.an; a[!(a %in% c("BL", "IM"))] = NA; a[cli$barPB!="BL"]=NA;a = factor(as.character(a));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="BL in PB", legend.title="Stroma",
  legendPos='bottomleft', test="perm", ppos=.1);
dev.off();

set.seed(123);
pdf(paste0(figDir, "survTSdec.final.1.pdf"), height=3, width=4);
par(mar=c(3,7,1.5,.5), mgp=c(1.5,.5,0))
a = cli$barS.an; a[!(a %in% c("MSL", "M"))] = NA; a[cli$barPB!="M"]=NA;
a = factor(as.character(a), levels=c("M", "MSL"), labels=c("'M-M'['STROMA']", "'M-MSL'['STROMA']"));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="M subtype in global PB",
 test="perm", ppos=.1, nperm=1e5, legendPos=NULL, parseCateg=TRUE, titleTable=NULL);
dev.off();

set.seed(123);
pdf(paste0(figDir, "survTSdec.final.2.pdf"), height=3, width=10);
par(mfrow=c(1,3), mar=c(3,9,1.5,.5), mgp=c(1.5,.5,0))
a = cli$barS.an; a[!(a %in% c("MSL", "M"))] = NA; a[cli$barT.an!="M"]=NA;
a = factor(as.character(a), levels=c("M", "MSL"), labels=c("paste(M['TUMOR'], '-M'['STROMA'])",
  "paste(M['TUMOR'], '-MSL'['STROMA'])"));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="M subtype in tumor PB",
  legendPos=NULL, test="perm", ppos=.1, parseCateg=TRUE);
a = cli$barS.an; a[!(a %in% c("BL", "IM"))] = NA; a[cli$barT.an!="BL"]=NA;
a = factor(as.character(a), levels=c("BL", "IM"), labels=c("paste(BL['TUMOR'], '-BL'['STROMA'])",
  "paste(BL['TUMOR'], '-IM'['STROMA'])"));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="BL subtype in tumor PB",
  legendPos=NULL, test="perm", ppos=.1, nperm=1e5, parseCateg=TRUE);
a = cli$barS.an; a[!(a %in% c("BL", "IM"))] = NA; a[cli$barPB!="BL"]=NA;
a = factor(as.character(a), levels=c("BL", "IM"), labels=c("'BL-BL'['STROMA']", "'BL-IM'['STROMA']"));
plotSurv(cli$DRFS, a, addTable=TRUE, xlab='Time (yr)', ylab="DRFS", main="BL subtype in global PB",
  legendPos=NULL, test="perm", ppos=.1, parseCateg=TRUE);
dev.off();

# Heatmap M / MSL
#####################
w = intersect(geneList$symbol, rownames(stromaAn)); x = t(stromaAn[w,]);
colnames(x) = geneList[w,"List.of.interesting.genes"];

x = list(csAnS[,intersect(names(sigInfo), colnames(csAnS))], xcAnS[,names(colXct)], x);
what = rep(c("Sig", "Xc", "Gene"), sapply(x, ncol))
x = do.call(cbind, x);
rownames(x) = rownames(cli);
p = moulinette(x, me2);
p2 = p.adjust(p, method='fdr');
y = x[, p2<.05]; wh = what[p2<.05]
sp = lapply(split(rownames(y), me2), function(i) as.character(sort(as.numeric(i))));
y = y[unlist(sp),];
rq = colQuantiles(y, probs=c(.05, .95));
y = (y-rep(rq[,1], each=nrow(y)))/rep(rowDiffs(rq), each=nrow(y))
y[y<0] = 0; y[y>1] = 1;

pdf(paste0(figDir, "TvsS.Mbp.hm.pdf"), height=5)
par(mar=c(3,10,1.5,2.5))
image(y[,ncol(y):1], col=colorRampPalette(c('blue', 'yellow'))(100), xaxt='n', yaxt='n');
wi = seq(0,1,length=ncol(y)); nd = which(!duplicated(wh))[-1];
mtext(rev(colnames(y)), side=2, at=wi, las=2, line=.5, font=1+2*(rev(wh)=="Gene")) 
mtext(rownames(y), side=1, at=seq(0,1,length=nrow(y)), cex=.5, line=-.3)
abline(h = wi[nd]-wi[2]/2, col='white', lwd=2) 
abline(v = (length(sp[[1]])-.5)/(nrow(y)-1), col='white', lwd=2)
text(c(length(sp[[1]])*.5, length(sp[[1]])+.5*length(sp[[2]])-1)/(nrow(y)-1),
  rep(par('usr')[4]+strheight("M"), 2), xpd=NA,
  parse(text=c("bold('M-M'['STROMA'])", "bold('M-MSL'['STROMA'])")), col=1:2)
sl = par('usr')[2]; em = strwidth("M");
for (i in 1:3)
{ lines(c(sl+em/2, sl+em/2), c(c(0, wi[nd])[i], c(wi[nd-1], 1)[i]), xpd=NA,
    col=c("#7031A0", "#3770C0", "#3CB050")[4-i], lwd=3)
}
text(rep(sl+em*1.5, 3),(c(0, wi[nd])+c(wi[nd-1], 1))/2, 
  c("Single genes", "Cell types", "Signatures"), font=2, xpd=NA, srt=270)
plotScale(colorRampPalette(c('blue', 'yellow'))(100), c("Low", "High"), c(0,1), .7,-strheight("M")*3,
  horizontal=TRUE,
  width=10, height=1)
dev.off();

################
## TLS
################

# TLS by gene
###############
pi = pinp;
ps = do.call(cbind, lapply(pi, function(i) i[,"p"]));
ms = do.call(cbind, lapply(pi, function(i) i[,"m"]));
rownames(ps) = rownames(ms) = rownames(xCan);
p = rowMaxs(ps); m = rowMins(ms); names(p) = names(m) = rownames(ps);

nfo = read.xlsx(paste0(dataDir, "misc/TLS signature.xlsx"), sheetName="TLS signature")
idG = nfo[,6]; names(idG) = nfo[,1]; idG=idG[1:30];
colG = palette.colors(8, "Tableau 10"); names(colG)=unique(idG); 

pdf(paste0(figDir, "tlsGene.pdf"), width=5, height=5)
g = names(idG); #g = names(ok)[ok]
m1 = rowMins(ms[,colnames(ms)!="Lymphocyte"])/log(2)
m2 = ms[,"Lymphocyte"]/log(2)
par(mar=c(3,3,3,3), mgp=c(1.5,.5,0));
z = smoothScatter(log2(1+2^m2), log2(1+2^m1), xlab="Fold-change TLS vs. lymphocyte compartment",
  ylab='Min fold-change TLS vs. another compartment',
  pch=NA, nrpoint=100, ret.selection = TRUE, xaxt='n', yaxt='n', transformation = function(x) {pmin(x,2) ^.25}) # 
a=c(0, .5, 1, 2, 4, 8, 16, 32); axis(1, at=log2(1+a), labels=a);
a=c(0, .5, 1, 2, 4, 8); axis(2, at=log2(1+a), labels=a);
abline(h=1, v=1, col='grey');

z =  union(z, which(names(p) %in% g));
co = rep('grey', length(z));
w = rownames(ms)[z] %in% g;
co[w] = colG[idG[rownames(ms)[z][w]]];

points(log2(1+2^m2[z]), log2(1+2^m1[z]), col=co)
#abline(h=-10, v=-10, col='red'); 
lbl=rownames(ms)[z]; lbl[co=="grey"]=NA;
#m1b = c(m1[z], rep(seq(from=-2, to=.9, len=10), 10), rep(seq(from=.5, to=3, len=10), 10))
#m2b = c(m2[z], rep(seq(from=-2, to=.9, len=10), each=10), rep(seq(from=-4, to=-1, len=10), each=10))
m1b = c(m1[z], rep(seq(from=-2, to=.9, len=10), 10), rep(seq(from=-3, to=1.5, len=10), 10))
m2b = c(m2[z], rep(seq(from=-2, to=.9, len=10), each=10), rep(seq(from=-2, to=5, len=10), each=10))
myAddTextLabels(log2(1+2^m2b), log2(1+2^m1b), c(lbl, rep(NA, 200)), cex.label=.8, col.label=c(co, rep(NA, 200)),
  heightPad=.3, widthPad=.1, figBorder=.03)

di = 1.6;
arrows(c(-.1, .1)+log2(2), par('usr')[4]+di*strheight("I"), par('usr')[1:2], par('usr')[4]+di*strheight("I"), xpd=NA, length=.17)
mtext(c("Lymphocytes", "TLS"), line=di/2+.7, at=(par('usr')[1:2]+1)/2, side=3)
arrows(a<-(par('usr')[2]+di*strwidth("N")), c(-.1, .1)+log2(2), a, par('usr')[3:4], xpd=NA, length=.17)
text(a+strwidth("N")*1.6, (par('usr')[3:4]+1)/2, xpd=NA, srt=270, c("Other non-lymphocytes","TLS")); 
legend('bottomright', legend=names(colG), fill=colG, bty='n', cex=.9);
dev.off();

# TLS by xCell
###############
ax = function(i)
{ a=axTicks(i); lbl = paste0("10^-", abs(a)); lbl[a==0] = 1;
  axis(i, at=a, labels=parse(text=lbl));
}
p = psGo$xCell; p=p[!grepl("Score", rownames(p)),];
p1 = -pmin(1, 2*rowMins(p[,c("Lymphocyte Higher", "Lymphocyte Lower")])) *
  sign(p[,"Lymphocyte Lower"]-p[,"Lymphocyte Higher"])
a = rowMins(p[,grepl("Lower", colnames(p)) & colnames(p) != "Lymphocyte Lower"])
p2a = -pmin(1, 2*pmin(a, p[,"Max other higher"])) * sign(a-p[,"Max other higher"])

pdf(paste0(figDir, "tlsXcl.pdf"));
par(mar=c(3,3,2.5,2.5), mgp=c(1.7,.5,0));
p2 = p2a;
plot(logAbs(p1), logAbs(p2), xlab="p-value TLS vs. lymphocyte compartment", col=colXct[names(p1)],
  ylab='p-value TLS vs. other non-lymphocytes', yaxt='n', xaxt='n')
ax(1); ax(2);
abline(h=0, v=0, col='grey')
di = 1.4;
arrows(c(-1, 1), par('usr')[4]+di*strheight("I"), par('usr')[1:2], par('usr')[4]+di*strheight("I"), xpd=NA)
mtext(c("Lymphocytes", "TLS"), line=di/2+.3, at=par('usr')[1:2]/2, side=3)
arrows(a<-(par('usr')[2]+di*strwidth("N")), c(-1, 1), a, par('usr')[3:4], xpd=NA)
text(rep(a+strwidth("N"), 2), par('usr')[3:4]/2, xpd=NA, srt=270,
  c("Other non-lymphocytes", "TLS"))

lbl = names(p1); lbl[abs(p1)>1e-4 & abs(p2)>1e-4] = NA;
myAddTextLabels(logAbs(p1), logAbs(p2), lbl, cex.label=.7, col.label=colXct[names(p1)],
  heightPad=0.25, widthPad=0.02)
legend('topleft', legend=names(colXc), fill=colXc, bty='n');
dev.off();

library(vioplot); library(colorspace);
co1 = colorspace::lighten('#405D92', .5); co2 = lighten('#E86E4D', .15)
p = psGo$xCell;
w = rownames(p)[p[,"Lymphocyte Higher"]<1e-4 | p[,"Lymphocyte Lower"]<1e-4];
w = w[!grepl("Score", w)];
w = w[order(-match(w, names(colXct)))]
pdf(paste0(figDir, "tlsViolin.pdf"), height=5);
par(mar=c(3,13,.5,.5), mgp=c(1.5,.5,0));
vioplot(xcAn[idAn[cli$hasTLS,"Lymphoid nodule"],w], horizontal=TRUE, side='right', las=1, col=co1,
  xlab="xCell enrichment score", plotCentre = "line", yaxt='n'); axis(1);
text(rep(par('usr')[1]-strwidth("n"), length(w)), seq_along(w), w, xpd=NA, adj=1, co=colXct[w])
vioplot(xcAn[idAn[,"Lymphocyte"],w], horizontal=TRUE, side='left', las=1, col=co2,
  add=TRUE, plotCentre = "line")
legend(.7, 19, legend=c("TLS", "Lymphocyte compartment"), fill=c(co1,co2), bty='n')
dev.off();

# tlsAllCli
###############
censor = function(s, t=10) { w = s[,1]>t; s[w,2]=0; s[w,1]=t; s; }
quartiles = function(x) { x = ceiling(4*rank(x)/length(x)); ordered(paste0("Q", x)); }
set.seed(123);
pdf(paste0(figDir, "tlsAllCli.pdf"), height=6, width=15)
par(mfcol=c(2,5), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), cex.main=1);
plotSurv(censor(cli$DRFS), quartiles(csPB[,'TLS ST']), addTable=TRUE, main="ST TNBC cohort", ylab="DRFS", xlab="Time (yr)",
  legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, args.legend=list(bty='n', inset=c(.02, 0)))
plotSurv(censor(ds$METABRIC$cli$DRFS), quartiles(ds$METABRIC$cs[,"TLS ST"]), addTable=TRUE, main="METABRIC - TNBC cohort", ylab="DRFS",
  xlab="Time (yr)", legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, args.legend=list(bty='n', inset=c(.02, 0)))
myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$wI & cliIspy$arm=="Pembro", ylab="TLS ST sig", subSide=3, subLine=-1.5,
  main="I-SPY2 - TNBC cohort - Pembrolizumab arm", addN=TRUE, alpha=.5)
myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$wI & cliIspy$arm=="Ctr", ylab="TLS ST sig", subSide=3, subLine=-1.5,
  main="I-SPY2 - TNBC cohort - control arm", addN=TRUE, alpha=.5)
for (i in c("pfs", "os"))
{ p = psSurv[[i]][[1]]; #p=p[rownames(p)!="TLS ST",]; # PFS, normal fig
  plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="HR", yaxt='n', xaxt='n',
    main=m<-paste("Other cancer types - immmunotherapy -", toupper(i)),
    col=1+(rownames(p)=="TLS ST"));
  abline(v=0, col='grey'); abline(h=-log10(0.05), col='grey');
  logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE, granularity=0)
  r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
  addAnnotArrows(side=1, xlab="HR", annotDir=c("Better", "Worse"), le=.05, cex=.6)
  nm = rownames(p); nm[p[,1]>.05]=NA; #nm[nm=="TLS_new"] = "TLS ST";
  myAddTextLabels(p[,2], -log10(p[,1]), labels=nm, col.label=1+grepl("TLS ST", nm), cex.label=.6, figBorder=.02);
  text(x=.0, y=par('usr')[4]-(2:4)*strheight("M")*1.5, pos=4, c(parse(text="bold('TLS ST')"),
    formatNiceP(p["TLS ST",1]), formatNiceP(exp(p["TLS ST",2]), pName="'HR'", italic=FALSE)))
  w = !is.na(suI[,i][,1]);
  plotSurv(suI[,i][w,], quartiles(dIm[w,"TLS ST"]), ylab=toupper(i), xlab="Time (yr)", marTable=1.5,
    addTable=TRUE, legend.title="TLS ST sig", leftTable=4, ppos=.35, main=m, args.legend=list(bty='n', inset=c(.02,0)))
}
p = pResp;
plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="OR", yaxt='n', xaxt='n',
  main="Other cancer types - immmunotherapy - RECIST",
  col=1+(rownames(pResp)=="TLS ST"));
nm = rownames(p); nm[p[,1]>.05]=NA; #nm[nm=="TLS ST"] = "TLS ST";
logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE)
r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
addAnnotArrows(side=1, xlab="OR", annotDir=c("Worse", "Better"), le=.05, cex=.6)
abline(v=0, col='grey'); abline(h=-log10(0.05), col='grey');
myAddTextLabels(p[,2], -log10(p[,1]), labels=nm, col.label=1+grepl("TLS ST", nm), cex.label=.6);
text(x=par('usr')[1]+strwidth("M")*3, y=par('usr')[4]-(2:4)*strheight("M")*1.5, pos=4, c(parse(text="bold('TLS ST')"),
    formatNiceP(p["TLS ST",1]), formatNiceP(exp(p["TLS ST",2]), pName="'OR'", italic=FALSE)))
myBoxplot(dIm[,"TLS ST"], c(R="Response", NR="No response")[as.character(resp)], ylab="TLS ST sig", subLine=-1,
  addN=TRUE, alpha=.5); title(main="Other cancer types - immmunotherapy - RECIST");
dev.off();

set.seed(123);
pdf(paste0(figDir, "tlsAllCliSmall.pdf"), height=6, width=9)
par(mfrow=c(2,3), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), cex.main=1);
plotSurv(censor(cli$DRFS), quartiles(csPB[,'TLS ST']), addTable=TRUE, main="ST TNBC cohort", ylab="DRFS", xlab="Time (yr)",
  legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, args.legend=list(bty='n', inset=c(.02, 0)))
plotSurv(censor(ds$METABRIC$cli$DRFS), quartiles(ds$METABRIC$cs[,"TLS ST"]), addTable=TRUE, main="METABRIC - TNBC cohort", ylab="DRFS",
  xlab="Time (yr)", legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, args.legend=list(bty='n', inset=c(.02, 0)))
plotSurv(censor(ds$`SCAN-B`$cli$DRFS), quartiles(ds$`SCAN-B`$cs[,"TLS ST"]), addTable=TRUE, main="SCAN-B - TNBC cohort", ylab="DRFS",
  xlab="Time (yr)", legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, , args.legend=list(bty='n', inset=c(.02, 0)))
myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$wI & cliIspy$arm=="Pembro", ylab="TLS ST sig", subSide=3, subLine=-1.5,
  main="I-SPY2 - TNBC cohort - Pembrolizumab arm", alpha=.9, addN=TRUE)
for (i in c("pfs"))#, "os"))
{ p = psSurv[[i]][[1]]; #p=p[rownames(p)!="TLS ST",]; # PFS, normal fig
  plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="HR", yaxt='n', xaxt='n',
    main=m<-paste("Other cancer types - immunotherapy -", toupper(i)),
    col=1+(rownames(p)=="TLS ST"));
  a = p.adjust(p[,1], 'fdr'); pl = max(p[a<.05,1])
  abline(v=0, col='grey'); abline(h=-log10(pl), col='grey');
  logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE, granularity=0)
  r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
  addAnnotArrows(side=1, xlab="HR", annotDir=c("Better", "Worse"), le=.05, cex=.6)
  nm = rownames(p); nm[p[,1]>pl]=NA; #nm[nm=="TLS_new"] = "TLS ST";
  o = order(p[,1]); adx = runif(50)*.1; ady = par('usr')[4]-runif(50)*strheight("M")*1.5*4
  myAddTextLabels(c(p[o,2],adx), c(-log10(p[o,1]),ady), labels=c(nm[o], rep(NA,50)),
    col.label=c(1+grepl("TLS ST", nm[o]), rep(NA, 50)), cex.label=.9,
    heightPad=0.3, widthPad=0.02);
  text(x=.0, y=par('usr')[4]-(1:3)*strheight("M")*1.5, pos=4, c(parse(text="bold('TLS ST')"),
    formatNiceP(p["TLS ST",1]), formatNiceP(exp(p["TLS ST",2]), pName="'HR'", italic=FALSE)))
}
p = pResp;
a = p.adjust(p[,1], 'fdr'); pl = max(p[a<.05,1])
plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="OR", yaxt='n', xaxt='n',
  main="Other cancer types - immunotherapy - RECIST",
  col=1+(rownames(pResp)=="TLS ST"));
nm = rownames(p); nm[p[,1]>pl]=NA; #nm[nm=="TLS ST"] = "TLS ST";
logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE)
r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
addAnnotArrows(side=1, xlab="OR", annotDir=c("Worse", "Better"), le=.05, cex=.6)
abline(v=0, col='grey'); abline(h=-log10(pl), col='grey');
o = order(p[,1]);
myAddTextLabels(p[o,2], -log10(p[o,1]), labels=nm[o], col.label=1+grepl("TLS ST", nm[o]), cex.label=.9,
    heightPad=0.3, widthPad=0.02);
text(x=par('usr')[1]+strwidth("M")*.5, y=par('usr')[4]-(1:3)*strheight("M")*1.5, pos=4, c(parse(text="bold('TLS ST')"),
    formatNiceP(p["TLS ST",1]), formatNiceP(exp(p["TLS ST",2]), pName="'OR'", italic=FALSE)))
dev.off();

set.seed(123);
pdf(paste0(figDir, "tlsAllCliRest.pdf"), height=3, width=15)
par(mfcol=c(1,5), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), cex.main=1);
p = pResp;
par(mar=c(2,3,2,.5))
myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$wI & cliIspy$arm=="Ctr", ylab="TLS ST sig", subSide=3, subLine=-1.5,
  main="I-SPY2 - TNBC cohort - control arm", addN=TRUE, alpha=.5)
myBoxplot(dIm[,"TLS ST"], c(R="Response", NR="No response")[as.character(resp)], ylab="TLS ST sig", subLine=-1,
  addN=TRUE, alpha=.5, main="Other cancer types - immmunotherapy - RECIST")
par(mar=c(3,3,2,.5))
for (i in c("pfs", "os"))
{ p = psSurv[[i]][[1]]; #p=p[rownames(p)!="TLS ST",]; # PFS, normal fig
  m = paste("Other cancer types - immmunotherapy -", toupper(i));
  w = !is.na(suI[,i][,1]);
  plotSurv(suI[,i][w,], quartiles(dIm[w,"TLS ST"]), ylab=toupper(i), xlab="Time (yr)",
    addTable=TRUE, legend.title="TLS ST sig", leftTable=4, ppos=.35, main=m, args.legend=list(bty='n', inset=c(.02,0)))
  if (i == "os")
  { plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="HR", yaxt='n', xaxt='n',
      main=m, col=1+(rownames(p)=="TLS ST"));
    abline(v=0, col='grey'); abline(h=-log10(0.05), col='grey');
    logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE, granularity=0)
    r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
    addAnnotArrows(side=1, xlab="HR", annotDir=c("Better", "Worse"), le=.05, cex=.6)
    nm = rownames(p); nm[p[,1]>.05]=NA; #nm[nm=="TLS_new"] = "TLS ST";
    myAddTextLabels(p[,2], -log10(p[,1]), labels=nm, col.label=1+grepl("TLS ST", nm), cex.label=.6, figBorder=.02);
    text(x=.0, y=par('usr')[4]-(2:4)*strheight("M")*1.5, pos=4, c(parse(text="bold('TLS ST')"),
      formatNiceP(p["TLS ST",1]), formatNiceP(exp(p["TLS ST",2]), pName="'HR'", italic=FALSE)))
  }
}
dev.off();

# Signatures from Charoentong on TLS
#####################################
w = idAn[,"Lymphoid nodule"]; w[!(rownames(cli) %in% idTLS[,"id"])]=NA;
z = data.frame(t(csX[,w]));
colnames(z) = sub("New.Immune_Charoentong_(.+)_CellRep.2017_PMID.28052254", "\\1", colnames(z));
colnames(z) = gsub("_", " ", colnames(z));

pdf(paste0(figDir, "sigCharoentongInTLS.FP.pdf"), width=7.5, height=6.5)
par(mar=c(.5,.5,.5,.5));
p = allForest(z, cli$DRFS, fdr=TRUE, annotDir=c("Better", "Worse"),sigOnFDR=TRUE)
dev.off();

pdf(paste0(figDir, "sigCharoentongInTLS.FPmv.pdf"), width=7.5, height=6.5)
par(mar=c(.5,.5,.5,.5));
p = allForest(data.frame(z, ctrl[rownames(cli),], check.names=FALSE), cli$DRFS, fdr=TRUE,
  annotDir=c("Better", "Worse"), control=ctrlForm,sigOnFDR=TRUE)
dev.off();

ww = colnames(z)[which(p$p<.05)];
f = function(x) { ordered(ceiling(rank(x, na.last="keep")/sum(!is.na(x))*3), labels=c("Low", "Medium", "High")) }
pdf(paste0(figDir, "sigCharoentongInTLS.pdf"), height=3.5, width=length(ww)*3)
par(mfrow=c(1,length(ww)), mar=c(3,4,.5,.5), mgp=c(1.5,.5,0))
for (i in ww)
{  plotSurv(cli[,"DRFS"], f(z[,i]), addTable=TRUE,
    xlab="Time (yr)", ylab="DRFS", legendPos='bottomleft', legend.title=i)
}
dev.off()

# AUCs of TLS sig
#####################
y = list(cli$annotClean[,"Lymphoid nodule"]>0, rownames(cli)%in%idTLS[,'id'])
pdf(paste0(figDir, "aucTLS.pdf"), width=4*2+.2, height=2*2+.4);
par(mfrow=c(2,4), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0), omi=c(0,.2,.4,0))
for (y2 in y)
{ y2 = factor(y2);
  for (i in 1:2)
  { x = list(csBulk, csPB)[[i]];
    x = x[,c("TLS ST", "TLS Lundeberg")]
    for (j in 1:2)
    { z = roc(response=y2, predictor=x[,j], direction="<", ci=TRUE);
      plot(z, mar=par('mar'), mgp=par('mgp'))
      mtext(paste("AUC:", round(z$auc, digits=2)), side=1, line=-3, adj=0, at=.6, cex=.8)
      mtext(paste0("(", round(z$ci[1], digits=2), "-",round(z$ci[3], digits=2), ")"),
        side=1, line=-1.4, adj=0, at=.6, cex=.8)
    }
  }
}
mtext(c("Annotation", "Regression-corrected"), side=2, at=c(.75, .25), outer=TRUE)
mtext(c("Bulk", "PB"), side=3, at=c(.25, .75), outer=TRUE, line=1.5)
mtext(rep(c("ST sig", "Lundeberg sig"), 2), side=3, at=(.5:3.5)/4, outer=TRUE, line=0.3, cex=.8)
dev.off();

# TLS in iSPY
##############
x = makeAllSubgroups(csI[,"TLS ST"], cliIspy[,c("Receptor.Subtype", "Arm..short.name.")]);
tr = sub(".+ ", "", colnames(x));
w = !grepl(" ", colnames(x)); tr[w] = colnames(x)[w]; tr[1:5]='';
li = c(1, 1+seq_along(tr)+unclass(factor(tr))*.5)*1.5; 
#m = tapply(li[-1], tr, mean);


cairo_pdf(paste0(figDir, "tlsIspy.pdf"), height=10, width=7.5)
par(mar=c(.7,.7,.7,.7));
allForest(x, cliIspy$pcr, useWilcox=TRUE,
  clip=c(.1, 10), boxsize=1, lineHeight=li, annotDir=c("Worse", "Better"), sigOnFDR=TRUE)
#ch = abs(strheight("!")); 
#mtext(names(m), side=2, at=m*ch)
dev.off();

# TLS vs. GO
#############
# Note: based on GOSS here...
f = function(i)
{ i=sub("GO_", "", i);
  i=tolower(gsub("_", " ", i));
  i=gsub("(^| )([bt]) ", "\\1\\U\\2\\E-", i, perl=TRUE);
  i=gsub("(^| )i ", "\\1I ", i);
  i=gsub("(^| )(igg|cd[0-9]+|adp|mhc|nadh|atp) ", "\\1\\U\\2\\E ", i, perl=TRUE)
  i=gsub("v d j","V/D/J", i);
  i=gsub("g2 m","G2/M", i);
  return(i);
}
cutText = function(x, cex=1)
{ s = strwidth(x, cex=cex)
  w = which(s>par('usr')[2])
  for (i in w)
  { z = x[i];
    while(strwidth(z, cex=cex) > par('usr')[2]-strwidth("\U2026", cex=cex))
    { z = sub(" [^ ]+$", "", z); } 
    x[i] = paste0(z, "\U2026");
  }
  return(x);
}
cols=c('lightblue', 'yellow')
p = cbind(p.adjust(psGo$GO[,"Lymphocyte Higher"], 'fdr'), p.adjust(psGo$GO[,"Lymphocyte Lower"], 'fdr'),
  p.adjust(psGo$GO[,"Max other higher"], 'fdr'), p.adjust(psGo$GO[,"Max other lower"], 'fdr'));
for (j in 1:2)
{ cairo_pdf(paste0(figDir, "tlsGO", c(".VSlympho", ".VSrest")[j], ".pdf"),
    height=12, width=5)
  par(mar=c(3,.5,3,1), mgp=c(1.5,.5,0))
  i = goss[[j]];
  at = barplot(-log10(abs(i[,2])), horiz=TRUE, yaxt='n', xlab=expression(paste(italic('P'), '-value')),
    col=cols[2-i[,3]],
    yaxs='i', ylim=c(0, 1.2*nrow(i)), xaxt='n');
  a = axTicks(1); lbl = paste0("10^-", a); lbl[a==0]=1; axis(1, at=a, label=parse(text=lbl));
  text(0, at, cutText(paste0(" ",f(i[,1])), cex=.9), adj=0, cex=.9)

  title(c("TLS vs. lymphocyte compartment", "TLS vs. other non-lymphocytes")[j], line=2);
  mtext(list(c("TLS", "Lymphocyte compartment"), c("TLS", "Other non-lymphocytes"))[[j]], side=3,
      col=colorspace::darken(cols,.4), line=.5, adj=c(.25,.75), cex=par('cex')) 
  dev.off();
}

## TLS immuno
##############
s = resp;
w = !is.na(s) & ty!="Liver";
st2=st[w]; ty2=ty[w]; x = dIm[w,"TLS ST"]; s=s[w];
N = table(ty2);
pR = do.call(rbind, tapply(seq_along(ty2), ty2, function(w)
{ if (length(unique(st2[w]))>1) { coef(summary(glmer(s ~ x+(1|st2), family=binomial, subset=w))) [2,] }
  else { coef(summary(glm(s ~ x, family=binomial, subset=w))) [2,] }
}))
colnames(pR) = c("coef", "se(coef)", "z", "Pr(>|z|)")
pR = list(pR, N);

pS = lapply(colnames(suI), function(i)
{ s = suI[,i];
  w = !is.na(s[,1]) & ty!="Liver";
  st2=st[w]; ty2=ty[w]; x = dIm[w,"TLS ST"]; s=s[w,];
  N = table(ty2);
  p = do.call(rbind, tapply(seq_along(ty2), ty2, function(w) coef(summary(coxph(s~x+strata(st2), subset=w)))[1,]))
  list(p, N);
}); names(pS) = colnames(suI);

pp = c(pS,RECIST=list(pR));

o = setdiff(rownames(pp[[3]][[2]])[order(-pp[[3]][[2]])], 'Unknown');
o = sub("_", " ", o);
td = sub("_", " ", setdiff(names(pp[[3]][[2]]), "Unknown"));
li = 1:(length(td)+1); names(li)=c("", o); li=li*1.8;

pdf(paste0(figDir, "fpTLSimmuno.pdf"), width=4.5, height=10)
par(mar=c(.5,1,2.5,1));
layout(matrix(1:3, ncol=1), height=c(.8,1,.9))
par(cex=1);
for (i in c(1, 3, 2))
{ p = pp[[i]][[1]]; N = pp[[i]][[2]];
  w = rownames(p) != "Unknown"; p=p[w,]; N=N[w];
  rownames(p) = names(N) = sub("_", " ", rownames(p));
  m = intersect(o, names(N)); N=N[m]; p=p[m,];
  m =cbind(lower=p[,"coef"]-1.96*p[,"se(coef)"], mean=p[,"coef"],
      upper=p[,"coef"]+1.96*p[,"se(coef)"])
  m = exp(m);
  col = sapply(p[,"Pr(>|z|)"], function(z) if (!is.na(z) && z < .05) { "#00AFBB" } else { "slategray4" })
  d = list(rownames(p), N=as.vector(N), P=parse(text=sapply(p[,"Pr(>|z|)"], formatNice, parse=FALSE)), #p[,"Pr(>|z|)"],
      HR=sapply(exp(p[,"coef"]), signif, 2))
  names(d)[1]="";
  ad = c("Better", "Worse"); if (i==3) { ad=rev(ad); }
  basicForest(d, m, xlim=c(.2, 5), xlog=TRUE, xlab=c("HR","HR", "OR")[i], col=col, li=2,
    annotDir=ad, cexAnnot=.75, titles=c("", "N", expression(bolditalic("p")), "HR"))
  title(main=toupper(names(pp)[i]))
}
dev.off();

#######################
## MC & ET
#######################

## Bareche by cluster
######################
t = table(idC$bar, idC[,1]); t = t[,rownames(cli)];
o = order(cli$barPB, -t[cbind(as.character(cli$barPB), colnames(t))]/colSums(t))
#colTime = colorRampPalette(c("white", "red"), space="Lab")(4)

pdf(paste0(figDir, "barByClust.pdf"), height=4, width=15)
par(mar=c(5,9,.5,.5), mgp=c(1.5,.5,0));
at = barplot(t[,o], col=colBar, xaxt="n", ylab="N clusters", xaxs='i', xlim=c(-1, ncol(t)*1.2+1))
mtext(rownames(cli)[o], line=0, cex=.6, side=1, at=at)
addAnnot("Subtype PB", colBar[cli$barPB[o]], line=1, side=1, at=at)

b = cli$barPBco[o,];
b = b/rowMaxs(b); b2 = apply(b,1, function(i) { o = order(i); c(o[4], i[o[4]]); })
b2[2,b2[2,]<0] = 0;
addAnnot("Subtype PB 2nd class", colBar[colnames(b)[b2[1,]]], line=2, side=1, at=at, heights=b2[2,])

f = cli$Immunophenotype.pathologist; f[f=="nd"]=NA; f=factor(f, levels=c("ID", "MR", "SR", "FI"));
addAnnot("TIME", colTIME[f[o]], line=3, side=1, at=at)
legend(-15, 7.5, names(colBar), fill=colBar, xpd=NA, bty='n', xjust=0, yjust=1)
text(-14, 7.7, 'Subtype', xpd=NA, adj=0);
legend(-15, 3, levels(f), fill=colTIME, xpd=NA, bty='n', xjust=0, yjust=1)
text(-14, 3.2, 'TIME', xpd=NA, adj=0);
dev.off();

# MC descriptions
###################
anByC = do.call(rbind, tapply(seq_along(K), K, function(i) colMeans(annot2[i,], na.rm=TRUE)));
anByC = anByC[14:1,c("Tumor", "Stroma", "Lymphocyte", 
  "Fat tissue", "Necrosis", "in situ", "Lymphoid nodule", "Vessels", "Lactiferous duct")]
co = colAnn2[colnames(anByC)]; co[2] = colAnn2["Stroma cell"]; names(co)[2]="Stroma";
descrs=c("High proliferation, low immune, Notch, Wnt, AR",
  "High proliferation, low immune, DNA repair, Notch",
  "High proliferation, low immune, DNA repair, Notch,
  high glycolysis, hypoxia",
  "Moderate proliferation, low immune, low PI3K, low KRAS",
  "High proliferation, moderate immune, high oxidative
   phosphorylation, low angiogenesis, low EMT",
  "p53, PI3K/AKT/mTOR",
  "Moderate/high immune, DNA repair",
  "AR, ER signaling, fatty metabolism, high glycolysis, 
  oxidative phosphorylation, immune low",
  "High immune +++, high apoptosis",
  "High immune, low proliferation",
  "High EMT, high angiogenesis, low proliferation",
  "Low immune, low Notch, moderate proliferation",
  "Low DNA repair, high apoptosis, high TGFb,
    high Hedgehog",
  "High apoptosis, high angiogenesis, high EMT,
  high JAK/STAT, low DNA repair, low proliferation");
descrs = read.xlsx("~/Data/Spatial/TNBC/Divers/megaclusters_characterization.xlsx", 1)$Summary;
descrs = descrs[!is.na(descrs)];
w = gregexpr(", ", descrs)
for (i in seq_along(w))
{ ww = c(w[[i]], nchar(descrs[i]));
  if (all(ww<50)) { next; }
  cu = ww[which(ww>=50)-1];
  descrs[i] = paste0(substr(descrs[i], 1, cu), "\n", substr(descrs[i], cu+2, nchar(descrs[i]))) 
}

sig2keep = unlist(read.xlsx("~/Data/Spatial/TNBC/Divers/megaclusters_characterization.xlsx", 2))
sig2keep = sub(" $", "", sig2keep[!is.na(sig2keep)]);
sig2keep[sig2keep=="VCpred_TN"] = "VCpred TN";
sig2keep[sig2keep=="TNFA signaling via NF−kB"] = "TNFA signaling via NF-kB"

x = xCn.cs;
colnames(x)[colnames(x)=="PI3K AKT mTOR signaling"] = "PI3K/AKT/mTOR signaling"
x = x[,sig2keep];
x = (x-rep(colMeans(x), each=nrow(x)))/rep(colSds(x), each=nrow(x))
x = x[order(K),]; Ko=sort(K);
for (i in 1:14)
{ w = Ko==i;
  h = hclust(as.dist(1-cor(t(x[w,]))), 'average')
  x[w,] = x[w,][h$order,];
}
x2 = x; x2[x2< -2] = -2; x2[x2>2] = 2;
col = colorRampPalette(c('blue', 'yellow'))(100);

pres = calcPres(K, idC);
addT = 20; til = apply(pres,2,function(i) cli$sTILs.pathologist.percentage[i]+addT)
#idC$bar = factor(TNBCclassif(y, version='bareche', shortName=TRUE), levels=names(colBar));

pdf(paste0(figDir, "tblsMC2.pdf"), width=14, height=4.5)
par(mar=c(7,3,.5,.5), mgp=c(1.5,.5,0));
layout(mat=matrix(1:6, nrow=1), width=c(.5,1,1,3,1.6,.5))
a = rev(colSums(pres)); names(a) = paste0("MC", names(a)); names(a) = sub("k", "", names(a));
at = barplot(a, horiz=TRUE, xlab="N patients", xaxt='n', las=2, col=rev(MC.colors))
axis(1)
par(mar=c(7,.5,.5,.5))
barplot(t(anByC)*100, xlab="% Annotation", horiz=TRUE, yaxt='n', col=co)
usr = par('usr');
t = table(idC$bar, K); t = t/rep(colSums(t), each=5); t=t[names(colBar),];
barplot(t[,14:1]*100, col=colBar, xlab="% TNBC subtype", horiz=TRUE, yaxt='n') 

plot.new(); plot.window(xlim=c(0,1), ylim=usr[3:4], yaxs='i', mar=c(3,0,.5,.5), xaxs='i');
for (i in 1:14)
{ subplot(image(t(x2[Ko==i,]), axes=FALSE, zlim=c(-2, 2), col=col), x=c(0,1), rev(at)[i]+c(-.5,+.5)); 
}
text(seq(0,1,len=ncol(x2)+1)[-1]-.5/ncol(x2), 0, colnames(x2), xpd=NA, srt=-40, adj=0, cex=.7,
  col=colSig[sigInfo[colnames(x2)]]) 

plot.new(); plot.window(xlim=c(0,1), ylim=usr[3:4], yaxs='i', mar=c(3,0,.5,.5));
par(lheight=.7)
text(rep(-.05, 14), at, rev(descrs), adj=0, xpd=NA)

subplot(image(matrix(1:25, ncol=1), axes=FALSE, col=col, xpd=NA, useRaster=TRUE), x=c(.2,.8), y=usr[3]-3+c(1,2))
text(x=c(.2, .8), y=rep(usr[3]-.7, 2),c("Low", "High"), xpd=NA, pos=c(4,2));

plot.new(); par(mar=c(0,0,.5,0)); plot.window(xlim=c(0,1), ylim=c(0,1));
legend('topleft', legend=colnames(anByC), fill=co, title="Annotations", xpd=NA, inset=c(-.4, .1))

dev.off();

# MC characterisations (signatures etc)
##########################################
K2 = factor(K, levels=1:14, labels=paste0("MC", 1:14));
for (what in c("Sigs", "xCell", "Genes"))
{ cc = list(Sigs=xCn.cs, Genes=t(xCn), xCell=xCn.xc)[[what]];
  co = NULL;
  if (what=="Genes")
  { o = intersect(rownames(geneList), colnames(cc)); cc = cc[,o];
    co = colGenes[o];
    names(co) = colnames(cc) = geneList[colnames(cc), "List.of.interesting.genes"]
  }
  if (what=="xCell")
  { o = order(-match(colnames(cc), names(colXct))); cc = cc[,o];
    co = colXct[colnames(cc)]; names(co) = colnames(cc);
  }
  if (what=="Sigs")
  { cc = cc[,intersect(names(sigInfo), colnames(cc))];
    co = colSig[sigInfo[colnames(cc)]]; names(co) = colnames(cc);
  }

  p1 = calcP(cc, K); 
  p = rowMins(p1$fdr);  
  
  w = colnames(cc)[which(p<=.05)];
  if (length(w)==0) { next; }
  
  z = p.adjust(p1$p[w,], method='fdr');
  cutOff = max(z[p1$fdr[w,]<=.05])+1e-10;
  
  cairo_pdf(paste0(figDir, "MCcomp.", what, ".pdf"),
      height=(length(w)+7)/7+1, width=6+10+max(nchar(w))/20)
  par(mar=c(8,max(nchar(w)+4)/2.5,.5,.5), mgp=c(1.5,.5,0), cex=1.1);
  layout(matrix(1:2, nrow=1), width=c(.4,.6));
  
  dotPlot(cc[,w,drop=FALSE], K2, col.lbl=co[w], maxP=cutOff, oma=NULL, cex.pch=.45, srt=45);
  legendDotPlot(.5,-6*strheight("M"), horizontal=TRUE, cex.pch=.45)

  
  a = cc[,w]; o = order(K2); ko = K2[o]
  a = t(a[o,]);
  annots = co[w];
  a = a-rowQuantiles(a, probs=.05); a=a/rowQuantiles(a, probs=.99); #rowMaxs(a);
  a[a<0] = 0; a[a>1]=1;
  a = a[names(annots),];
  b=t(a); 
   
  #par(mar=c(3,max(nchar(rownames(a))+4)/2.5,.5,.5));#a = a[nrow(a):1,]; 
  o = c(3, 4, 1, 2);
  image(b, xaxt='n', yaxt='n', col=colorRampPalette(c('blue', 'yellow'))(100), xaxs='i', yaxs='i')
  pU = par('usr')[o];
  mtext(rownames(a), at=midPoint(pU[1:2], nrow(a)), adj=1, side=2, line=.5,
      col=annots[rownames(a)], las=2, cex=par('cex'))
  z = rle(as.integer(ko))
  cu = cumsum(z$lengths); cu=cu/cu[length(cu)]
  abline(v=cu, col='white', lwd=2)

  text(x = (c(0, cu[-length(cu)])+cu)/2, y = rep(par('usr')[1], length(z$values))-strheight("M"),
    paste0("MC", z$values), srt=45, adj=1, xpd=NA)
  
  plotScale(colorRampPalette(c('blue', 'yellow'))(100), c("Lower", "Higher"), 0:1,
    .5-strwidth("M")*10, -strheight("M")*10, horizontal=TRUE, width=20, height=2)
  
  dev.off();
}

# Number of MC
################
ploo = lapply(d<-dir(paste0(dataDir, "clustering/looMC"), pattern="RData"), function(i)
{ load(paste0(dataDir, "clustering/looMC/", i)); return(pLoo[,2]);
}); names(ploo) = sub("k([0-9]+).RData", "\\1", d);
ploo = ploo[order(as.integer(names(ploo)))];

pdf(paste0(figDir, "meanAUC.pdf"), height=3.5, width=3.5)
par(mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
plot(sapply(ploo, length), sapply(ploo, mean), ylab="Mean AUC", xlab="N megaclusters", col=(names(ploo)=="15")+1)
abline(v=14, col='grey')
dev.off();


## Number of ET
###############
pdf(paste0(figDir, "Ncluster.pdf"), height=3.5, width=4);
par(mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
nb = NbClust2(log10(100*(faN)), distance = "euclidean", min.nc = 5, max.nc = 15,
    method="ward.D2", index ="allNonGraph"); 
barplot(table(nb$Best.nc[1,]), ylab="N indices", xlab='N ecotypes');

dev.off();


## Heatmap ET
################
hc = hclust(dist(log10(100*(faN))), method='ward.D2');
cu = which(diff(cutree(hc, max(ecot))[hc$order])!=0);

ff = faN; colnames(ff) = paste0("MC", colnames(ff));
pdf(paste0(figDir, "heatmapMC.pdf"), width=15, height=10)
h = heatmap.3(log10(100*(t(ff))), scale='none', col=colorRampPalette(c('blue', 'yellow'), space='Lab')(100),
  Rowv=FALSE, 
  Colv=as.dendrogram(hc), colsep=cu,  dendrogram='column',
  zlim=c(0,2), ColSideColors=cbind(Subtype=colBar[cli$barPB],
    #`TIME (expression)`=colTIME[c(`Margin restricted`='ID', `Fully Inflamed`='FI', `Stroma Restricted`='SR')[cli$TIMEpb]],
    `TIME`=colTIME[cli$Immunophenotype.pathologist],
    ET.colors[ecot], ET.colors[ecot]),
    lhei=c(2,10), lwid=c(2,10), cexRow=1.8,
  margins=c(2,5), trace='none', sepcolor="white", side.height.fraction=1)#,
par(new=TRUE); plot.new();
layout(matrix(1:4, ncol=2), height=c(2,10), width=c(2,10))
par(cex=0.2 + 1/log10(94)); par(mar=c(0,0,0,5))#, cex=10);
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1), xaxs='i');
w = c(0, cumsum(tabulate(ecot))); w = (w[-1]+w[-length(w)])/(2*w[length(w)])
w = w*94.5/94;
mtext(paste0("ET", 1:max(ecot)), side=3, cex=1, col='black', at=w, font=2, line=-.8)
mtext("Ecotype   ", side=3, cex=1, col='red', at=0, font=2, line=-.5, adj=1);
par(new=TRUE, mar=c(0,0,0,0), xpd=NA, cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1), xaxs='i');
#legend(-.15, .9, legend=c("Relapse", "No relapse"), fill=c('red', 'black'), xjust=0, title="    Relapse", bty='n',
#  title.adj=0)
legend(-.15, .9, legend=names(colTIME), fill=colTIME, xjust=0, title="    TIME", bty='n', title.adj=0)
legend(-.15, .6, legend=names(colBar), fill=colBar, xjust=0, title='    Subtype', bty='n', title.adj=0)
dev.off();       

## Number of ET
################
nb = NbClust2(log10(100*(faN)), distance = "euclidean", min.nc = 5, max.nc = 15,
    method="ward.D2", index ="allNonGraph"); 
pdf(paste0(figDir, "N_ET.pdf"), height=3, width=4);
par(mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
barplot(table(nb$Best.nc[1,]), ylab="N indices", xlab="N Ecotypes");
dev.off();

a = split(colnames(nb$Best.nc), nb$Best.nc[1,])
a["0"]=NULL
invisible(sapply(names(a), function(i) cat(i, ":", a[[i]], "\n")))

# Content ET
#############
n = do.call(cbind, tapply(seq_along(ecot), ecot, function(i) colMeans(fa[i,])))
n = 100*n/rep(colSums(n), each=nrow(n));
colnames(n) = paste0("ET", colnames(n));

pdf(paste0(figDir, "barplotMCvsET.pdf"), height=5, width=6)
par(mar=c(2,3,1,5), mgp=c(1.5,.5,0));
at = barplot(n, col=MC.colors, ylab="% MC")
legend(x=par('usr')[2], y=par('usr')[4], paste0("MC", 1:nrow(n)), fill=MC.colors, xpd=NA)
dev.off();

pdf(paste0(figDir, "barEcot.pdf"), width=12, height=4);
par(mfcol=c(2,4), mar=c(2,3,1.5,.5), mgp=c(1.5,.5,0))
for (i in 1:4)
{ e = paste0("ET", list(ecot, ecotRec[stud=="S"], ecotRec[stud=="M"], ecotRec[stud=="SC"])[[i]]);
  b = list(cli$barPB, ds$`Spatial TNBC`$cli$bar, ds$METABRIC$cli$bar, ds$`SCAN-B`$cli$bar)[[i]]
  nm = c("ST global pseudobulk", "ST bulk", "METABRIC", "SCAN-B")[i]
  b = factor(b, levels=names(colBar));
  t = table(b,e);
  h = barplot(t, col=colBar, main=nm, ylab="N samples", xaxt='n')
  mtext(colnames(t), side=1, at=h, line=par("mgp")[2], cex=par('cex'))#, col=ET.colors);
  t=t/rep(colSums(t), each=nrow(t));
  barplot(100*t, col=colBar, main=nm, ylab="% samples", xaxt='n')
  mtext(colnames(t), side=1, at=h, line=par("mgp")[2], cex=par('cex'))#, col=ET.colors);
}
dev.off();


## Forest plot survival MC and ET
#####################################
ctr = list(METABRIC=ctrlMB, `Spatial TNBC`=ctrl[rownames(ds[["Spatial TNBC"]]$cli),], `SCAN-B`=ctrlSB)
ctrTot = do.call(rbind, ctr[names(ds)])
ct = update(ctrlForm, ~.+strata(stud))

colWidth = c(1.5,1.5,4,2.5,2,5);
for (i in c("PB", "All3", names(ds)))
{ if (i=="All3")
  { sus = susDs;
    x = do.call(rbind, lapply(ds, function(i) i[["15"]]$H));
    e = ecotRec;
  }
  if (i=="PB")
  { sus = cli[,sapply(cli, is, "Surv")]; sus$iDFS=NULL;
    x = fa;
    e = ecot;
  }
  if (i %in% names(ds))
  { x = ds[[i]]$cli;
    sus = x[,sapply(x, is, "Surv")];
    if (ncol(sus)>4) { sus = sus[,2:5]; }  # Remove RFS for Bulk
    x = ds[[i]][["15"]]$H;
    e = ds[[i]]$ecot[,1];
  }
    
  x[x<.01] = .01; 
  colnames(x) = paste0("MC", 1:ncol(x));
  x = scale(log10(x)); 
  
  x2 = matrix(FALSE, nrow=nrow(x), ncol=max(e));
  x2[cbind(1:nrow(x2), e)]=TRUE;
  colnames(x2) = paste0("ET", 1:ncol(x2));
  
  ctr = list(All3=ctrTot, PB=ctrl[rownames(cli),], METABRIC=ctrlMB, `Spatial TNBC`=ctrl[rownames(ds[[i]]$cli),], `SCAN-B`=ctrlSB)[[i]];

  nm = c(All3="All3", PB="PB", METABRIC="MB", `Spatial TNBC`="Bulk", `SCAN-B`="SCANB")[i]
  nm2 = c(All3="All 3 studies", PB="Pseudo-bulk", METABRIC="METABRIC", `Spatial TNBC`="ST cohort (bulk RNA-seq)", `SCAN-B`="SCAN-B")[i]
  
  if (i=="All3") { x = data.frame(x, stud=stud); x2 = data.frame(x2, stud=stud); ct1 = ~strata(stud); ct2=ct; }
  else { ct1 = ~1; ct2 = ctrlForm; }
  if (nm=="Bulk") { clip = c(.2, 5) } else { clip = c(.5, 2); }
  
  cairo_pdf(paste0(figDir, "forestMC.", nm, ".pdf"), width=8.5, height=13.5)
  par(mar=c(.1,1,1,1));
  layout(mat=matrix(1:8, ncol=2, byrow=TRUE), heights=rep(c(1.05, 1), 4))
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(x, sus[[j]], control=ct1, clip=clip, annotDir=c("Better", "Worse"),
      lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Univariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  if (length(sus)==3) { plot.new(); }
  
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(data.frame(x, ctr), sus[[j]], control=ct2,
      clip=clip, annotDir=c("Better", "Worse"), lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Multivariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  dev.off();
  
  cairo_pdf(paste0(figDir, "forestET.", nm, ".pdf"), width=8.5, height=11.5)
  par(mar=c(.1,1,1,1));
  layout(mat=matrix(1:8, ncol=2, byrow=TRUE), heights=rep(c(1.05, 1), 4))
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(x2, sus[[j]], control=ct1, clip=clip, annotDir=c("Better", "Worse"),
      lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Univariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  if (length(sus)==3) { plot.new(); }
  
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(data.frame(x2, ctr), sus[[j]], control=ct2,
      clip=clip, annotDir=c("Better", "Worse"), lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Multivariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  dev.off();
}

## KM Survival Ecotypes (OK)
###############################
bar = unlist(lapply(ds, function(i) as.character(i$cli$bar)));
a = ifelse(ecotRec!=4 & bar=="IM", "IM not ET4", "Others") 
et4im = ifelse(ecotRec==4 & bar=="IM", "IM ET4", "Others")
et4 = ifelse(ecotRec==4, "ET4", "Others")
b = a; b[et4im!="Others"]=et4im[et4im!="Others"];
pdf(paste0(figDir, "survEcotypeNew2.pdf"), height=7.5, width=7.5);
par(mfrow=c(3,3), mar=c(2.5,4.5,1.5,.5), mgp=c(1.5,.5,0), oma=c(0,1,0,0))
for (i in names(susDs))
{ su = susDs[[i]];
  plotSurv(su, et4, ylab=toupper(i), xlab="Time (yrs)",
    legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET4"], "black"),
    pval = getCoxP(coxph(su ~ et4+strata(stud))))
  #plotSurv(su, ifelse(ecotRec==3, "ET3", "Others"), ylab=toupper(i), xlab="Time (yrs)",
  #  legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET3"], "black"),
  #  pval = getCoxP(coxph(su ~ (ecotRec==3)+strata(stud))))  
  #plotSurv(su, a, ylab=toupper(i), xlab="Time (yrs)",
  #  legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=2:1,
  #  pval = getCoxP(coxph(su ~ a+strata(stud))))
  #plotSurv(su, ifelse(ecotRec==3, "ET3", "ET4"), subset=ecotRec%in%c(3,4), ylab=toupper(i), xlab="Time (yrs)",
  #  legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=ET.colors[c("ET3", "ET4")],
  #  pval = getCoxP(coxph(su ~ ecotRec+strata(stud), subset=ecotRec%in%3:4)))
  plotSurv(su, b, subset=b!="Others", ylab=toupper(i), xlab="Time (yrs)",
    legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET4"], palette()[2]),
    pval = getCoxP(coxph(su ~ b+strata(stud), subset=b!="Others")))
  plotSurv(su, ifelse(ecotRec==8, "ET8", "Others"), ylab=toupper(i), xlab="Time (yrs)",
    legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET8"], "black"),
    pval = getCoxP(coxph(su ~ (ecotRec==8)+strata(stud))))
}
dev.off();

for (i in c("ST", names(ds)))
{ if (i!="ST")
  { z = ds[[i]]; nm = c(METABRIC="MB", `Spatial TNBC`="Bulk", `SCAN-B`="SCANB")[i];
    e = z$ecot[,1]; cl = z$cli;
  } else { e=ecot; cl=cli; nm="ST"; }
  a = ifelse(e!=4 & cl$bar=="IM", "IM not ET4", "Others") 
  et4 = ifelse(e==4, "ET4", "Others")
  et4im = ifelse(e==4 & cl$bar=="IM", "IM ET4", "Others")
  b = a; b[et4im!="Others"]=et4im[et4im!="Others"];
  w = names(cl)[sapply(cl, is, "Surv")];
  pdf(paste0(figDir, "survEcotype.", nm, ".pdf"), height=2.5*length(w), width=7.5);
  par(mfrow=c(length(w),3), mar=c(2.5,4.5,1.5,.5), mgp=c(1.5,.5,0), oma=c(0,1,0,0))
  for (j in w)
  { su = cl[[j]];
    plotSurv(su, et4, ylab=toupper(j), xlab="Time (yrs)",
      legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET4"], "black"))
    #plotSurv(su, a, ylab=toupper(j), xlab="Time (yrs)",
    #  legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=2:1)
    plotSurv(su, b, subset=b!="Others", ylab=toupper(j), xlab="Time (yrs)",
      legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET4"], palette()[2]))   
    plotSurv(su, ifelse(e==8, "ET8", "Others"), ylab=toupper(j), xlab="Time (yrs)",
      legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["ET8"], "black"))   
  }
  dev.off()
}

for (i in c("ST", names(ds)))
{ if (i=="ST") { nm='ST'; cl=cli; e=ecot; cl$iDFS=NULL; }
  else
  { nm = c(METABRIC="MB", `Spatial TNBC`="Bulk", `SCAN-B`="SCANB")[i]
    z = ds[[i]]; cl = z$cli; e=z$ecot[,1]; ;
  }
  pdf(paste0(figDir, "kmET.",nm, ".pdf"), height=10, width=10)
  par(mfrow=c(2,2), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0));
  for (what in colnames(cl)[sapply(cl, is, "Surv")])
  { plotSurv(cl[[what]], paste0("ET", e), col=ET.colors, addTable=TRUE, legendPos='bottomleft',
      ppos=.1, ylab=what, xlab="Time (yr)")
  }
  dev.off();
}

for (what in c("RFS", "iBCFS", "DRFS", "OS")) #colnames(cli)[sapply(cli, is, "Surv")])
{ pdf(paste0(figDir, "survEco1vAll.", what, ".pdf"), width=12, height=6);
  par(mfrow=c(2,5), mar=c(2.5,3,1.5,.5), mgp=c(1.5,.5,0))
  for (i in 1:max(ecot))
  { plotSurv(cli[[what]], ifelse(ecot==i, paste0("ET", i), "Other"), xlab="Time (yrs)", ylab=what,
      ppos=.1, legendPos=NULL, addTable=TRUE, col=c(ET.colors[i], 'black'));
  }
  dev.off();
}

for (nm in names(ds))
{ nm2 = c(METABRIC="MB", `Spatial TNBC`="Bulk", `SCAN-B`="SCANB")[nm];
  cl = ds[[nm]]$cli; e = ds[[nm]]$ecot[,1];
  for (what in colnames(cl)[sapply(cl, is, "Surv")])
  { pdf(paste0(figDir, "survEco1vAll.", nm2, ".", what, ".pdf"), width=12, height=6);
    par(mfrow=c(2,5), mar=c(2.5,3,1.5,.5), mgp=c(1.5,.5,0))
    for (i in 1:max(ecot))
    { plotSurv(cl[, what], ifelse(e==i, paste0("ET", i), "Other"),
        xlab="Time (yrs)", ylab=what, ppos=.1,
        legendPos=NULL, addTable=TRUE, col=c(ET.colors[i], 'black'));
    }
    dev.off();
  }
}

## Recovery of ET
###################
b = ds[["Spatial TNBC"]]$ecot[,1]; names(b) = rownames(ds[["Spatial TNBC"]]$cli); b = b[rownames(cli)]
sum(diag(table(ecot, b)))
t = matrix(as.character(table(ecot, factor(b, levels=1:max(ecot)))), ncol=max(ecot));
t[t=="0"] = "''"; colnames(t) = rownames(t) = paste0("ET", 1:max(ecot));
diag(t) = paste0("bold('", diag(t), "')");
t[t=="bold('''')"] = "''"; 
cairo_pdf(paste0(figDir, "confMatEcot.pdf"), width=4, height=4)
par(mar=c(0,1.5,1.5,0))
plotTable(t, parseCells=TRUE,
  tableLabels=c("Real ecotype", "Recovered ecotype"))
dev.off();

## Eco comp
############
for (version in c("MB", "Bulk", "ScanB"))
{ wg = intersect(rownames(PB), intersect(geneList$symbol,
    rownames(list(MB=ds$METABRIC$dn, Bulk=ds[["Spatial TNBC"]]$dn, ScanB=ds$`SCAN-B`$dn)[[version]])))
  
  for (what in c("Sigs", "xCell", "Genes", "Annot", "Cli"))
  { if (what %in% c("Annot", "Cli") && version!="Bulk") { next; }
    cc = list(Sigs=csPB, Genes=t(PB[wg,]), xCell=xcPB, Annot=cli$annotClean[,1:11], Cli=cliShort)[[what]];
    d = ds[[c(MB="METABRIC", Bulk="Spatial TNBC", ScanB="SCAN-B")[version]]];
    ccM = list(Sigs=d$cs, Genes=t(d$dn[wg,]), xCell=d$xc)[[what]];
    e2 = factor(d$ecot[,1]);
    
    levels(e2) = paste0("ET", levels(e2));
    e = factor(ecot); levels(e) = paste0("ET", levels(e));
    
    if (!is.null(ccM)) { w = intersect(colnames(cc), colnames(ccM)); cc=cc[,w]; ccM=ccM[,w]; }
  
    co=NULL;
    if (what=="Genes")
    { o = order(-match(colnames(cc), geneList$symbol)); cc = cc[,o]; ccM = ccM[,o];
      co = colGenes[colnames(cc)];
      names(co) = colnames(cc) = colnames(ccM) = geneList[colnames(cc), "List.of.interesting.genes"]
    }
    if (what=="xCell")
    { cc = cc[,!grepl("Score", colnames(cc))];
      o = order(-match(colnames(cc), names(colXct))); cc = cc[,o]; ccM=ccM[,o];
      co = colXct[colnames(cc)]; names(co) = colnames(cc);
    }
    if (what=="Sigs")
    { cc = cc[,intersect(names(sigInfo), colnames(cc))];
      ccM = ccM[,intersect(names(sigInfo), colnames(ccM))];
      co = colSig[sigInfo[colnames(cc)]]; names(co) = colnames(cc);
    }
  
    p1 = calcP(cc, ecot); 
    p = rowMins(p1$fdr);  
    if (!is.null(ccM)) { p2 = calcP(ccM, e2); p = pmin(p, rowMins(p2$fdr))  }
    
    w = colnames(cc)[which(p<=.05)];
    if (length(w)==0) { next; }
    
    z = p.adjust(p1$p[w,], method='fdr');
    cutOff = max(z[p1$fdr[w,]<=.05])+1e-10;
    
    if (is.null(ccM))
    { cairo_pdf(paste0(figDir, "EcoComp.", what, ".pdf"),
        height=(length(w)+7)/7, width=6)
    } else 
    { cairo_pdf(paste0(figDir, "EcoComp.", what, ".", version, ".pdf"),
        height=(length(w)+7)/7, width=6)
    }
    par(mar=c(3,18,.5,.5), mgp=c(1.5,.5,0), cex=.5, omi=c(0,0,0,1));
    dotPlot(cc[,w,drop=FALSE], e, col.lbl=co[w], pch=c("full", "left")[(!is.null(ccM))+1], maxP=cutOff,
      oma=NULL);
    if (!is.null(ccM))
    { z = p.adjust(p2$p[w,], method='fdr');
      cutOff = max(z[p2$fdr[w,]<=.05])+1e-10;
      dotPlot(ccM[,w,drop=FALSE], e2, pch='right', add=TRUE, maxP=cutOff, oma=NULL);
      legendDotPlot(x=par('usr')[2]+strwidth("M")*6, y=par('usr')[4]-strheight("J")*1, interline=2,
        double=TRUE, doubleAnnot=c("ST PB", c(ScanB="SCAN-B", MB="METABRIC", Bulk="ST Bulk")[version]))
    } else { legendDotPlot(x=par('usr')[2]+strwidth("M")*6, y=par('usr')[4]-strheight("J")*3, interline=2) }

    dev.off()
  }
}

# ET by ... heatmaps
#####################
sp = split(seq_along(ecot), ecot);
ecot2 = factor(paste0("ET", ecot));
for (i in c("xc", "sig", "genes"))
{ if (i=="xc") { x = xcPB; annots=colXct; colAnnots=colXc; }
  if (i=="sig") { x = csPB; annots=sigInfo; colAnnots=colSig; }
  if (i=="genes")
  { x = t(PB[intersect(geneList$symbol, rownames(PB)),]); 
    colnames(x) = geneList[colnames(x), "List.of.interesting.genes"];
    annots=geneList$Role.of.genes; names(annots) = geneList$List.of.interesting.genes;
    u = unique(geneList$Role.of.genes);
    colAnnots=palette.colors(palette="Polychrome 36")[-2][seq_along(u)]; names(colAnnots)=u;
  }
  x = x[, rev(names(annots)[names(annots)%in%colnames(x)])];
  if (all(annots%in%names(colAnnots))) { z=colAnnots[annots]; names(z)=names(annots); annots=z; }
  
  p1 = calcP(x, ecot); 
  p = rowMins(p1$fdr);  
  
  w = colnames(x)[which(p<=.05)];
  x = x[,w];
  
  z = p.adjust(p1$p[w,], method='fdr');
  cutOff = max(z[p1$fdr[w,]<=.05])+1e-10;
  
  cairo_pdf(paste0(figDir, "ET.HM.", i, ".pdf"),
      height=1.2*((length(w)+7)/7+1), width=5+8+max(nchar(w))/10)
  par(mar=c(6,max(nchar(w)+4)/2.5,.5,.5), mgp=c(1.5,.5,0));
  layout(matrix(1:2, nrow=1), width=c(.4,.6));
  
  dotPlot(x, ecot2, col.lbl=annots[colnames(x)], maxP=cutOff, oma=NULL, cex.pch=.45);
  legendDotPlot(.5,-3*strheight("M"), horizontal=TRUE, cex.pch=.45, interline=.8)

  if (i=="xc") { x = log(x+.01);  }

  a = t(x[unlist(sp),]);
  a = a-rowMins(a); a=a/rowMaxs(a);

  par(mar=c(6,max(nchar(rownames(a))+4)/2.5,.5,.5));

  image(t(a), xaxt='n', yaxt='n', col=colorRampPalette(c('blue', 'yellow'))(100), xaxs='i', yaxs='i')
  mtext(rownames(a), at=((1:nrow(a))-1)/(nrow(a)-1), side=2, line=.5,
      col=annots[rownames(a)], las=2)
  ww = cumsum(sapply(sp, length));
  abline(v = (ww[-9]-.5)/(ncol(a)-1), col='white', lwd=2)
  mtext(paste0("ET", 1:9), side=1, at=((ww+c(0, ww[-9]))/2-.5)/(ncol(a)-1),line=0.5)
  plotScale(colorRampPalette(c('blue', 'yellow'))(100), c('Low', 'High'), atV=c(0,1), horizontal=TRUE,
    width=20, height=2, posx=0, posy=par('usr')[3]-6*strheight("M"))
  dev.off();
}

## ET single plot
####################
toDisp = read.xlsx(paste0(dataDir, "misc/ET-characteristics.xlsx"), 2, header=FALSE)[,1:2]
toDisp[,2] = sub(" +$", "", toDisp[,2])
toDisp[,2] = sub(" gene$", "", toDisp[,2]);
toDisp[,2] = sub("Myc", "MYC", toDisp[,2]);
toDisp[which(toDisp[,2]=="VCpredTN"),2] = "VCpred TN";
toDisp=toDisp[rowSums(!is.na(toDisp))>0,];
toDisp[which(toDisp[,2]=="GENE70"),2]="GGI";
at = (1:nrow(toDisp))[!is.na(rev(toDisp[,2]))];
colDisp = c(Immune="#68339A", Stroma="#2F6EBA", Oncogenic="#AF2318", `Stress response`="#4EAD5B",
  Metabolism="#F8DA78");
ETnames = read.xlsx(paste0(dataDir, "misc/ET-characteristics.xlsx"), 1, header=TRUE)[,1]
ETnames = ETnames[!is.na(ETnames)];
w = which(!is.na(toDisp[,1])); cl = rep(toDisp[w,1],c(w[-1], nrow(toDisp)+1)-w-1)

wg = intersect(rownames(PB), intersect(geneList$symbol, rownames(bulk)))
x = t(PB[wg,]); colnames(x) = geneList[wg,"List.of.interesting.genes"];
cc = cbind(csPB, x, xcPB[,1:48])
pi = apply(cc, 2, function(i) tapply(1:nrow(cc), ecot, function(s) wilcox.test(i[s], i[-s])$p.value))
pi[] = p.adjust(pi, method='fdr')
w = rev(toDisp[!is.na(toDisp[,2]),2])
z = cc[,w]; a = colnames(z) %in% geneList[wg,"List.of.interesting.genes"];
#colnames(z)[a] = paste(colnames(z)[a], "(gene)");
lbls = paste0("'", colnames(z), "'"); lbls[a] = paste0("italic(", lbls[a], ")"); 
#lbls=sub("^(.+) \\(gene\\)$", "paste(italic('\\1'), ' gene')", colnames(z));
#ww=grepl("^paste", lbls); lbls[!ww]=paste0("'", lbls[!ww], "'");
lbls = parse(text=lbls);
wxc = colnames(z) %in% colnames(xcPB);
colLbls = c('grey45', 'black', 'darkblue'); names(colLbls) = c("Cell types", "Signatures", "Single genes");
#colLbls=colLbls[c(2,1,3)];
col = colLbls["Signatures"]; col[a] = colLbls["Single genes"]; col[wxc] = colLbls["Cell types"];
o = order(factor(rev(cl), levels=unique(rev(cl))), factor(col,levels=colLbls[c(1,3,2)]));
lbls=lbls[o]; z=z[,o]; col=col[o];

cairo_pdf(paste0(figDir, "EcoSingle.horiz.pdf"), width=8.6, height=2.5)
par(mar=c(1,18,9,3), oma=c(0,0,0,4), mgp=c(1.5,.5,0), cex=.5, xpd=NA);
dotPlot(z, factor(ETnames[ecot], levels=rev(ETnames)), toDisp=t(pi[nrow(pi):1,w])<.05, at=max(at)-at+1,
  horizontal=TRUE, axPos=2, srt=30, inMa=1, lbls=lbls, col.lbl=col, oma=NULL)
wi = (which(!is.na(toDisp[,1]))-1)[-1];
rect(wi-.9, rep(9.7,4), wi+.9, rep(10.5, 4), col='white', xpd=NA, border=NA);
wi = c(0, wi, nrow(toDisp));
mtext(txt<-toDisp$X1[!is.na(toDisp$X1)], side=1, line=-.3, at=(wi[-1]+wi[-length(wi)])/2, cex=.5)
for (i in 1:(length(wi)-1))
{ lines(wi[(0:1)+i]+c(.5,-.5), c(.3,.3), xpd=NA, col=colDisp[txt[i]], lwd=2)
}
legendDotPlot(par('usr')[2]+2, par('usr')[4]+1, interline=2)
legend(par('usr')[2]+2, par('usr')[4]-6, names(colLbls)[c(2,3,1)], fill=colLbls[c(2,3,1)], bty='n', xjust=.5,
  text.col=colLbls[c(2,3,1)], border=NA)
dev.off();