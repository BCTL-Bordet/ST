figDir = "~/Documents/BCTL/ST/TNBC/figsSP/";
figDataDir = "~/Documents/BCTL/ST/TNBC/figsDataSP/";

library(vioplot); library(colorspace);
library(alluvial);
library(NbClust)
library(scales)
colDisp = c(Immune="#68339A", Stroma="#2F6EBA", Oncogenic="#AF2318", `Stress response`="#4EAD5B",
  Metabolism="#F8DA78");

# Note: originally spatial Archetypes were called ecotypes, hence the "ET" and "ecot" all around
# the code and the file names

# Figure files
#
# x Fig 3c - morphoDistUp.pdf
# x Fig 4b - alluvialDec.2.pdf
# x Fig 4c - tumorCompartment.pdf
# x Fig 4d - stromaCompartment.pdf
# x Fig 4e - TvsS.Mbp.hm.pdf
# Fig 4f - survTSdec.final.1.pdf
# x Fig 5b - tlsViolin.pdf
# x Fig 5c - tlsGO.VSlympho.blueOnly.pdf
# x Fig 5d - tlsGene.pdf
# x Fig 6a-d - tlsAllCliSmall.pdf
# Fig 6e - controlTLSimmunoT.pfs.pdf
# Fig 6f - controlTLSby_immunoT.pfs.pdf
# x Fig 7c - barByClust.pdf
# x Fig 7d - tblsMC2.pdf
# x Fig 7e - forestMC.PB.pdf
# x Fig 7f - forestMC.All3.pdf
# x Fig 8a - heatmapMC.pdf
# x Fig 8b - barplotMCvsET.pdf
# x Fig 8c - barplotETs.pdf
# x Fig 8d - EcoSingle.horiz.pdf
# x Fig 8e - hmTargets.pdf
# x Fig 8f - forestET.All3.pdf
# x Fig 8g - survEcotypeNew2.pdf
#
# Supplementary figures
# x Fig S2a - aucs.pdf
# x Fig S3a - TvsS.bar.dec.Tumor.*.horiz.pdf
# x Fig S3b - TvsS.bar.dec.Stroma.*.horiz.pdf
# Fig S4 -
# x Fig S5a - tlsXcl.pdf
# x Fig S5b - tlsGO.VSlympho.pdf
# x Fig S5c - tlsGO.VSrest.pdf
# x Fig S6a - aucTLS2.pdf
# x Fig S6b - aucTLS.pdf
# x Fig S6c - compTLSsigs.pdf
# x Fig S7a - tlsIspy.pdf
# Fig S7b - 
# Fig S7c - tlsAllCliRest.pdf
# Fig S7d - fpTLSimmuno.pdf
# x Fig S8a - controlTLS_TNBC.pdf
# x Fig S8b - controlTLSby_TNBC.pdf
# x Fig S9a - controlTLSispy.pdf
# x Fig S9b - controlTLSbyispy.pdf
# x Fig S9cd - controlTLSimmunoT.all.pdf
# x Fig S9ef - controlTLSby_immunoT.all.pdf
# x Fig S11a - meanAUC.pdf
# x Fig S11bc - MCcomp.Sigs.pdf
# x Fig S12 - MCcomp.Genes.pdf
# x Fig S13 - MCcomp.xCell.pdf
# x Fig S14 - forestMC.PB.pdf
# x Fig S15 - forestMC.Bulk.pdf
# x Fig S16 - forestMC.MB.pdf
# x Fig S17 - forestMC.SCANB.pdf
# x Fig S18 - forestMC.All3.pdf
# x Fig S19a - N_ET.pdf
# x Fig S19b - barEcot.pdf
# x Fig S19c - EcoComp.Annot.pdf
# x Fig S19d - EcoComp.Cli.pdf
# x Fig S20ac - ET.HM.xc.pdf
# x Fig S20bd - ET.HM.sig.pdf
# x Fig S21ab - ET.HM.genes.pdf
# x Fig S22a - EcoComp.Sigs.Bulk.pdf
# x Fig S22b - EcoComp.Sigs.MB.pdf
# x Fig S23 - EcoComp.Sigs.ScanB.pdf
# x Fig S24 - EcoComp.xCell.Bulk.pdf
# x Fig S25 - EcoComp.xCell.MB.pdf
# x Fig S26 - EcoComp.xCell.ScanB.pdf
# x Fig S27a - EcoComp.Genes.Bulk.pdf
# x Fig S27b - EcoComp.Genes.MB.pdf
# x Fig S28 - EcoComp.Genes.ScanB.pdf
# x Fig S29ab - forestET.PB.pdf
# x Fig S30ab - forestET.Bulk.pdf
# x Fig S31ab - forestET.MB.pdf
# x Fig S32ab - forestET.SCANB.pdf
# x Fig S33ab - forestMC.All3.pdf

# Fig S4a/b/c - survTSdec.final.2.pdf
# Fig S4d - TvsS.bar.dec.Mpb.Cli.horiz.pdf, TvsS.bar.dec.Mpb.Annot.horiz.pdf

# Fig S5b - sigCharoentongInTLS.FP.pdf
# Fig S5c - sigCharoentongInTLS.pdf

# Fig S6c - tlsAllCliRest.pdf

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
par(mfcol=c(1,2), mar=c(3,3,1,.5), mgp=c(1.5,.5,0), oma=c(0,0,0,12), cex=1.3);
for (j in 1:2)
{ if (j==1) { yl = c(0,100); } else { yl = c(95, 100); }
  h = barplot(t(n3), col=co, ylab="% pixels", border=NA, ylim=yl, xpd=FALSE)
  par(xpd=TRUE)
  if (j==1)
  { lines(c(h[length(h)]+.5, par('usr')[2]+2.5*strwidth("M")), c(100, 100), xpd=NA, lwd=2)
    lines(c(h[length(h)]+.5, par('usr')[2]+2.5*strwidth("M")), c(95, 0), xpd=NA, lwd=2)
  }
 
  if (j==2) { legend('topright', inset=c(-max(strwidth(colnames(n2), units="figure"))*1.4,0),
    legend=rev(colnames(n2)), fill=rev(co), xpd=NA, bty='n')} 
}
dev.off();

x = data.frame(ID=rownames(n2), n2, subtype=cli$barPB, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "morphoDistUp.xlsx"), row.names=FALSE);

# Full version
pdf(paste0(figDir, "morphoDistDetails.pdf"), height=7.5, width=12)
par(mfrow=c(3,4), mar=c(2.5,3,.5,.5), mgp=c(1.5,.5,0))
for (i in colnames(n2))
{ ylim = range(n2[,i]); if (i=="Stroma") { ylim[2]=ylim[2]+5; }
  myBoxplot(n2[,i], cli$barPB, ylab=paste(i, " (% pixels)"), subLine=-1.5, colPoints=colBar[cli$barPB],
    subSide=3, ylim=ylim)#, Plog=1)
}
dev.off()

# Other datasets
###################
for (i in c("Tumor", "Stroma"))
{ x=cli[[paste0("patches.", i)]][,c("N", "Np", "evenness", "N0.5")];
  x[,2] = x[,2]/x[,1];
  colnames(x) = c("Number", "Size", "Evenness", "50% surface recovery");
  x = data.frame(ID=rownames(cli), Subtype=cli$barPB, x, check.names=FALSE);
  write.xlsx2(x, file=paste0(figDataDir, "morphoPatch.", i, ".xlsx"), row.names=FALSE)
}

x = data.frame(ID=rownames(cli), cli$annotations, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "morphoTotal.xlsx"), row.names=FALSE)

x = csAnT;
y = cli$patches.Tumor; y = y[,"Np"]/y[,"N"];
x = data.frame(ID=rownames(cli), `Mean tumor patch size`=y, x, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "patchSizeVsSig.Tumor.xlsx"), row.names=FALSE)

x = csAnS;
y = cli$patches.Stroma; y = y[,"Np"]/y[,"N"];
x = data.frame(ID=rownames(cli), `Mean stroma patch size`=y, x, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "patchSizeVsSig.Stroma.xlsx"), row.names=FALSE)

#####################
## Tumor / stroma PB
#####################
load(paste0(dataDir, 'classification/loosReg.RData'));
a = c("Tumor", "Stroma", "Necrosis", "Fat tissue", "Vessels", "in situ", "Lymphoid nodule",
  "Lymphocyte", "Lactiferous duct");
nm2 = c(`in situ`="In situ", `Lymphoid nodule`="Tertiary lymphoid structures",
  Lymphocyte="Lymphocytes", `Lactiferous duct`="Lactiferous ducts");
pdf(paste0(figDir, "aucs.pdf"), height=2*3, width=2*3)
par(mfrow=c(3,3), mar=c(3,3,1.5,.5), mgp=c(1.5,.5,0))
for (i in a)
{ r = roc(n[,i]>.25, loos[[i]], direction="<", plot=TRUE, mar=par("mar"), mgp=c(1.5,.5,0),
    print.auc=FALSE, print.auc.cex=1.3, smooth=FALSE, print.auc.y=.2, print.auc.x=.6, ci=TRUE)
  mtext(paste("AUC:", format(r$auc, digits=3)), side=1, line=-2.5, adj=0, cex=.75, at=.7);
  mtext(paste("CI:", format(r$ci[1], digits=3), "to", format(r$ci[3], digits=3)), side=1,
    line=-1.25, adj=0, cex=.75, at=.7);
  if (any(names(nm2)==i)) { title(main=nm2[i], xpd=NA); } else { title(main=i); }
}
dev.off();

x = data.frame(n[,a]>0,do.call(cbind, loos[a]), check.names=FALSE);
w = colnames(x) %in% names(nm2); colnames(x)[w] = nm2[colnames(x)[w]];
colnames(x) = paste(colnames(x), rep(c("Presence (TRUE/FALSE)", "LPO regression"), each=length(a)))
write.xlsx2(x, file=paste0(figDataDir, "aucs.xlsx"), row.names=FALSE);

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
x = data.frame(Tumor=f(cli$barT.an), PB=f(cli$barPB), Stroma=f(cli$barS.an),
    TIME=factor(cli$TIME_classes.bypathologist, levels=rev(c('ID', 'MR', 'SR', 'FI'))))
pdf(paste0(figDir, "alluvialDec.2.pdf"));
par(mfrow=c(1,1), mar=c(3,3,0.5,.5), mgp=c(1.5,.5,0), xpd=NA)
#t = plyr::count(data.frame(Tum=f(tT), Tumor=f(cli$barTS.T), Full=f(cli$barPB), Stroma=f(cli$barTS.S),
#  Str=f(tS)));
t = plyr::count(x);
t = t[complete.cases(t),];
alluvial(t[,1:4], freq=t$freq, col=colBar[as.character(t$PB)], blockCol=c(colBar,  colTIME),
  draw_ticks=FALSE)
dev.off();

x = cbind(ID=rownames(cli), x);
write.xlsx2(x, file=paste0(figDataDir, "alluvialDec.2.xlsx"), row.names=FALSE)

pdf(paste0(figDir, "alluvialDec.3.pdf"));
par(mfrow=c(1,1), mar=c(3,3,0.5,.5), mgp=c(1.5,.5,0), xpd=NA)
t = plyr::count(data.frame(Tumor=f(cli$barT.an), PB=f(cli$barPB), Stroma=f(cli$barS.an),
  TIME=factor(cli$TIME_classes.bypathologist, levels=rev(c('ID', 'MR', 'SR', 'FI')))));
t = t[complete.cases(t),];
alluvial(t[,c(2,1,3,4)], freq=t$freq, col=colBar[as.character(t$PB)], blockCol=c(colBar,  colTIME),
  draw_ticks=FALSE)
dev.off();

pdf(paste0(figDir, "alluvialDec.4.pdf"));
par(mfrow=c(1,1), mar=c(3,3,0.5,.5), mgp=c(1.5,.5,0), xpd=NA)
t = plyr::count(data.frame(Tumor=f(cli$barT.an), PB=f(cli$barPB), Stroma=f(cli$barS.an),
  TIME=factor(cli$TIME_classes.bypathologist, levels=rev(c('ID', 'MR', 'SR', 'FI')))));
t = t[complete.cases(t),];
alluvial(t[,c(2,1,4)], freq=t$freq, col=colBar[as.character(t$PB)], blockCol=c(colBar,  colTIME),
  draw_ticks=FALSE)
dev.off();

nfo = read.xlsx(paste0(dataDir, "misc/Stroma and tumor distribution of molecular feature.xlsx"), 1)[,1:4];
#nfo = lapply(1:3, function(i) read.xlsx(paste0(dataDir, "misc/List of signatures, genes of interest, cell type enrichment analysis computed in the study.xlsx"), sheetIndex=i))
nfo = tapply(1:nrow(nfo), nfo[,1], function(i) nfo[i,2:4]);
nfo = nfo[c("Signature", "Sinngle_gene", "Cell_type_xcell")];
names(nfo) = c("Sigs", "Genes", "xCell");
for (i in 1:3) { nfo[[i]][,2] = sub("^ *([^ ].+[^ ]) *$", "\\1", nfo[[i]][,2]); }
nfo[[1]][nfo[[1]][,2]=="VCpred_TN",2]="VCpred TN";
nfo[[1]][,2] = sub("^TNFA signaling.+", "TNFA signaling via NF-kB", nfo[[1]][,2]);
nfo[[2]][nfo[[2]][,2]=="mTOR",2]="MTOR";
nfo[[1]][nfo[[1]][,2]=="TGFbeta-myCAF",2]="TGFβ-myCAF";
nfo[[1]][nfo[[1]][,2]=="IFNgamma-iCAF",2]="IFNγ-iCAF";
nfo[[1]][nfo[[1]][,2]=="Normal fibroblast",2]="Normal Fibroblast";
for (i in 1:3) { rownames(nfo[[i]]) = nfo[[i]][,2]; }
rownames(nfo[[2]]) = rownames(geneList)[match(rownames(nfo[[2]]),
  geneList$List.of.interesting.genes)];


f = function(x) factor(x, levels=names(colBar));
bl = paste(cli$barPB, cli$TIME_classes.bypathologist);
bl[!(bl %in% c("BL SR", "IM FI", "IM SR"))] = NA; bl[cli$barT.an!="BL"] = NA;
bl = factor(bl, levels=c("BL SR", "IM SR", "IM FI"))
me = paste(cli$barPB, cli$barS.an);
me[!(me %in% c("MSL MSL", "M MSL", "M M"))] = NA; me[cli$barT.an!="M"] = NA;
me2 = paste("Stroma", cli$barS.an); me2[cli$barPB!="M" | !(cli$barS.an %in% c("M", "MSL"))] = NA;
rownames(geneList) = geneList$symbol;
for (what in c("Sigs", "xCell", "Genes", "Annot", "Cli"))
{ for (version in c("StromaByTumor", "Tumor", "Stroma", "BL", "M", "Mpb")[2:3])
  { v2 = c(StromaByTumor="Stroma", Tumor2="Tumor", Tumor="Tumor", Stroma="Stroma", Stroma2="Stroma",
      BL="Stroma", M.Tumor="Tumor", M="Stroma", Mpb="Stroma")[version]
    cc = list(Tumor=list(Sigs=csAnT, SigsH=csAnT, xCell=xcAnT, Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol),]), Annot=cli$annotClean, Cli=ccn),
      Stroma=list(Sigs=csAnS, SigsH=csAnS, xCell=xcAnS, Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol),]), Annot=cli$annotClean, Cli=ccn))[[v2]][[what]]
    cmps = list(StromaByTumor=cli$barT.an, Tumor2=f(cli$barT.an), Tumor=f(cli$barT.an), Stroma=f(cli$barS.an), Stroma2=f(cli$barS.an), BL=bl,
      M=me, Mpb=me2)[[version]]
    if (version%in%c("StromaByTumor", "Tumor")) { cmps[cmps=="IM"]=NA; }
    if (version=="Stroma") { cmps[cmps=="BL"]=NA; }
    #if (what=="Genes" && v2=="Stroma")
    #{ cc = cc[,!(geneList[colnames(cc), "Role.of.genes"] %in% c("Methylation", "HRD", "Cycling", "MMR"))]
    #}
    
    n = nfo[[what]];
    if (!is.null(n))
    { cc = cc[,colnames(cc) %in% (rownames(n)[n[,3] %in% list(Tumor=c("T", "T and S"), Stroma=c("S", "T and S"))[[v2]]]), drop=FALSE]
    }
    
    bubullePlot(cc, cmps, what=what,
      fileName=paste0(figDir, "TvsS.bar.dec.",version, ".pdf"), dataDir=figDataDir, dataNfo=paste(version,"subtype"))
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
{ a = read.xlsx(paste0(dataDir, "misc/Molecular features selection.xlsx"), i);
  a = a[,!is.na(a[1,])];
  a = a[!is.na(a$List_selected),]; a$List_selected = sub("^ *([^ ].+[^ ]) *$", "\\1", a$List_selected)
  a = a[order(factor(a$Main_classes, levels=c("Immune", "Stress response", "Stroma", "Pathway",
    "Metabolism", "Proliferation"))),]
  return(a);
}); names(aa) = c("T", "S");
col = c(colDisp, Proliferation=colDisp[["Oncogenic"]], Pathway=MC.colors[[1]]);

for (what in c("T", "S"))
{ a = aa[[what]];
  if (what=="T")
  { x = list(Sigs=csAnT, xCell=xcAnT[,names(colXct)], Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol), ]))
    cmps = f(cli$barT.an); cmps[cmps=="IM"]=NA; 
  } else
  { x = list(Sigs=csAnS, xCell=xcAnS[,names(colXct)], Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol), ]))
    cmps = f(cli$barS.an); cmps[cmps=="BL"]=NA;
  }
  #x = lapply(names(x), function(i) x[[i]][,grepl(what, nfo[[i]][colnames(x[[i]]), 3]),drop=FALSE]) # Tumor/stroma features only
  colnames(x[[3]]) = geneList[colnames(x[[3]]), "List.of.interesting.genes"]
  cc = do.call(cbind, x);

  cmps=factor(cmps);
  
  x2 = data.frame(ID=rownames(cli), cmps, cc, check.names=FALSE);
  colnames(x2)[2] = paste(c(T="Tumor", S="Stroma")[what], "subtype");
  write.xlsx2(x2, file=paste0(figDataDir, c(T="tumor", S="stroma")[what], "Compartment.xlsx"),
    row.names=FALSE)

  p1 = calcP(cc, cmps); 

  ww = a$List_selected; ww[ww=="Macrophages M1"]="Macrophages";
  lbls=paste0("'", ww, "'");
  w = a$Feature_type=="Single genes"; lbls[w] = paste0('italic(', lbls[w], ')');
  lbls = parse(text=lbls);
  colLbls = c('grey45', 'black', 'darkblue'); names(colLbls) = c("Cell types", "Signatures", "Single genes");

  z = p.adjust(p1$p[ww,], method='fdr');
  cutOff = max(z[p1$fdr[ww,]<.05])+1e-10;

  at = 1:nrow(a);
  for (i in which(!duplicated(a$Main_classes))[-1]) { at[i:length(at)] = at[i:length(at)]+1; }

  leftMar = 1+max(nchar(ww))/2; 
  sl = 0.06*max(nchar(colnames(p1$p))+8)/2; sr = 0.06*(max(nchar(ww))+3)*.6; sb=(0.06*leftMar*.7+2)/3;
  if (what=="T") { sr = sr+1; }
  cairo_pdf(paste0(figDir, c(T="tumor", S="stroma")[what], "Compartment.pdf"),
        height=(ncol(p1$p))*.23 + sb + .2, width=.2*(length(ww)) + sl + sr)
  par(mai=c(.2, sl, sb, sr), mgp=c(1.5,.5,0)+.2, cex=.6);
  dotPlot(cc[,ww,drop=FALSE], cmps, maxP=cutOff, legend=FALSE, cex.pch=1, srt=30, horizontal=TRUE, axPos = 3, lbls=lbls,
    col.lbl=colLbls[a$Feature_type], at=at);
  wi = at[which(!duplicated(a$Main_classes))[-1]]-1
  ba = par('usr')[4]+.2*strheight("M");
  rect(wi-.95, rep(ba-.05,length(wi)), wi+.95, rep(ba+.1, length(wi)), col='white', xpd=NA, border=NA);
  wi2 = at[which(a$Feature_type[-1] != a$Feature_type[-nrow(a)])]+1;
  rect(wi2-.7, rep(ba-.05,length(wi2)), wi2-.3, rep(ba+.1, length(wi2)), col='white', xpd=NA, border=NA);
  wi = c(0, wi, at[nrow(a)]+1);
  mtext(txt<-a$Main_classes[wi[-length(wi)]+1], side=1, line=.4, at=(wi[-1]+wi[-length(wi)])/2, cex=.6)
  for (i in 1:(length(wi)-1))
  { lines(wi[(0:1)+i]+c(.7,-.7), c(.3,.3), xpd=NA, col=col[txt[i]], lwd=2)
  }
  if (what == "T")
  { xs = max(at)+3; ys = par('usr')[4]+2.5;
    legendDotPlot(xs, ys, horizontal=FALSE, interline=1.6, significantLabel=NULL)
    li = strheight("M", cex=par('cex'))*1.5*1.6;
    em = strwidth("M", cex=par('cex'));
    xs = xs+5;
    xs2 = rep(xs, 2);
    ys2 = ys-li*1.5*c(0,2.5)
    text(xs2, ys2, c("FDR < 5%", "FDR ≥ 5%"), xpd=NA, font=2, adj=c(.5, 1))
    points(xs2[1]+c(-1.5*em, 1.5*em)/par('cex'), rep(ys2[1]-li*2,2), pch=21, col='black', lwd=.5,
      cex=3, xpd=NA, bg = c("#0072B2", "#F27052"));
    points(xs2[2]+c(-1.5*em, 1.5*em)/par('cex'), rep(ys2[2]-li*2,2), pch=16,
      col=lighten(c("#0072B2", "#F27052"), .6), cex=3, xpd=NA);

    legend(xs,ys-li*1.5*5, legend=sub("s$", "", names(colLbls)), fill=colLbls, bty='n', xjust=.5, border=NA,
      title="Molecular\nfeature type", title.font=2, xpd=NA);
    #le = le$text$y; le = le[1]-1.6*diff(le)[1]
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

x = data.frame(ID=rownames(cli), Class=a, DRFS.time=cli$DRFS[,1], DRFS.event=cli$DRFS[,2])
x = x[!is.na(x[,2]),];
x$Class = ifelse(grepl("MSL", x$Class), "M-MSLstroma", "M-Mstroma")
write.xlsx2(x, file=paste0(figDataDir, "survTSdec.final.xlsx"), row.names=FALSE);

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
SA = lapply(split(rownames(y), me2), function(i) as.character(sort(as.numeric(i))));
y = y[unlist(SA),];
rq = colQuantiles(y, probs=c(.05, .95));
y = (y-rep(rq[,1], each=nrow(y)))/rep(rowDiffs(rq), each=nrow(y))
y[y<0] = 0; y[y>1] = 1;

cairo_pdf(paste0(figDir, "TvsS.Mbp.hm.pdf"), height=5)
par(mar=c(3,10,1.5,2.5))
image(y[,ncol(y):1], col=colorRampPalette(c('blue', 'yellow'))(100), xaxt='n', yaxt='n');
wi = seq(0,1,length=ncol(y)); nd = which(!duplicated(rev(wh)))[-1];
mtext(rev(colnames(y)), side=2, at=wi, las=2, line=.5, font=1+2*(rev(wh)=="Gene")) 
mtext(rownames(y), side=1, at=seq(0,1,length=nrow(y)), cex=.5, line=-.3)
abline(v = (length(SA[[1]])-.5)/(nrow(y)-1), col='white', lwd=2)
abline(h = wi[nd]-wi[2]/2, col='white', lwd=2) 
xi = c(length(SA[[1]])*.5, length(SA[[1]])+.5*length(SA[[2]]))/(nrow(y));
xi = xi*diff(par('usr')[1:2])+par('usr')[1]
text(xi, rep(par('usr')[4]+strheight("M"), 2), xpd=NA,
 # parse(text=c("bold(paste('M'['TUMOR'], '-M'['STROMA']))", "bold(paste('M'['TUMOR'], '-MSL'['STROMA']))")), col=1:2)
  parse(text=c("bold(paste('M-M'['STROMA']))", "bold(paste('M-MSL'['STROMA']))")), col=1:2)
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

w =  !is.na(me2);
x2 = data.frame(ID=rownames(x)[w], `Stroma subtype`=me2[w],
  `PB subtype`=cli$barPB[w], x[w,], check.names=FALSE)
write.xlsx2(x2, file=paste0(figDataDir, "TvsS.Mbp.hm.xlsx"), row.names=FALSE)

cmps = StromaByTumor=cli$barT.an
cmps[cmps=="IM"]=NA;
cmps = factor(cmps);
for (what in c("Sigs", "xCell", "Genes"))
{ for (ver in c("StromaByTumor", "Tumor")[2])
  { x = list(stromaByTumor=list(Sigs=csAnS, xCell=xcAnS, Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol),])),
        Tumor=list(Sigs=csAnT, xCell=xcAnT, Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol),])))[[ver]][[what]];
    if (what=="Genes")
    { if (ver=="StromaByTumor") { x = x[,!(geneList[colnames(cc), "Role.of.genes"] %in% c("Methylation", "HRD", "Cycling", "MMR"))] }
      co = factor(geneList[colnames(x), "Role.of.genes"], levels=unique(geneList[, "Role.of.genes"]));
      co = palette.colors(palette="Polychrome 36")[-2][unclass(co)];
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
    rownames(x) = rownames(cli);


    p = moulinette(x, cmps);
    p2 = p.adjust(p, method='fdr');
    y = x[, p2<.05];# wh = what[p2<.05]
    co = co[p2<.05];
    SA = lapply(split(rownames(y), cmps), function(i) as.character(sort(as.numeric(i))));
    y = y[unlist(SA),];
    rq = colQuantiles(y, probs=c(.05, .95));
    y = (y-rep(rq[,1], each=nrow(y)))/rep(rowDiffs(rq), each=nrow(y))
    y[y<0] = 0; y[y>1] = 1;

    cairo_pdf(paste0(figDir, "TvsS.", ver, ".", what, ".hm.pdf"), height=.2+.15*ncol(y))
    par(mar=c(.5,10,1.5,.5))
    image(y[,ncol(y):1], col=colorRampPalette(c('blue', 'yellow'))(100), xaxt='n', yaxt='n');
    wi = seq(0,1,length=ncol(y)); #nd = which(!duplicated(rev(wh)))[-1];
    mtext(rev(colnames(y)), side=2, at=wi, las=2, line=.5, col=rev(co), cex=.7)#, font=1+2*(rev(wh)=="Gene")) 
   #mtext(rownames(y), side=1, at=seq(0,1,length=nrow(y)), cex=.5, line=-.3)
    ni = (cumsum(sapply(SA, length))-.5)/(nrow(y)-1)
    abline(v = ni[-3], col='white', lwd=2)
    xi = (c(0, ni[-3]) + ni)/2;
    #xi = xi*diff(par('usr')[1:2])+par('usr')[1]
    mtext(names(SA), side=3, at=xi, col=colBar[names(SA)])
    dev.off();
  }
}

for (what in c("Sigs", "xCell", "Genes"))
{ for (ver in c("StromaByTumor", "Tumor")[2])
  { cc = list(stromaByTumor=list(Sigs=csAnS, xCell=xcAnS, Genes=t(stromaAn[intersect(rownames(stromaAn), geneList$symbol),])),
        Tumor=list(Sigs=csAnT, xCell=xcAnT, Genes=t(tumorAn[intersect(rownames(tumorAn), geneList$symbol),])))[[ver]][[what]];
    if (what=="Genes" && ver=="StromaByTumor")
    { cc = cc[,!(geneList[colnames(cc), "Role.of.genes"] %in% c("Methylation", "HRD", "Cycling", "MMR"))]
    }
    for (out in c("BL", "M", "LAR"))
    { cmps = StromaByTumor=cli$barT.an
      cmps[cmps%in%c("IM", out)]=NA;
      cmps = factor(cmps);
      bubullePlot(cc, cmps, what=what,
          fileName=paste0(figDir, "TvsS.", ver, ".",paste(levels(cmps), collapse=".vs."), ".pdf"))
    }
  }
}

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

x = do.call(cbind, pi);
rownames(x) = rownames(xCan);
colnames(x) = paste(c("Log FC", "P-value"), rep(names(pi), each=2), "vs TLS")
w = grep("Log FC", colnames(x)); x[,w] = exp(x[,w]); colnames(x)[w]=sub("Log ", "", colnames(x)[w])
x = data.frame(Gene=rownames(x), x, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "tlsGene.xlsx"), row.names=FALSE);

# TLS by xCell
###############
ax = function(i)
{ a=axTicks(i); lbl = paste0("10^-", abs(a)); lbl[a==0] = 1;
  axis(i, at=a, labels=parse(text=lbl));
}
p = psGo$xCell[names(colXct), ]; p=p[!grepl("Score", rownames(p)),];
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

x = data.frame(t(psGo$xCell[names(colXct), -(1:2)]), check.names=FALSE)
x = cbind(comparison=rownames(x), x)
write.xlsx2(x, file=paste0(figDataDir, "tlsXcl.xlsx"), row.names=FALSE);

co1 = colorspace::lighten('#405D92', .5); co2 = lighten('#E86E4D', .15)
p = psGo$xCell[names(colXct),];
p = 2*pmin(p[,"Lymphocyte Higher"], p[,"Lymphocyte Lower"])
q = p.adjust(p, 'fdr');
w = rev(names(q)[q<.001]);
star = symnum(q[w], corr = FALSE, na = FALSE, 
                  cutpoints = c(0, 1e-4, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("****", "***", "**", "*", ".", " "))
#w = rownames(p)[p[,"Lymphocyte Higher"]<1e-4 | p[,"Lymphocyte Lower"]<1e-4];
#w = w[!grepl("Score", w)];
#w = rev(intersect(names(colXct), w))
xl=c(0, .2+max(xcAn[idAn[,"Lymphocyte"],w]));
pdf(paste0(figDir, "tlsViolin.pdf"), height=5);
par(mar=c(3,13,3,.5), mgp=c(1.5,.5,0));
vioplot(xcAn[idAn[cli$hasTLS,"Lymphoid nodule"],w], horizontal=TRUE, side='right', las=1, col=co1,
  xlab="xCell enrichment score", plotCentre = "line", yaxt='n', xaxt='n',
  ylim=xl); axis(1);
text(rep(par('usr')[1]-strwidth("n"), length(w)), seq_along(w), w, xpd=NA, adj=1, co=colXct[w])
vioplot(xcAn[idAn[,"Lymphocyte"],w], horizontal=TRUE, side='left', las=1, col=co2,
  add=TRUE, plotCentre = "line", ylim=xl)
legend(.7, 17, legend=c("TLS", "Lymphocyte compartment"), fill=c(co1,co2), bty='n')
mtext(star, at=seq_along(w), side=4, line=-.1, las=2, adj=1);
dev.off();

w = idCan[,1] %in% c("Lymphoid nodule", "Lymphocyte");
x = data.frame(ID=idCan[w,2], Compartment=
  c(`Lymphoid nodule`="TLS", Lymphocyte="Lymphocytes")[idCan[w,1]], xcAn[w,names(colXct)], check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "tlsViolin.xlsx"), row.names=FALSE)

# tlsAllCli
###############
censor = function(s, t=10) { w = s[,1]>t; s[w,2]=0; s[w,1]=t; s; }
quartiles = function(x) { x = ceiling(4*rank(x)/length(x)); ordered(paste0("Q", x)); }
set.seed(123);
cairo_pdf(paste0(figDir, "tlsAllCli.pdf"), height=6, width=15)
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

colIHC = colsAnnot("IHC"); cliIspy$Type = c(`HR-HER2+`="HER2+", `HR+HER2-`="HR+/HER2-", `HR+HER2+`="HER2+", TN="TNBC")[cliIspy$Receptor.Subtype]
set.seed(123);
#cairo_pdf(paste0(figDir, "tlsAllCliSmall.pdf"), height=10, width=15)
#par(mfrow=c(2,3), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), cex.main=1, cex=1);
cairo_pdf(paste0(figDir, "tlsAllCliSmall2.pdf"), height=5, width=20)
par(mfrow=c(1,4), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), cex.main=1, cex=1);
#plotSurv(censor(cli$DRFS), quartiles(csPB[,'TLS ST']), addTable=TRUE, main="ST TNBC cohort", ylab="DRFS", xlab="Time (yr)",
#  legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, args.legend=list(bty='n', inset=c(.02, 0)), leftTable=5)
#plotSurv(censor(ds$METABRIC$cli$DRFS), quartiles(ds$METABRIC$cs[,"TLS ST"]), addTable=TRUE, main="METABRIC - TNBC cohort", ylab="DRFS",
#  xlab="Time (yr)", legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, args.legend=list(bty='n', inset=c(.02, 0)), leftTable=5)
#plotSurv(censor(ds$`SCAN-B`$cli$DRFS), quartiles(ds$`SCAN-B`$cs[,"TLS ST"]), addTable=TRUE, main="SCAN-B - TNBC cohort", ylab="DRFS",
#  xlab="Time (yr)", legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5, , args.legend=list(bty='n', inset=c(.02, 0)), leftTable=5)
qrt = unlist(lapply(ds, function(i) quartiles(i$cs[,"TLS ST"])))
plotSurv(censor(sus$DRFS), qrt, addTable=TRUE, main="TNBC - all cohorts", ylab="DRFS",
  pval=getCoxP(coxph(censor(sus$DRFS) ~ qrt +strata(stud)) ),
  xlab="Time (yr)", legendPos='bottomleft', legend.title="TLS ST sig", ppos=.2, marTable=1.5,
  args.legend=list(bty='n', inset=c(.02, 0)), leftTable=5)
#myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$wI & cliIspy$arm=="Pembro", ylab="TLS ST sig", subSide=3, subLine=-1.5,
#  main="I-SPY2 - TNBC cohort - Pembrolizumab arm", alpha=.9, addN=TRUE)
myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$arm=="Pembro", ylab="TLS ST sig", subSide=3, subLine=-1.5,
  main="I-SPY2 - Pembrolizumab arm", alpha=.9, addN=TRUE, colPoints=colIHC[cliIspy$Type], ylim=c(min(csI[,"TLS ST"]), 1.2))
legend('topleft', legend=c("TNBC", "HR+/HER2-"), fill=colIHC[c("TNBC", "HR+/HER2-")])
for (i in c("pfs"))#, "os"))
{ p = psSurv[[i]][[1]]; #p=p[rownames(p)!="TLS ST",]; # PFS, normal fig
  plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="", yaxt='n', xaxt='n',
    main=m<-paste("Other cancer types - immunotherapy -", toupper(i)),
    col=1+(rownames(p)=="TLS ST"));
  a = p.adjust(p[,1], 'fdr'); pl = max(p[a<.05,1])
  abline(v=0, col='grey'); abline(h=-log10(pl), col='grey');
  logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE, granularity=0)
  r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
  addAnnotArrows(side=1, xlab="HR", annotDir=c("Better", "Worse"), le=.05, cex=.8, posXlab=0, plotXlab=TRUE)
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
plot(p[,2], -log10(p[,1]), ylab="P-value", xlab="", yaxt='n', xaxt='n',
  main="Other cancer types - immunotherapy - RECIST",
  col=1+(rownames(pResp)=="TLS ST"));
nm = rownames(p); nm[p[,1]>pl]=NA; #nm[nm=="TLS ST"] = "TLS ST";
logAxis(ceiling(-log10(min(p[,1], na.rm=TRUE))), 0, side=2, inverted=TRUE)
r = trunc(range(exp(p[,2])*10))/10; r = seq(r[1], r[2], by=.1); axis(side=1, at=log(r), labels=r)
addAnnotArrows(side=1, xlab="OR", annotDir=c("Worse", "Better"), le=.05, cex=.8, posXlab=0, plotXlab=TRUE)
abline(v=0, col='grey'); abline(h=-log10(pl), col='grey');
o = order(p[,1]);
myAddTextLabels(p[o,2], -log10(p[o,1]), labels=nm[o], col.label=1+grepl("TLS ST", nm[o]), cex.label=.9,
    heightPad=0.3, widthPad=0.02);
text(x=par('usr')[1]+strwidth("M")*.5, y=par('usr')[4]-(1:3)*strheight("M")*1.5, pos=4, c(parse(text="bold('TLS ST')"),
    formatNiceP(p["TLS ST",1]), formatNiceP(exp(p["TLS ST",2]), pName="'OR'", italic=FALSE)))
dev.off();

x = data.frame(ID=unlist(lapply(ds, function(i) colnames(i$dn))), Cohort=rep(names(ds),
  sapply(ds, function(i) ncol(i$dn))), `DRFS - time (yrs)`=sus$DRFS[,1], `DRFS - event`=sus$DRFS[,2],
  `TLS sig`=unlist(lapply(ds, function(i) i$cs[,"TLS ST"])), `TLS sig (quartile by study)`=qrt, check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "tlsAllCliSmallA.xlsx"), row.names=FALSE)

w = which(cliIspy$arm=="Pembro")
x = data.frame(`TLS sig`=csI[w,"TLS ST"], pCR=cliIspy$pcr[w], type=cliIspy$Type[w], check.names=FALSE)
write.xlsx2(x, paste0(figDataDir, "tlsAllCliSmallB.xlsx"), row.names=FALSE)

st = paste(rep(names(dsIm), sapply(dsIm, function(i) nrow(i$g))), ty)
for (i in c("pfs", "os"))
{ s = suI[[i]]
  x = data.frame(ID=unlist(lapply(dsIm, function(i) rownames(i$cl))), study=st,
    s[,1], s[,2],  dIm, check.names=FALSE)
  colnames(x)[3:4] = paste(i, "-", c("time (yrs)", "event"))
  write.xlsx2(x, paste0(figDataDir, "tlsAllCliSmall", c(pfs="C", os="D")[i], ".xlsx"), row.names=FALSE)
}

d = dIm;
ps = lapply(colnames(suI), function(what)
{ s = suI[,what];
  p = t(apply(d, 2, function(i) getCoxP(coxph(s~i+ strata(st)), HR=TRUE))) #, subset=cancer=="Bladder")))
  p2 = t(apply(d, 2, function(i) coef(summary(coxph(s~d[,"TLS ST"]+i+ strata(st))))[, c('Pr(>|z|)', "coef")]))
  return(list(p, p2));
}); names(ps) = colnames(suI);
psSurv=ps;

set.seed(123);
cairo_pdf(paste0(figDir, "tlsAllCliRest.pdf"), height=3, width=15)
par(mfcol=c(1,5), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), cex.main=1);
p = pResp;
par(mar=c(2,3,2,.5))
#myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$wI & cliIspy$arm=="Ctr", ylab="TLS ST sig", subSide=3, subLine=-1.5,
#  main="I-SPY2 - TNBC cohort - control arm", addN=TRUE, alpha=.5)
myBoxplot(csI[,"TLS ST"], cliIspy$pcr, subset=cliIspy$arm=="Ctr", ylab="TLS ST sig", subSide=3, subLine=-1.5,
  main="I-SPY2 - control arm", addN=TRUE, alpha=.5, colPoints=colIHC[cliIspy$Type], ylim=c(min(csI[,"TLS ST"]), 1.35))
legend('topleft', legend=names(colIHC), fill=colIHC)
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

x = data.frame(ID=rownames(dIm), Study=st, as.matrix(suI$os), dIm, check.names=FALSE)
colnames(x)[3:4] = paste("OS", colnames(x)[3:4]);
write.xlsx2(x, file=paste0(figDataDir, "TLSimmunoOS.xlsx"), row.names=FALSE);

# Signatures from Charoentong on TLS
#####################################
w = idAn[,"Lymphoid nodule"]; w[!cli$hasTLS]=NA;
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
x = readRDS("~/tmp/regsFull.RDS");
cl= x[,"ann.Lymphoid.nodule"]>.25;

pdf(paste0(figDir, "aucTLS2.pdf"), width=7.5, height=2.5);
par(mfrow=c(1,3), mar=c(3,3,1.5,.5), mgp=c(1.5,.5,0))
for (i in c("loos.Lymphoid.nodule", "pred.Lymphoid.nodule", "sig.TLS.ST"))
{ z = roc(response=cl, predictor=x[,i], direction="<", ci=TRUE);
  plot(z, mar=par('mar'), mgp=par('mgp'))
  mtext(paste("AUC:", round(z$auc, digits=2)), side=1, line=-3, adj=0, at=.6, cex=.8)
  mtext(paste0("(", round(z$ci[1], digits=2), "-",round(z$ci[3], digits=2), ")"),
      side=1, line=-1.4, adj=0, at=.6, cex=.8)
  title(main=c(loos.Lymphoid.nodule="Regression (LOO)", pred.Lymphoid.nodule="Regression (Final)",
    sig.TLS.ST='TLS ST signature')[i])
}
dev.off();

z = data.frame(cl, x[,c("loos.Lymphoid.nodule", "pred.Lymphoid.nodule", "sig.TLS.ST")]);
colnames(z) = c("Presence TLS", "Regression (LOO)", "Regression (Final)", 'TLS ST signature')
write.xlsx2(z, file=paste0(figDataDir, "aucTLS2.xlsx"), row.names=FALSE);

# AUC TLS
###################
y = list(cli$annotClean[,"Lymphoid nodule"]>0, cli$hasTLS)
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

x = data.frame(ID=rownames(cli), `TLS by annotation`=y[[1]], `TLS by regression`=y[[2]],
  csPB[,c("TLS ST", "TLS Lundeberg")], csBulk[,c("TLS ST", "TLS Lundeberg")], check.names=FALSE)
colnames(x)[4:7] = paste(colnames(x)[4:7], c("PB", "PB", "Bulk", "Bulk"));
write.xlsx2(x, file=paste0(figDataDir, "aucTLS.xlsx"), row.names=FALSE);


# Genes in other TLS sigs
##########################
si = c(Lundeberg=list(sigs[["TLS Lundeberg"]]$name), Cabrita=unname(sigH["TLS Cabrita"]),
  Meylan=list(sigs[["TLS Meylan"]]$name))
pdf(paste0(figDir, "compTLSsigs.pdf"), width=15, height=5)
gB = names(idG); #g = names(ok)[ok]
m1 = rowMins(ms[,colnames(ms)!="Lymphocyte"])/log(2)
m2 = ms[,"Lymphocyte"]/log(2); names(m1) = names(m2);
par(mfrow=c(1,3), mar=c(3,3,3,3), mgp=c(1.5,.5,0));
for (i in names(si))
{ z = smoothScatter(log2(1+2^m2), log2(1+2^m1), xlab="Fold-change TLS vs. lymphocyte compartment",
    ylab='Min fold-change TLS vs. another compartment',
    pch=NA, nrpoint=100, ret.selection = TRUE, xaxt='n', yaxt='n', transformation = function(x) {pmin(x,2) ^.25}) # 
  a=c(0, .5, 1, 2, 4, 8, 16, 32); axis(1, at=log2(1+a), labels=a);
  a=c(0, .5, 1, 2, 4, 8); axis(2, at=log2(1+a), labels=a);
  abline(h=1, v=1, col='grey');

  z = intersect(si[[i]], rownames(ms));
  co = ifelse(z %in% gB, "black", "red");

  points(log2(1+2^m2[z]), log2(1+2^m1[z]), col=co)
  #abline(h=-10, v=-10, col='red'); 

  #m1b = c(m1[z], rep(seq(from=-2, to=.9, len=10), 10), rep(seq(from=.5, to=3, len=10), 10))
  #m2b = c(m2[z], rep(seq(from=-2, to=.9, len=10), each=10), rep(seq(from=-4, to=-1, len=10), each=10))
  m1b = c(m1[z], rep(seq(from=-2, to=.9, len=10), 10), rep(seq(from=-3, to=1.5, len=10), 10))
  m2b = c(m2[z], rep(seq(from=-2, to=.9, len=10), each=10), rep(seq(from=-2, to=5, len=10), each=10))
  myAddTextLabels(log2(1+2^m2b), log2(1+2^m1b), c(z, rep(NA, 200)), cex.label=.8, col.label=c(co, rep(NA, 200)),
    heightPad=.3, widthPad=.1, figBorder=.03)

  di = 1.6;
  arrows(c(-.1, .1)+log2(2), par('usr')[4]+di*strheight("I"), par('usr')[1:2], par('usr')[4]+di*strheight("I"), xpd=NA, length=.17)
  mtext(c("Lymphocytes", "TLS"), line=di/2+.7, at=(par('usr')[1:2]+1)/2, side=3)
  arrows(a<-(par('usr')[2]+di*strwidth("N")), c(-.1, .1)+log2(2), a, par('usr')[3:4], xpd=NA, length=.17)
  text(a+strwidth("N")*1.6, (par('usr')[3:4]+1)/2, xpd=NA, srt=270, c("Other non-lymphocytes","TLS")); 
  mtext(i, side=3, line=-1.2, font=2)
}
dev.off();

x = data.frame(Gene=rownames(ms), 2^cbind(m1, m2));
for (i in c("TLS ST", names(si)))
{ x[paste("Present in signature",i)] = x$Gene %in% c(`TLS ST`=list(gB), si)[[i]];
}
colnames(x)[2:3] = c('Min fold-change TLS vs. another compartment', "Fold-change TLS vs. lymphocyte compartment")
write.xlsx2(x, file=paste0(figDataDir, "compTLSsigs.xlsx"), row.names=FALSE);

# TLS in TNBC (controlled)
########################
ws = c("Immune1", "Immune2", "Inflammatory response", "Interferon alpha response",
  "Interferon gamma response", "IL2 STAT5 signaling", "IL6 JAK STAT3 signaling",
  "TLS Lundeberg", "TLS Cabrita", "TLS Meylan")
  
x = do.call(rbind, lapply(ds, function(i) scale(i$cs[,c("TLS ST", ws)])))
colnames(x)[colnames(x)=="Interferon alpha response"] = "IFN-α response";
colnames(x)[colnames(x)=="Interferon gamma response"] = "IFN-γ response";
x = data.frame(x, `IM subtype`=unlist(lapply(ds, function(i) i$cli$bar=="IM")), check.names=FALSE)
y = cbind(x, ctrTot, stud=stud)

ctr = list(METABRIC=ctrlMB, `Spatial TNBC`=ctrl[rownames(ds[["Spatial TNBC"]]$cli),], `SCAN-B`=ctrlSB)
ctrTot = do.call(rbind, ctr[names(ds)])
ct = update(ctrlForm, ~.+strata(stud))

f = update(ct, sus[[i]] ~ . + loc);
f2 = update(ct, sus[[i]] ~ . + loc + `TLS ST`);

nf=list();
for (iter in 1:2) {
cairo_pdf(paste0(figDir, "controlTLS_TNBC.pdf"), width=10, height=3);
par(mfrow=c(1,3), mar=c(.5,1,1.5,1), mgp=c(1.5,.5,0))
for (i in names(sus))
{ co = do.call(rbind, lapply(colnames(x)[-1], function(j)
  { y$loc = x[,j];
    fm = coxph(f, data=y)
    fm2 = coxph(f2, data=y)
    p = anova(fm, fm2)[2,"Pr(>|Chi|)"]
    c(p=p, coef(summary(fm2))["`TLS ST`", c("coef", "se(coef)")]);
  }))
  
  p1 = coef(summary(coxph(update(f2, ~ . -loc), data=y)))["`TLS ST`",c(5, 1, 3)]
  
  r = list();
  r[[1]] = c(list(expression(bold('TLS ST')), expression(bold('TLS ST controlled for'))), 
    colnames(x)[-1])
  r[[2]] = c(list(formatNice(p1[1]), NA), lapply(co[,'p'], formatNice))
  
  col = ifelse(c(p1[1], NA, co[,'p'])<.05, "#00AFBB", "slategray4") 

  m = rbind(p1[2:3], NA, co[,2:3])
  m = exp(m[,1] +1.96*cbind(-m[,2], 0, +m[,2]));
  
  nfo = nf[[i]]$figInfo;
  nf[[i]] = basicForest(r, m, lineHeight=cumsum(c(1.5,2,2,rep(1.5, nrow(m)-2))),
      xlab="HR", xlog=TRUE, xlim=c(.5, 1.1), annotDir=c("Better", "Worse"),
      titles=c("", expression(bolditalic(P))), col=col, beforePlot=bf)
  title(main=i);
}
dev.off();
}

cairo_pdf(paste0(figDir, "controlTLSby_TNBC.pdf"), width=10, height=2.6);
par(mfrow=c(1,3), mar=c(.5,1,1.5,1), mgp=c(1.5,.5,0))
for (i in names(sus))
{ b = allForest(y, sus[[i]], control=update(ct, ~.+TLS.ST), annotDir=c("Better", "Worse"), 
    clip=c(.5, 1.3), columns="P", lineHeight=c(1.5, (2:ncol(x))*1.5+1.5), colWidth=c(1,.3,.6))
  b = b$basicForest$linePos;
  text(1, (b[1]+b[2])/2, "Signatures controlled for TLS ST", adj=0, font=2)#, cex=par('cex'))
  title(main=i);
}
dev.off();

z = data.frame(ID=rownames(y), Study=rep(names(ds), sapply(ds, function(i) nrow(i$cli))),
  do.call(cbind, sus), y, check.names=FALSE);
colnames(z)[2+(1:6)] = paste(rep(names(sus), each=2), colnames(z)[2+(1:6)])
write.xlsx2(z, file=paste0(figDataDir, "controlTLS_TNBC.xlsx"), row.names=FALSE);

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

x = data.frame(ID=rownames(cliIspy), `Receptor`=cliIspy$Receptor.Subtype, Arm=cliIspy$Arm..short.name.,
  pCR=cliIspy$pcr, `TLS ST`=csI[,"TLS ST"], check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "tlsIspy.xlsx"), row.names=FALSE);

w = cliIspy$arm=="Pembro";
x = scale(csI[w,c("TLS ST", ws)])
colnames(x)[colnames(x)=="Interferon alpha response"] = "IFN-α response";
colnames(x)[colnames(x)=="Interferon gamma response"] = "IFN-γ response";

pcr = cliIspy$pcr[w]
p1 = wilcox.test(x[,1] ~ pcr)$p.value;
m1 = coef(summary(glm(pcr ~ x[,1], family=binomial)))[2,1:2]
p2 = t(apply(x[,-1], 2, function(i) coef(summary(glm(pcr ~ i + x[,1], family=binomial)))["x[, 1]", c(4,1,2)]))

r = list();
r[[1]] = c(list(expression(bold('TLS ST')), expression(bold('TLS ST controlled for'))), 
  colnames(x)[-1])
r[[2]] = c(list(formatNice(p1), NA), lapply(p2[,1], formatNice))

m = rbind(m1, NA, p2[,2:3])
m = exp(m[,1] +1.96*cbind(-m[,2], 0, +m[,2]));

bf = function()
{ if (is.null(nfo)) { return(); }
  rect(m[1,1], nfo$li[1], m[1,3], nfo$li[length(nfo$li)]+.5/length(nfo$li), col=alpha('lightgrey', .5), border=NA, xpd=NA)
}
nfo=NULL;

for (iter in 1:2)
{ cairo_pdf(paste0(figDir, "controlTLSispy.pdf"), height=3.5, width=5);
  par(mar=c(.5, .5, .5, 1));
  nfo = basicForest(r, m, lineHeight=cumsum(c(2,1.5,2,rep(1.5, nrow(m)-2))),
      xlab="OR", xlog=TRUE, xlim=c(.9, 5), annotDir=c("Worse", "Better"), titles=c("", expression(bolditalic(p))),
      beforePlot=bf, col = ifelse(c(p1[1], NA, p2[,1])<.05, "#00AFBB", "slategray4"))
  nfo = nfo$figInfo
  dev.off();
}

cairo_pdf(paste0(figDir, "controlTLSbyispy.pdf"), height=3.5, width=5);
par(mar=c(.5, .5, .5, 1));
bf = allForest(x, cliIspy$pcr[w], control=~TLS.ST, columns=c("P"), annotDir=c("Worse", "Better"),
  lineHeight=c(1.5, (2:ncol(x))*1.5+1.5), colWidth=c(1,.3,.6))
bf = bf$basicForest$linePos;
text(1, (bf[1]+bf[2])/2, "Signatures controlled for TLS ST", adj=0, font=2, cex=par('cex'))
dev.off();

z = data.frame(ID=rownames(cliIspy)[w], pCR=cliIspy$pcr[w], x, check.names=FALSE);
z = z[!is.na(z[,1]),]
write.xlsx2(z, file=paste0(figDataDir, "controlTLSispy.xlsx"), row.names=FALSE);

# TLS vs. GO
#############
# Note: based on selected GOs...
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

cutText = function(x, cex=1, maxSz=par('usr')[2]-par('usr')[1])
{ maxSz = abs(maxSz);
  s = strwidth(x, cex=cex)
  w = which(s>maxSz)
  for (i in w)
  { z = x[i];
    while(strwidth(z, cex=cex) > maxSz-strwidth("\U2026", cex=cex))
    { z = sub(" [^ ]+$", "", z); } 
    x[i] = paste0(z, "\U2026");
  }
  return(x);
}

gos = lapply(2:5, function(i)
{ read.xlsx(paste0(dataDir, "misc/GO_TLS_selection.xlsx"), sheetIndex=i) })
gos = lapply(1:4, function(i) { x = gos[[i]]; data.frame(x[,1], x[,2], i%%2) } )
goss = list(rbind( gos[[1]], gos[[2]]), rbind(gos[[3]], gos[[4]]) )
goss = lapply(goss, function(x) x[order(-x[,2]),])

a = psGo$GO;
a[,"Max other lower"] = rowMins(a[,grepl("Lower", colnames(a)) & colnames(a) != "Lymphocyte Lower"])
cols=c('lightblue', 'yellow')
p = cbind(a[,"Lymphocyte Higher"], a[,"Lymphocyte Lower"],
  a[,"Max other higher"], a[,"Max other lower"]);
  
for (j in 1:2)
{ cairo_pdf(paste0(figDir, "tlsGO", c(".VSlympho", ".VSrest")[j], ".pdf"),
    height=12, width=5)
  par(mar=c(3,.5,3,1), mgp=c(1.5,.5,0))
  #par(mar=c(3,20,3,1), mgp=c(1.5,.5,0))
  nm = goss[[j]][,1];
  ii = p[,j*2 + c(-1, 0)];
  i = cbind(2*rowMins(ii), ii[,1]<ii[,2]);
  i[,1] = p.adjust(i[,1], method='fdr'); i = i[nm,];
  o = order(-i[,1]); i=i[o,]; nm = nm[o];
  at = barplot(-log10(abs(i[,1])), horiz=TRUE, yaxt='n', xlab="FDR", #expression(paste(italic('P'), '-value')),
    col=cols[2-i[,2]], border=NA,
    yaxs='i', ylim=c(0, 1.2*nrow(i)), xaxt='n');
  a = axTicks(1); lbl = paste0("10^-", a); lbl[a==0]=1; axis(1, at=a, label=parse(text=lbl));
  text(0, at, cutText(paste0(" ",f(nm)), cex=.9), adj=0, cex=.9)

  title(c("TLS vs. lymphocyte compartment", "TLS vs. other non-lymphocytes")[j], line=2);
  mtext(list(c("TLS", "Lymphocyte compartment"), c("TLS", "Other non-lymphocytes"))[[j]], side=3,
      col=colorspace::darken(cols,.4), line=.5, adj=c(.25,.75), cex=par('cex')) 
  dev.off();
  
  x = data.frame(rownames(ii), ii, rownames(ii) %in% nm);
  colnames(x) = c("GO", paste(c("FDR TLS >", "FDR TLS <"), c("Lymphocyte", "Rest")[j]), "Displayed?")
  write.xlsx2(x, file=paste0(figDataDir, "tlsGO", c(".VSlympho", ".VSrest")[j], ".xlsx"), row.names=FALSE); 
  
  cairo_pdf(paste0(figDir, "tlsGO", c(".VSlympho", ".VSrest")[j], ".blueOnly.pdf"),
    height=5, width=4+1)
  par(mar=c(3,15,2,1), mgp=c(1.5,.5,0), cex.axis=3/4)
  #par(mar=c(3,.5,3,1), mgp=c(1.5,.5,0))
  i = i[i[,2]==1,];
  at = barplot(-log10(abs(i[,1])), horiz=TRUE, yaxt='n', xlab="FDR", #expression(paste(italic('P'), '-value')),
    col=cols[1], border=NA,
    yaxs='i', ylim=c(0, 1.2*nrow(i)), xaxt='n');
  a = axTicks(1); if (length(a)>5) { a = a[seq(1, length(a), by=2)]; }
  lbl = paste0("10^-", a); lbl[a==0]=1; axis(1, at=a, label=parse(text=lbl));
  #text(0, at, cutText(paste0(" ",f(rownames(i))), cex=.75), adj=0, cex=.75)
  mtext(cutText(f(rownames(i)), maxSz=par('mai')[2]*diff(par('usr')[1:2])/par('pin')[1], cex=.75),
    at=at, side=2, adj=1, cex=.75, las=2)

  #title(c("Higher in TLS vs. lymphocyte compartment", "Higher in TLS vs. other non-lymphocytes")[j],
  #  cex.main=1);
  mtext(c("Higher in TLS vs. lymphocyte compartment", "Higher in TLS vs. other non-lymphocytes")[j],
    outer=TRUE, line=-1.5, font=2);
  dev.off();
}

x = data.frame(rownames(p), p[,1:2]); colnames(x) = c("GO", c("FDR TLS > Lymphocyte", "FDR TLS < Lymphocyte"))
write.xlsx(x[,c(1,2)], file=paste0(figDataDir,"tlsGO.VSlympho.blueOnly.xlsx"), row.names=FALSE); 

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

#pdf(paste0(figDir, "fpTLSimmuno.pdf"), width=4.5, height=10)
pdf(paste0(figDir, "fpTLSimmuno.pdf"), width=7, height=2.5)
par(mar=c(.5,1,2.5,1));
#layout(matrix(1:3, ncol=1), height=c(.8,1,.9))
layout(matrix(1:3, nrow=1))
par(cex=2/3);
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

su = do.call(cbind, suI); colnames(su) = paste(rep(toupper(names(suI)), each=2), colnames(su));
x = data.frame(ID=unlist(lapply(dsIm, function(i) rownames(i$cl)))[w], Study=st2, `Cancer origin`=ty2,
  RECIST=s, su[w,], `TLS ST`=dIm[w,"TLS ST"], check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "fpTLSimmuno.xlsx"), row.names=FALSE);

## TLS correlogram
##################
load(paste0(dataDir, "classification/LOOstart.RData"))

n = do.call(rbind, annots);
n = cbind(n, Stroma=NA);
n = n[,!(colnames(n) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))]
n[,"Tumor"] = n[,"Tumor"] + n[,"Tumor region"];
n[,"Stroma"] = n[,"Stroma cell"] + n[,"Low TIL stroma"] + n[,"High TIL stroma"]/2 + n[,"Acellular stroma"];
n[,"Lymphocyte"] = n[,"Lymphocyte"] + n[,"High TIL stroma"]/2
n = n[,!(colnames(n) %in% c("Low TIL stroma", "High TIL stroma", "Stroma cell", "Acellular stroma",
  "Tumor region"))];
  
pr2 = pr; pts2 = pts;
w = rowSums(n)>1000; n=n[w,]; pr2=pr2[w,]; pts2=pts2[w];

load(paste0(dataDir, 'classification/loosReg.RData'));

y = n[,"Lymphoid nodule"] > .25
plot(roc(y, loos[["Lymphoid nodule"]]), print.auc=TRUE)

# 

tlsReg = rowsum(0+(loos[["Lymphoid nodule"]]>.01), pts2)[rownames(cli),1];

ist = paste0(dataDir, "classification/classifAll.RDS")
ii = do.call(rbind, lapply(ist, colMeans2))

tlsSig = grep("TLS", colnames(csPB), value=TRUE)
x = cbind(Annot=cli$annotClean[,"Lymphoid nodule"], Reg=ii[,"Lymphoid nodule"], Sig=csPB[,"TLS ST"])

pdf(paste0(figDir, "TLScor.pdf"), width=3, height=2.5)
par(xpd=NA);
myCorrplot(x, tl.srt=0, mar=c(.5,1,.5,1), tl.col='black', cl.ratio=.34, cl.cex=.7, cl.align.text='l',
  cl.offset=0.2)
dev.off();

## TLS vs. immune sig survival
###############################
getCoxP2 = function(x, HR)
{ s = summary(x)
  p=s$logtest[3]
  return(c(p, coef(s)[c(1, 3)]))
  #if (HR) { return(c(p, coef(x))); } else { return(p); }
}

ws = c("Immune1", "Immune2", "Inflammatory response", "Interferon alpha response",
  "Interferon gamma response", "IL2 STAT5 signaling", "IL6 JAK STAT3 signaling",
  "TLS Lundeberg", "TLS Cabrita", "TLS Meylan")#, "Hou")#, "TLS ST") 
ps = lapply(colnames(suI), function(what)
{ s = suI[,what];
  p = t(apply(dIm[,c("TLS ST", ws)], 2, function(i) getCoxP2(coxph(s~i+ strata(st)), HR=TRUE))) #, subset=cancer=="Bladder")))
  p2 = t(apply(dIm[,ws], 2, function(i) t(coef(summary(coxph(s~dIm[,"TLS ST"]+i+ strata(st))))[, c('Pr(>|z|)', "coef", "se(coef)")])))
  colnames(p2) = paste(c('Pr(>|z|)', "coef", "se(coef)"), rep(c("TLS ST", "Other"), each=3))
  return(list(p, p2));
}); names(ps) = colnames(suI);
pR = t(apply(dIm[,c("TLS ST", ws)], 2, function(i) coef(summary(glmer(resp ~ i+(1|st), family=binomial)))["i", c("Pr(>|z|)", "Estimate", "Std. Error")]))
p2R = t(apply(dIm[,setdiff(ws, "TLS ST")], 2, function(i)
  t(coef(summary(glmer(resp ~ dIm[,"TLS ST"]+i+(1|st), family=binomial)))[2:3, c("Pr(>|z|)", "Estimate", "Std. Error")])))
colnames(p2R) = paste(c('Pr(>|z|)', "coef", "se(coef)"), rep(c("TLS ST", "Other"), each=3))
p2S = lapply(ps, function(i) i[[2]]);
pS = lapply(ps, function(i) i[[1]]);
pS = c(pS, RECIST=list(pR)); p2S = c(p2S, RECIST=list(p2R));
  
bf = function()
{ if (is.null(nfo)) { return(); }
  rect(m[1,1], nfo$li[1], m[1,3], nfo$li[length(nfo$li)]+.5/length(nfo$li), col=alpha('lightgrey', .5), border=NA, xpd=NA)
}

ws2 = ws;
ws2[ws=="Interferon gamma response"] = "IFN-γ response";
ws2[ws=="Interferon alpha response"] = "IFN-α response";

nfoL = list();
cairo_pdf(paste0(figDir, "controlTLSimmunoT.all.pdf"), height=4, width=5*3);
par(mar=c(.5, .5, 2, 1), mfrow=c(1,3), cex=1);
for (i in names(pS))
{ p = pS[[i]]; p2 = p2S[[i]];
  if (i!="RECIST") { xlim=c(.6, 1.1); hr="HR"; ad=c("Better", "Worse"); }
  else { xlim = c(.8, 2); hr="OR"; ad=c("Worse", "Better"); }
  r = list();
  r[[1]] = c(list(expression(bold('TLS ST')), expression(bold('TLS ST controlled for'))), 
    ws2)

  r[[2]] = c(list(formatNice(p["TLS ST", 1]), NA), lapply(p2[,1], formatNice))

  m = rbind(p["TLS ST", 2:3], NA, p2[,2:3])
  m = exp(m[,1] + 1.96*cbind(-m[,2], 0, +m[,2]));
  nfo = NULL;
  #for (iter in 1:2)
  #{ pdf(paste0(figDir, "controlTLSimmunoT.", i, ".pdf"), height=3.5, width=5);
    #par(mar=c(.5, .5, .5, 1));
    nfo = nfoL[[i]]$figInfo
    nfoL[[i]] = basicForest(r, m, lineHeight=cumsum(c(1.5,2,2,rep(1.5, nrow(m)-2))),
        xlab=hr, xlog=TRUE, xlim=xlim, annotDir=ad, titles=c("", expression(bolditalic(p))),
        beforePlot=bf, col = ifelse(c(p["TLS ST", 1], NA, p2[,1])<.05, "#00AFBB", "slategray4"), colWidth=c(.8,.2,.5))
    
    #dev.off();
 # }
  title(main=toupper(i));
}
dev.off();

cairo_pdf(paste0(figDir, "controlTLSby_immunoT.all.pdf"), height=4, width=5*3);
par(mar=c(.5, .5, 2, 1), mfrow=c(1,3), cex=1);
x = data.frame(dIm[,c("TLS ST", ws)], st=st, check.names=FALSE)
colnames(x)[colnames(x)=="Interferon gamma response"] = "IFN-γ response";
colnames(x)[colnames(x)=="Interferon alpha response"] = "IFN-α response";

for (i in names(pS))
{ if (i!="RECIST") { xlim=c(.7, 1.3); hr="HR"; ad=c("Better", "Worse"); y=suI[[i]]}
  else { xlim = c(.6, 1.5); hr="OR"; ad=c("Worse", "Better"); y=resp; }
  #pdf(paste0(figDir, "controlTLSby_immunoT.", i, ".pdf"), height=3.5, width=5);
  #par(mar=c(.5, .5, .5, 1));
  bf = allForest(x, y, control=~TLS.ST+st, columns="P", clip=xlim, annotDir=ad, colWidth=c(.8,.2,.5),
    lineHeight=c(1.5, (2:ncol(x))*1.5+1.5));
  bf = bf$basicForest$linePos;
  text(1, (bf[1]+bf[2])/2, "Signatures controlled for TLS ST", adj=0, font=2, cex=par('cex'))
  #dev.off();
  title(main=toupper(i));
}
dev.off();

x = data.frame(id=rownames(dIm), `PFS - time`=suI$pfs[,1], `PFS - event`=suI$pfs[,2],
  dIm[,c("TLS ST", ws)], check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "controlTLSby_immunoT.all.xlsx"), row.names=FALSE);

x = data.frame(id=rownames(dIm), `OS - time`=suI$os[,1], `OS - event`=suI$os[,2], RECIST=resp,
  dIm[,c("TLS ST", ws)], check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "controlTLSby_immunoT.all.SI.xlsx"), row.names=FALSE);

## Comparison of TLS signatures
#################################
pdf(paste0(figDir, "compTLS.pdf"), width=7, height=7.5)
f = factor(idCan[,1], levels=c('Lymphoid nodule','Lymphocyte','Fat tissue','in situ',
  'Lactiferous duct','Necrosis','Stroma','Tumor','Vessels'));
par(mfrow=c(4,1), mar=c(6,3,1.5,1.5), mgp=c(1.5,.5,0))
for (i in c("TLS ST", "TLS Lundeberg", "TLS Cabrita", "TLS Meylan"))
{ myBoxplot(csAn[,i], f, leg45=TRUE, ylab=i)
  for (j in 2:length(levels(f)))
  { w = f %in% c('Lymphoid nodule',levels(f)[j]);
    p = wilcox.test(csAn[w,i] ~ f[w])$p.value
    mtext(formatNiceP(p), side=3, at=j, cex=par('cex'))
  }
}
dev.off()

x = data.frame(idCan[,2:1], csAn[,c("TLS ST", "TLS Lundeberg", "TLS Cabrita", "TLS Meylan")], check.names=FALSE);
colnames(x)[2:1] = c("Compartment", "Sample ID")
write.xlsx2(x, file=paste0(figDataDir, "compTLS.xlsx"), row.names=FALSE);

## Datasets for other figs
#############################
x = data.frame(ID=rownames(cli), Subtype=cli$barPB, `TLS ST signature`=csPB[,"TLS ST"], check.names=FALSE)
write.xlsx(x, file=paste0(figDataDir, "TLSsigVsBar.xlsx"), row.names=FALSE);

x = data.frame(ID=rownames(cli), TIME=cli$TIME_classes.bypathologist, `TLS ST signature`=csPB[,"TLS ST"], check.names=FALSE)
write.xlsx(x, file=paste0(figDataDir, "TLSsigVsTIME.xlsx"), row.names=FALSE);

x = data.frame(`Sample ID`=idCan[,2], `Compartment`=idCan[,1], csAn[,grep("^TLS", colnames(csAn))], check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "TLSbyCompartment.xlsx"), row.names=FALSE);

#######################
## MC & ET
#######################

## Bareche by cluster
######################
t = table(idC$bar, idC[,1]); t = t[,rownames(cli)];
o = order(cli$barPB, -t[cbind(as.character(cli$barPB), colnames(t))]/colSums(t))

pdf(paste0(figDir, "barByClust.pdf"), height=2.6, width=15)
par(mar=c(4,3,.5,9), mgp=c(1.5,.5,0));
at = barplot(t[,o], col=colBar, xaxt="n", ylab="N clusters", xaxs='i', xlim=c(-1, ncol(t)*1.2+1))
mtext(rownames(cli)[o], line=0, cex=.6, side=1, at=at)
addAnnot("Subtype PB", colBar[cli$barPB[o]], line=1, side=1, at=at, textAtStart=FALSE)

b = cli$barPBco[o,];
b = b/rowMaxs(b); b2 = apply(b,1, function(i) { o = order(i); c(o[4], i[o[4]]); })
b2[2,b2[2,]<0] = 0;
#addAnnot("Subtype PB 2nd class", colBar[colnames(b)[b2[1,]]], line=2, side=1, at=at, heights=b2[2,], textAtStart=FALSE)

f = cli$TIME_classes.bypathologist; f[f=="nd"]=NA; f=factor(f, levels=c("ID", "MR", "SR", "FI"));
addAnnot("TIME", colTIME[f[o]], line=2, side=1, at=at, textAtStart=FALSE)
legend(-15, 7.5, names(colBar), fill=colBar, xpd=NA, bty='n', xjust=0, yjust=1)
text(-14, 7.7, 'Subtype', xpd=NA, adj=0);
legend(-15, 3, levels(f), fill=colTIME, xpd=NA, bty='n', xjust=0, yjust=1)
text(-14, 3.2, 'TIME', xpd=NA, adj=0);
dev.off();

x = data.frame(ID=idC[,1], `ID cluster`=as.integer(idC[,2]), `Subtype cluster`=idC$bar,
 `Subtype sample PB`=idC$barPB, `TIME sample`=cli[idC$id, "TIME_classes.bypathologist"], check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "barByClust.xlsx"), row.names=FALSE);

# MC descriptions
###################
anByC = do.call(rbind, tapply(seq_along(K), K, function(i) colMeans(annot2[i,], na.rm=TRUE)));
anByC = anByC[14:1,c("Tumor", "Stroma", "Lymphocyte", 
  "Fat tissue", "Necrosis", "in situ", "Lymphoid nodule", "Vessels", "Lactiferous duct")]
co = colAnn2[colnames(anByC)]; co[2] = colAnn2["Stroma cell"]; names(co)[2]="Stroma";

#descrs = read.xlsx(paste0(dataDir, "misc/Signature summary - tblsMC2.xlsx"), 1)$Characteristics;
#descrs = descrs[!is.na(descrs)];
#w = gregexpr(", ", descrs); lim=40;
#for (i in seq_along(w)) # Cut in lines of max lim characters
#{ ww = c(w[[i]], nchar(descrs[i]));
#  if (all(ww<40)) { next; }
#  cu = ww[which(ww>=40)-1];
#  descrs[i] = paste0(substr(descrs[i], 1, cu), "\n", substr(descrs[i], cu+2, nchar(descrs[i]))) 
#}

#sig2keep = unlist(read.xlsx(paste0(dataDir, "misc/megaclusters_characterization.xlsx"), 2))
#sig2keep = sub(" $", "", sig2keep[!is.na(sig2keep)]);
#sig2keep[sig2keep=="VCpred_TN"] = "VCpred TN";
#sig2keep[sig2keep=="TNFA signaling via NF−kB"] = "TNFA signaling via NF-kB"
#sig2keep = setdiff(sig2keep, c("Complement", "Immune2", "Allograft rejection", "CIN70", "Myogenesis"));

tmp = read.xlsx(paste0(dataDir, "misc/tblsMC2.annot.xlsx"), 1)
tmp[,1] = sub(" $", "", tmp[,1]); 
tmp[,1][tmp[,1]=="EMT"] = "Epithelial mesenchymal transition";
w = tmp[,1]=="Macrophages M2"; tmp[w,1:2] = "Macrophages";
rownames(tmp) = tmp[,1];
tmp = tmp[order(tmp$Category), ]
ct = intersect(tmp[,1], colnames(xCn.xc))
gns = intersect(tmp[,1], rownames(xCn));
sig2keep = setdiff(rownames(tmp), c(ct, gns));
nok = setdiff(sig2keep, colnames(xCn.cs));
ad = adist(nok, colnames(xCn.cs))
tmp[nok,1] = colnames(xCn.cs)[apply(ad, 1, which.min)];
rownames(tmp) = tmp[,1];
sig2keep = setdiff(rownames(tmp), c(ct, gns));

#sig2keep = setdiff(sig2keep, c("ESR1 gene", "Estrogen response late", "NOTCH signaling",
#  "WNT beta catenin signaling", "E2F targets", "Xenobiotic metabolism", "Heme metabolism",
#  "TNFA signaling via NF-kB"))

#gns = c("NT5E", "PALB2", "RAD51C", "TRDMT1", "MSH6", "MSH2", "CD47");
#ct = c("Th2 cells", "DC", "Monocytes", "Macrophages M2", "Neutrophils", "Adipocytes");

x = xCn.cs;
x = x[,sig2keep];
x = cbind(x, xCn.xc[,ct], t(xCn[gns,]))
x = x[, rownames(tmp)];
xb = x;

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
addT = 20; til = apply(pres,2,function(i) cli$sTILs.percentage.bypathologist[i]+addT)
#idC$bar = factor(TNBCclassif(y, version='bareche', shortName=TRUE), levels=names(colBar));

lbl = tmp[colnames(x),2];
#lbl[lbl=="Epithelial mesenchymal transition"]="EMT";
#lbl = sub("signaling", "sig.", lbl);
#lbl = sub("response", "resp.", lbl);

pdf(paste0(figDir, "tblsMC2.pdf"), width=8.5, height=4)
par(mar=c(5,3,.5,.3), oma=c(0,0,0,5), mgp=c(1.5,.5,0), xpd=NA);
#layout(mat=matrix(1:5, nrow=1), width=c(.5,.8,.8,3.3,1.6))
layout(mat=matrix(1:4, nrow=1), width=c(.5,.8,.8,3.3))
a = rev(colSums(pres)); names(a) = paste0("MC", names(a)); names(a) = sub("k", "", names(a));
at = barplot(a, horiz=TRUE, xlab="N samples", xaxt='n', las=2, col=rev(MC.colors), xpd=NA)
axis(1)
par(mar=c(5,.5,.5,.5))
barplot(t(anByC)*100, xlab="% Annotation", horiz=TRUE, yaxt='n', col=co)
usr = par('usr');
t = table(idC$bar, K); t = t/rep(colSums(t), each=5); t=t[names(colBar),];
barplot(t[,14:1]*100, col=colBar, xlab="% TNBC subtype", horiz=TRUE, yaxt='n') 

plot.new(); plot.window(xlim=c(0,1), ylim=usr[3:4], yaxs='i', mar=c(3,0,.5,.5), xaxs='i');
for (i in 1:14)
{ subplot(image(t(x2[Ko==i,]), axes=FALSE, zlim=c(-2, 2), col=col), x=c(0,1), rev(at)[i]+c(-.5,+.5)); 
}
text(seq(0,1,len=ncol(x2)+1)[-1]-.5/ncol(x2), 0, lbl, xpd=NA, srt=-40, adj=0, cex=.9,
  col=colSig[unclass(factor(tmp$Category))])# colSig[sigInfo[colnames(x2)]]) 

plotScale(cols=col, v=c("Low", "High"), atV=c(0,1), posx=par('usr')[2]+strwidth("N")*1,
  posy=5, horizontal=FALSE, width=2, height=15)

#plot.new(); plot.window(xlim=c(0,1), ylim=usr[3:4], yaxs='i', mar=c(3,0,.5,.5));
#par(lheight=.7)
#text(rep(-.05, 14), at, rev(descrs), adj=0, xpd=NA)

#subplot(image(matrix(1:25, ncol=1), axes=FALSE, col=col, xpd=NA, useRaster=TRUE), x=c(.2,.8), y=usr[3]-3+c(1,2))
#text(x=c(.2, .8), y=rep(usr[3]-.7, 2),c("Low", "High"), xpd=NA, pos=c(4,2));

#plot.new(); par(mar=c(0,0,.5,0)); plot.window(xlim=c(0,1), ylim=c(0,1));
#legend('topleft', legend=colnames(anByC), fill=co, title="Annotations", xpd=NA, inset=c(-.4, .1))

dev.off();

xx = data.frame(`ID sample`=idC[,1], `ID cluster`=idC[,2], Megacluster=paste0("MC", K),
  `Subtype cluster`=idC$bar,
  annot2, xb, check.names=FALSE);
write.xlsx2(xx, file=paste0(figDataDir, "tblsMC2.xlsx"), row.names=FALSE)

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
  { o = rev(intersect(names(colXct), colnames(cc))); cc = cc[,o];
    #o = order(-match(colnames(cc), names(colXct))); cc = cc[,o];
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
  cutOff = max(z[p1$fdr[w,]<.05])+1e-10;
  
  cairo_pdf(paste0(figDir, "MCcomp.", what, ".pdf"),
      height=(length(w)+7)/7+1, width=6+10+max(nchar(w))/20)
  par(mar=c(8,max(nchar(w)+4)/2.5,.5,.5), mgp=c(1.5,.5,0), cex=1.1);
  layout(matrix(1:2, nrow=1), width=c(.4,.6));
  
  dotPlot(cc[,w,drop=FALSE], K2, col.lbl=co[w], maxP=cutOff,legend=FALSE, cex.pch=.45, srt=45);
  legendDotPlot(-1,-6*strheight("M"), horizontal=TRUE, cex.pch=.45, interline=1.7, significantLabel=c("FDR < 5%", "FDR ≥ 5%"))

  
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
  
  x = data.frame(`Sample ID`=idC[,1], ClustNr=idC[,2], Megacluster=K2, cc, check.names=FALSE)
  write.xlsx2(x, file=paste0(figDataDir, "MCcomp.", what, ".xlsx"), row.names=FALSE);
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

x=matrix(ncol=length(ploo), nrow=max(sapply(ploo, length)))
for (i in seq_along(ploo)) { x[seq_along(ploo[[i]]),i] = ploo[[i]]; }
colnames(x) = paste(names(ploo), "MCs");
write.xlsx(x, file=paste0(figDataDir, "meanAUC.xlsx"), row.names=FALSE);


## Number of SA
###############
pdf(paste0(figDir, "Ncluster.pdf"), height=3.5, width=4);
par(mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
#nb = NbClust2(log10(100*(faN)), distance = "euclidean", min.nc = 5, max.nc = 15,
#    method="ward.D2", index ="allNonGraph"); 
nb = NbClust(log10(100*(faN)), distance = "euclidean", min.nc = 5, max.nc = 15,
    method="ward.D2"); 
barplot(table(nb$Best.nc[1,]), ylab="N indices", xlab='N Spatial Archetypes');

dev.off();


## Heatmap ET
################
hc = hclust(dist(log10(100*(faN))), method='ward.D2');
cu = which(diff(cutree(hc, max(ecot))[hc$order])!=0);

ff = faN; colnames(ff) = paste0("MC", colnames(ff));
pdf(paste0(figDir, "heatmapMC.pdf"), width=13, height=10)
par(cex.axis=1.3, mgp=c(1,.5,0), oma=c(0,1,0,0))
h = heatmap.3(log10(100*(t(ff))), scale='none', col=colorRampPalette(c('blue', 'yellow'), space='Lab')(100),
  Rowv=FALSE, key=FALSE,
  Colv=as.dendrogram(hc), colsep=cu,  dendrogram='column',
  zlim=c(0,2), ColSideColors=cbind(colBar[cli$barPB],
    #`TIME (expression)`=colTIME[c(`Margin restricted`='ID', `Fully Inflamed`='FI', `Stroma Restricted`='SR')[cli$TIMEpb]],
    colTIME[cli$TIME_classes.bypathologist],
    ET.colors[ecot], ET.colors[ecot]),
    lhei=c(2,10), cexRow=1.8, lwid=c(.75,10), labRow=NA,
  margins=c(2,.5), trace='none', sepcolor="white", side.height.fraction=1)#,
par(new=TRUE); plot.new();
layout(matrix(1:4, ncol=2), height=c(2,10), width=c(.75,10))
par(cex=0.2 + 1/log10(94)); par(mar=c(0,0,0,.5), oma=c(0,1,0,0));
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1), xaxs='i');
w = c(0, cumsum(tabulate(ecot))); w = (w[-1]+w[-length(w)])/(2*w[length(w)])
w = w*94.5/94;
mtext(paste0("SA", 1:max(ecot)), side=3, cex=1.3, col='black', at=w, font=2, line=-.95)
mtext("Spatial", side=3, cex=1.2, col='red', at=-.003, font=1, line=-.3, adj=1);
mtext("Archetype", side=3, cex=1.2, col='red', at=-.003, font=1, line=-1.3, adj=1);
mtext(c("TIME", "Subtype"), side=3, cex=1.2, at=-.003, line=c(-2.8,-4.1), adj=1);

mtext(paste0("MC", 1:14), side=3, cex=1.8, at=-.003, line=-7 - (0:13)*3.8, adj=1 )

#par(new=TRUE, mar=c(0,0,0,0), xpd=NA, cex=1)
#plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1), xaxs='i');
#legend(-.15, .9, legend=c("Relapse", "No relapse"), fill=c('red', 'black'), xjust=0, title="    Relapse", bty='n',
#  title.adj=0)
#legend(-.15, .9, legend=names(colTIME), fill=colTIME, xjust=0, title="    TIME", bty='n', title.adj=0)
#legend(-.15, .6, legend=names(colBar), fill=colBar, xjust=0, title='    Subtype', bty='n', title.adj=0)
dev.off();     

x = data.frame(ID=rownames(cli), Subtype=cli$barPB, TIME=cli$TIME_classes.bypathologist, ff)
write.xlsx(x, file=paste0(figDataDir, "heatmapMC.xlsx"), row.names=FALSE)

pdf(paste0(figDir, "heatmapMC.legend.pdf"), width=5, height=4)
par(mar=c(.5,.5,.5,.5))
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1), xaxs='i');
tw = strwidth("MC12");
a = legend(0, 1, legend=names(colTIME), fill=colTIME, xjust=0, title="   TIME", bty='n', title.adj=0, title.font=2,
  ncol=5, text.width=tw)
h = 1-a$rect$h-strheight("M")
a = legend(0, h, legend=rev(names(colBar)), fill=rev(colBar), xjust=0, title='   TNBC molecular subtypes',
  bty='n', title.adj=0, title.font=2, ncol=5, text.width=tw)
h = h-a$rect$h-strheight("M")
a = legend(0, h, legend=names(MC.colors), fill=MC.colors, xjust=0, title='   Megaclusters',
  bty='n', title.adj=0, title.font=2, ncol=5, text.width=tw)
h = h-a$rect$h-5*strheight("M")
text(0,h+3.5*strheight("M"), "   Color scale", adj=0, font=2)
plotScale(cols=colorRampPalette(c('blue', 'yellow'), space='Lab')(100), c("Zero", "High"), atV=0:1,
  posx=.05, posy=h, horizontal=TRUE, width=14, height=1.5)
dev.off();

## Number of ET
################
nb = NbClust2(log10(100*(faN)), distance = "euclidean", min.nc = 5, max.nc = 15,
    method="ward.D2", index ="allNonGraph"); 
pdf(paste0(figDir, "N_ET.pdf"), height=3, width=4);
par(mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
barplot(table(nb$Best.nc[1,]), ylab="N indices", xlab="N Spatial Archetypes");
dev.off();

x = nb$Best.nc[1,]; x=x[!is.na(x)];
x = data.frame(`Index name`=names(x), `Best number of clusters`=x, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "N_ET.xlsx"), row.names=FALSE);

# Content ET
#############
n = do.call(cbind, tapply(seq_along(ecot), ecot, function(i) colMeans(fa[i,])))
n = 100*n/rep(colSums(n), each=nrow(n));
colnames(n) = paste0("SA", colnames(n));

pdf(paste0(figDir, "barplotMCvsET.pdf"), height=2, width=3)
par(mar=c(2,3,1.5,.5), mgp=c(1.5,.5,0), cex=2/3);
at = barplot(n, col=MC.colors, ylab="% MC", xaxt='n')
#legend(x=par('usr')[2], y=par('usr')[4], paste0("MC", 1:nrow(n)), fill=MC.colors, xpd=NA)
mtext(colnames(n), side=1, at=at, line=par("mgp")[2], cex=par('cex'))#, col=ET.colors);
title(main="ST TNBC cohort", cex.main=1)
dev.off();

x = fa; colnames(x) = paste0("MC", colnames(x));
x = data.frame(ID=rownames(cli), `Spatial Archetype`=paste0("SA", ecot), x, check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "barplotMCvsET.xlsx"), row.names=FALSE);

pdf(paste0(figDir, "barplotETs.pdf"), height=2, width=3)
par(mar=c(2,3,1.5,.5), mgp=c(1.5,.5,0), cex=2/3)
e = paste0("SA", ecotRec);
b = unlist(lapply(ds, function(i) as.character(i$cli$bar)))
#nm = c("ST global pseudobulk", "TNBC cohorts")[i]
b = factor(b, levels=names(colBar));
t = table(b,e);
t=t/rep(colSums(t), each=nrow(t));
at = barplot(100*t, col=colBar, ylab="% samples", xaxt='n')
title(main="ST TNBC (bulk), METABRIC and SCAN-B cohorts", cex.main=1)
mtext(colnames(t), side=1, at=at, line=par("mgp")[2], cex=par('cex'))#, col=ET.colors);
dev.off();

x = data.frame(ID=unlist(lapply(ds, function(i) rownames(i$cli))), 
  study=rep(names(ds), sapply(ds, function(i) nrow(i$cli))), `Spatial Archetype`=e, subtype=b, check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "barplotETs.xlsx"), row.names=FALSE);

pdf(paste0(figDir, "barEcot.pdf"), width=12, height=2);
par(mfcol=c(1,4), mar=c(2,3,1.5,.5), mgp=c(1.5,.5,0))
for (i in 1:4)
{ e = paste0("SA", list(ecot, ecotRec[stud=="S"], ecotRec[stud=="M"], ecotRec[stud=="SC"])[[i]]);
  b = list(cli$barPB, ds$`Spatial TNBC`$cli$bar, ds$METABRIC$cli$bar, ds$`SCAN-B`$cli$bar)[[i]]
  nm = c("ST global pseudobulk", "ST bulk", "METABRIC", "SCAN-B")[i]
  b = factor(b, levels=names(colBar));
  t = table(b,e);
  h = barplot(t, col=colBar, main=nm, ylab="N samples", xaxt='n')
  mtext(colnames(t), side=1, at=h, line=par("mgp")[2], cex=par('cex'))#, col=ET.colors);
  t=t/rep(colSums(t), each=nrow(t));
  #barplot(100*t, col=colBar, main=nm, ylab="% samples", xaxt='n')
  #mtext(colnames(t), side=1, at=h, line=par("mgp")[2], cex=par('cex'))#, col=ET.colors);
}
dev.off();

x = data.frame(ID=c(rownames(cli), unlist(lapply(ds, function(i) rownames(i$cli)))),
  Study=rep(c("Global ST PB", names(ds)), c(nrow(cli), sapply(ds, function(i) nrow(i$cli)))),
  Subtype=c(as.character(cli$barPB), unlist(lapply(ds, function(i) as.character(i$cli$bar)))),
  `Spatial Archetype`=paste0("SA", c(ecot, ecotRec)))
write.xlsx2(x, file=paste0(figDataDir, "barEcot.xlsx"), row.names=FALSE);
  
pdf(paste0(figDir, "barEcot2.pdf"), width=6, height=4);
par(mfcol=c(2,2), mar=c(2,3,1.5,.5), mgp=c(1.5,.5,0), cex=2/3)
for (i in 1:2)
{ e = paste0("SA", list(ecot, ecotRec)[[i]]);
  b = list(cli$barPB, unlist(lapply(ds, function(i) i$cli$bar)))[[i]]
  nm = c("ST global pseudobulk", "TNBC cohorts")[i]
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
for (i in c("PB", "All2", "All3", names(ds)))
{ if (i %in% c("All2", "All3"))
  { sus = susDs;
    x = do.call(rbind, lapply(ds, function(i) i[["15"]]$H));
    e = ecotRec;
    id = unlist(lapply(ds, function(i) rownames(i$cli)))
  }
  if (i=="PB")
  { sus = cli[,sapply(cli, is, "Surv")]; sus$iDFS=NULL;
    x = fa;
    e = ecot;
    id = rownames(cli);
  }
  if (i %in% names(ds))
  { x = ds[[i]]$cli;
    id = rownames(x);
    sus = x[,sapply(x, is, "Surv")];
    if (ncol(sus)>4) { sus = sus[,2:5]; }  # Remove RFS for Bulk
    x = ds[[i]][["15"]]$H;
    e = ds[[i]]$ecot[,1];
  }
  colnames(x) = paste0("MC", 1:ncol(x));
  xb = x;
    
  x[x<.01] = .01;
  x = scale(log10(x)); 
  
  x2 = matrix(FALSE, nrow=nrow(x), ncol=max(ecot));
  x2[cbind(1:nrow(x2), e)]=TRUE;
  colnames(x2) = paste0("SA", 1:ncol(x2));
  
  ctr = list(All3=ctrTot, All2=ctrTot, PB=ctrl[rownames(cli),], METABRIC=ctrlMB, `Spatial TNBC`=ctrl[rownames(ds[[i]]$cli),], `SCAN-B`=ctrlSB)[[i]];

  nm = c(All3="All3", All2="All2", PB="PB", METABRIC="MB", `Spatial TNBC`="Bulk", `SCAN-B`="SCANB")[i]
  nm2 = c(All3="All 3 studies", All2="METABRIC and SCAN-B", B="Pseudo-bulk", METABRIC="METABRIC",
  	`Spatial TNBC`="ST cohort (bulk RNA-seq)", `SCAN-B`="SCAN-B")[i]
  
  if (i%in%c("All2", "All3")) { x = data.frame(x, stud=stud); x2 = data.frame(x2, stud=stud); ct1 = ~strata(stud); ct2=ct; }
  else { ct1 = ~1; ct2 = ctrlForm; }
  if (nm=="Bulk") { clip = c(.2, 5) } else { clip = c(.5, 2); }
  
  sbs = rep(TRUE, length(e)); if (i=="All2") { sbs = stud!="S"; }
  cairo_pdf(paste0(figDir, "forestMC.", nm, ".pdf"), width=8.5, height=14.5)
  par(mar=c(.1,1,1,1));
  layout(mat=matrix(1:8, ncol=2, byrow=TRUE), heights=rep(c(1.05, 1), 4))
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(x, sus[[j]], control=ct1, clip=clip, annotDir=c("Better", "Worse"),
      lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE, subset=sbs)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Univariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  if (length(sus)==3) { plot.new(); }
  
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(data.frame(x, ctr), sus[[j]], control=ct2,
      clip=clip, annotDir=c("Better", "Worse"), lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE, subset=sbs)
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
      lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE, subset=sbs)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Univariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  if (length(sus)==3) { plot.new(); }
  
  for (j in seq_along(sus))
  { if (j<=2) { par(mar=c(.1,1,4,1));} else { par(mar=c(.1,1,1,1));}
    allForest(data.frame(x2, ctr), sus[[j]], control=ct2,
      clip=clip, annotDir=c("Better", "Worse"), lineHeight=2, colWidth=colWidth, sigOnFDR=TRUE, subset=sbs)
    mtext(names(sus)[[j]], side=3, line=0, cex=par('cex'), font=2);
    if (j==1)
    { mtext(paste0("Multivariable - ", nm2), side=3, line=2, cex=par('cex')*1.3, font=2, at=11);
    }
  }
  dev.off();
  
  z = data.frame(ID=id, ctr, check.names=FALSE);
  if (i%in%c("All2", "All3")) { z = cbind(z, Study=rep(names(ds), sapply(ds, function(i) nrow(i$cli)))); x$stud=NULL; }
  ss = do.call(cbind, sus); colnames(ss) = paste(rep(names(sus), 2), colnames(ss))
  z1 = cbind(z, ss, xb)[sbs,];
  write.xlsx2(z1, file=paste0(figDataDir, "forestMC.", nm, ".xlsx"), row.names=FALSE);
  z1 = cbind(z, ss, `Spatial Archetype`=paste0("SA", e)); z1 = z1[sbs,];
  write.xlsx2(z1, file=paste0(figDataDir, "forestET.", nm, ".xlsx"), row.names=FALSE);
}

# Tailored iBCFS
colWidth = c(1.5,1.4,3.2,2.5,2,8);
cairo_pdf(paste0(figDir, "forestMCtaylored.pdf"), width=10, height=3)
par(mfrow=c(1,2), mar=c(1,1,2,1), mgp=c(1.5,.5,0), cex=.6)
ctr= ctrl[rownames(cli),];
x = fa;
colnames(x) = paste0("MC", 1:ncol(x));
x[x<.01] = .01;
x = scale(log10(x))  
allForest(data.frame(x, ctr), cli$iBCFS, control=ctrlForm, colWidth=colWidth, sigOnFDR=TRUE,
	lineHeight=1.7, annotDir=c("Better", "Worse"), cex.axis=.9, cex.annot=.9)
title(main="iBCFS – multivariate analysis – ST TNBC cohort")

x = do.call(rbind, lapply(ds, function(i) i[["15"]]$H));
ctr = ctrTot;
ct2=update(ctrlForm, ~.+strata(study))
colnames(x) = paste0("MC", 1:ncol(x));
x[x<.01] = .01;
x = scale(log10(x))  
x = data.frame(x, stud=stud); ct2=update(ctrlForm, ~.+strata(stud)); 
allForest(data.frame(x, ctr), susDs$iBCFS, control=ct2, colWidth=colWidth, sigOnFDR=TRUE,
	subset=stud != "S", clip=c(.45, 2.1), lineHeight=1.7, annotDir=c("Better", "Worse"), cex.axis=.9, cex.annot=.9)
title(main="iBCFS – multivariate analysis – METABRIC and SCAN-B cohorts")
dev.off();

x = fa; colnames(x) = paste0("MC", 1:14);
x = data.frame(ID=rownames(cli), `iBCFS - time`=cli$iBCFS[,1], `iBCFS - event`=cli$iBCFS[,2],
  x, ctrl[rownames(cli),], check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "forestMC.PB.xlsx"), row.names=FALSE);

x = do.call(rbind, lapply(ds, function(i) i[["15"]]$H));
colnames(x) = paste0("MC", 1:14);
x = data.frame(ID=unlist(lapply(ds, function(i) rownames(i$cli))),
  Study=rep(names(ds), sapply(ds, function(i) nrow(i$cli))), x,
  `iBCFS - time`=susDs[["iBCFS"]][,1], `iBCFS - event`=susDs[["iBCFS"]][,2], ctrTot, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "forestMC.All3.xlsx"), row.names=FALSE);

x = data.frame(ID=unlist(lapply(ds, function(i) rownames(i$cli))),
  Study=rep(names(ds), sapply(ds, function(i) nrow(i$cli))), `Spatial Archetype`=paste0("SA", ecotRec),
  `iBCFS - time`=susDs[["iBCFS"]][,1], `iBCFS - event`=susDs[["iBCFS"]][,2], ctrTot, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "forestET.All3.xlsx"), row.names=FALSE);

## KM Survival Ecotypes
###############################
bar = unlist(lapply(ds, function(i) as.character(i$cli$bar)));
a = ifelse(ecotRec!=4 & bar=="IM", "IM not SA4", "Others") 
et4im = ifelse(ecotRec==4 & bar=="IM", "IM SA4", "Others")
et4 = ifelse(ecotRec==4, "SA4", "Others")
b = a; b[et4im!="Others"]=et4im[et4im!="Others"];
pdf(paste0(figDir, "survEcotypeNew2.pdf"), height=7.5, width=7.5);
par(mfrow=c(3,3), mar=c(2.5,4.5,1.5,.5), mgp=c(1.5,.5,0), oma=c(0,1,0,0))
for (i in names(susDs))
{ su = susDs[[i]];
  plotSurv(su, et4, ylab=toupper(i), xlab="Time (yrs)",
    legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c( "black", ET.colors["SA4"]),
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
    legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["SA4"], palette()[2]),
    pval = getCoxP(coxph(su ~ b+strata(stud), subset=b!="Others")))
  plotSurv(su, ifelse(ecotRec==8, "SA8", "Others"), ylab=toupper(i), xlab="Time (yrs)",
    legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c("black", ET.colors["SA8"]),
    pval = getCoxP(coxph(su ~ (ecotRec==8)+strata(stud))))
}
dev.off();

x = data.frame(ID=unlist(lapply(ds, function(i) rownames(i$cli))),
  Study=rep(names(ds), sapply(ds, function(i) nrow(i$cli))), Bareche=bar,
  `Spatial Archetype`=paste0("SA", ecotRec), `Spatial Archetype comp 1`=et4, `Spatial Archetype comp 2`=b,
  `Spatial Archetype comp 3`=ifelse(ecotRec==8, "SP8", "Others"),
  `iBCFS - time`=susDs[["iBCFS"]][,1], `iBCFS - event`=susDs[["iBCFS"]][,2], check.names=FALSE)
write.xlsx2(x, file=paste0(figDataDir, "survEcotypeNew2.xlsx"), row.names=FALSE);

for (i in c("ST", names(ds)))
{ if (i!="ST")
  { z = ds[[i]]; nm = c(METABRIC="MB", `Spatial TNBC`="Bulk", `SCAN-B`="SCANB")[i];
    e = z$ecot[,1]; cl = z$cli;
  } else { e=ecot; cl=cli; nm="ST"; }
  a = ifelse(e!=4 & cl$bar=="IM", "IM not SA4", "Others") 
  et4 = ifelse(e==4, "SA4", "Others")
  et4im = ifelse(e==4 & cl$bar=="IM", "IM SA4", "Others")
  b = a; b[et4im!="Others"]=et4im[et4im!="Others"];
  w = names(cl)[sapply(cl, is, "Surv")];
  pdf(paste0(figDir, "survEcotype.", nm, ".pdf"), height=2.5*length(w), width=7.5);
  par(mfrow=c(length(w),3), mar=c(2.5,4.5,1.5,.5), mgp=c(1.5,.5,0), oma=c(0,1,0,0))
  for (j in w)
  { su = cl[[j]];
    plotSurv(su, et4, ylab=toupper(j), xlab="Time (yrs)",
      legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c("black", ET.colors["SA4"]))
    #plotSurv(su, a, ylab=toupper(j), xlab="Time (yrs)",
    #  legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=2:1)
    plotSurv(su, b, subset=b!="Others", ylab=toupper(j), xlab="Time (yrs)",
      legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c(ET.colors["SA4"], palette()[2]))   
    plotSurv(su, ifelse(e==8, "SP8", "Others"), ylab=toupper(j), xlab="Time (yrs)",
      legendPos=NULL, addTable=TRUE, leftTable=4, ppos=.2, col=c("black", ET.colors["SA8"]))   
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
  { plotSurv(cl[[what]], paste0("SA", e), col=ET.colors, addTable=TRUE, legendPos='bottomleft',
      ppos=.1, ylab=what, xlab="Time (yr)")
  }
  dev.off();
}

for (what in c("RFS", "iBCFS", "DRFS", "OS")) #colnames(cli)[sapply(cli, is, "Surv")])
{ pdf(paste0(figDir, "survEco1vAll.", what, ".pdf"), width=12, height=6);
  par(mfrow=c(2,5), mar=c(2.5,3,1.5,.5), mgp=c(1.5,.5,0))
  for (i in 1:max(ecot))
  { plotSurv(cli[[what]], ifelse(ecot==i, paste0("SA", i), "Other"), xlab="Time (yrs)", ylab=what,
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
t[t=="0"] = "''"; colnames(t) = rownames(t) = paste0("SA", 1:max(ecot));
diag(t) = paste0("bold('", diag(t), "')");
t[t=="bold('''')"] = "''"; 
cairo_pdf(paste0(figDir, "confMatEcot.pdf"), width=4, height=4)
par(mar=c(0,1.5,1.5,0))
plotTable(t, parseCells=TRUE,
  tableLabels=c("Real Spatial Archetype", "Recovered Spatial Archetype"))
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
    
    levels(e2) = paste0("SA", levels(e2));
    e = factor(ecot); levels(e) = paste0("SA", levels(e));
    
    if (!is.null(ccM)) { w = intersect(colnames(cc), colnames(ccM)); cc=cc[,w]; ccM=ccM[,w]; }
  
    co=NULL;
    if (what=="Genes")
    { o = order(-match(colnames(cc), geneList$symbol)); cc = cc[,o]; ccM = ccM[,o];
      co = colGenes[colnames(cc)];
      names(co) = colnames(cc) = colnames(ccM) = geneList[colnames(cc), "List.of.interesting.genes"]
    }
    if (what=="xCell")
    { cc = cc[,!grepl("Score", colnames(cc))];
      o = rev(intersect(names(colXct), colnames(cc))); cc = cc[,o]; ccM=ccM[,o];
      co = colXct[colnames(cc)]; names(co) = colnames(cc);
    }
    if (what=="Sigs")
    { cc = cc[,intersect(names(sigInfo), colnames(cc))];
      ccM = ccM[,intersect(names(sigInfo), colnames(ccM))];
      co = colSig[sigInfo[colnames(cc)]]; names(co) = colnames(cc);
    }
    
    if (is.null(ccM))
    { z = data.frame(ID=rownames(cli), Ecotype=e, cc, check.names=FALSE);
      write.xlsx2(z, file=paste0(figDataDir, "EcoComp.", what, ".xlsx"), row.names=FALSE);
    } else
    { z = data.frame(ID=c(rownames(cli), rownames(d$cli)),
        study=rep(c("ST",c(MB="METABRIC", Bulk="Spatial TNBC", ScanB="SCAN-B")[version]), c(nrow(cli), nrow(d$cli))),
        `Spatial Archetype`=c(e, e2), rbind(cc, ccM), check.names=FALSE);
      write.xlsx2(z, file=paste0(figDataDir, "EcoComp.", what, ".", version, ".xlsx"), row.names=FALSE);
    }
  
    p1 = calcP(cc, ecot); 
    p = rowMins(p1$fdr);  
    if (!is.null(ccM)) { p2 = calcP(ccM, e2); p = pmin(p, rowMins(p2$fdr))  }
    
    w = colnames(cc)[which(p<.05)];
    if (length(w)==0) { next; }
    
    z = p.adjust(p1$p[w,], method='fdr');
    cutOff = max(z[p1$fdr[w,]<.05])+1e-10;
    
    if (is.null(ccM))
    { cairo_pdf(paste0(figDir, "EcoComp.", what, ".pdf"),
        height=(length(w)+7)/7, width=6)
    } else 
    { cairo_pdf(paste0(figDir, "EcoComp.", what, ".", version, ".pdf"),
        height=(length(w)+7)/7, width=6)
    }
    par(mar=c(3,18,.5,.5), mgp=c(1.5,.5,0), cex=.5, omi=c(0,0,0,1));
    dotPlot(cc[,w,drop=FALSE], e, col.lbl=co[w], pch=c("full", "left")[(!is.null(ccM))+1], maxP=cutOff,
      legend=FALSE);
    if (!is.null(ccM))
    { z = p.adjust(p2$p[w,], method='fdr');
      cutOff = max(z[p2$fdr[w,]<=.05])+1e-10;
      dotPlot(ccM[,w,drop=FALSE], e2, pch='right', add=TRUE, maxP=cutOff, legend=FALSE);
      legendDotPlot(x=par('usr')[2]+strwidth("M")*6, y=par('usr')[4]-strheight("J")*1, interline=2,
        double=TRUE, doubleAnnot=c("ST PB", c(ScanB="SCAN-B", MB="METABRIC", Bulk="ST Bulk")[version]),
        significantLabel=c("FDR < 5%", "FDR ≥ 5%"))
    } else
    { legendDotPlot(x=par('usr')[2]+strwidth("M")*6, y=par('usr')[4]-strheight("J")*3, interline=2,
        significantLabel=c("FDR < 5%", "FDR ≥ 5%")) }

    dev.off()
  }
}

# ET by ... heatmaps
#####################
SA = split(seq_along(ecot), ecot);
ecot2 = factor(paste0("SA", ecot));
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
  
  z = data.frame(ID=rownames(cli), `Spatial Archetype`=paste0("SA", ecot), x, check.names=FALSE);
  write.xlsx2(z, file=paste0(figDataDir, "ET.HM.", i, ".xlsx"), row.names=FALSE);
  
  p1 = calcP(x, ecot); 
  p = rowMins(p1$fdr);  
  
  w = colnames(x)[which(p<.05)];
  x = x[,w];
  
  z = p.adjust(p1$p[w,], method='fdr');
  cutOff = max(z[p1$fdr[w,]<.05])+1e-10;
  
  cairo_pdf(paste0(figDir, "ET.HM.", i, ".pdf"),
      height=1.2*((length(w)+7)/7+1), width=5+8+max(nchar(w))/10)
  par(mar=c(6,max(nchar(w)+4)/2.5,.5,.5), mgp=c(1.5,.5,0));
  layout(matrix(1:2, nrow=1), width=c(.4,.6));
  
  dotPlot(x, ecot2, col.lbl=annots[colnames(x)], maxP=cutOff, legend=FALSE, cex.pch=.45);
  legendDotPlot(par("usr")[1]-2*strwidth("M"),-3*strheight("M"), horizontal=TRUE, cex.pch=.45, interline=1.6,
    significantLabel=c("FDR < 5%", "FDR ≥ 5%"))

  if (i=="xc") { x = log(x+.01);  }

  a = t(x[unlist(SA),]);
  a = a-rowMins(a); a=a/rowMaxs(a);

  par(mar=c(6,max(nchar(rownames(a))+4)/2.5,.5,.5));

  image(t(a), xaxt='n', yaxt='n', col=colorRampPalette(c('blue', 'yellow'))(100), xaxs='i', yaxs='i')
  mtext(rownames(a), at=((1:nrow(a))-1)/(nrow(a)-1), side=2, line=.5,
      col=annots[rownames(a)], las=2)
  ww = cumsum(sapply(SA, length));
  abline(v = (ww[-9]-.5)/(ncol(a)-1), col='white', lwd=2)
  mtext(paste0("SA", 1:9), side=1, at=((ww+c(0, ww[-9]))/2-.5)/(ncol(a)-1),line=0.5)
  plotScale(colorRampPalette(c('blue', 'yellow'))(100), c('Low', 'High'), atV=c(0,1), horizontal=TRUE,
    width=20, height=2, posx=0, posy=par('usr')[3]-6*strheight("M"))
  dev.off();
}

## ET single plot
####################
#toDisp = read.xlsx(paste0(dataDir, "misc/ET-characteristics.xlsx"), 2, header=FALSE)[,1:2]
toDisp = read.xlsx(paste0(dataDir, "misc/Selected molecular feature Fig 5D.xlsx"), 1)
toDisp[,3] = sub(" +$", "", toDisp[,3])
toDisp = toDisp[!is.na(toDisp[,3]),];

#toDisp[,2] = sub(" gene$", "", toDisp[,2]);
toDisp[,3] = sub("Myc", "MYC", toDisp[,3]);
toDisp[which(toDisp[,3]=="VCpredTN"),3] = "VCpred TN";

at = cumsum(1+!duplicated(toDisp[,1]))-1;
ETnames = read.xlsx(paste0(dataDir, "misc/ET-characteristics.xlsx"), 1, header=TRUE)[,1]
ETnames = ETnames[!is.na(ETnames)];
ETnames = sub(" (.)$", "\\1", ETnames);
#w = which(!is.na(toDisp[,1])); cl = rep(toDisp[w,1],c(w[-1], nrow(toDisp)+1)-w-1)
cl = toDisp[,1];

wg = intersect(rownames(PB), intersect(geneList$symbol, rownames(bulk)))
x = t(PB[wg,]); colnames(x) = geneList[wg,"List.of.interesting.genes"];
cc = cbind(csPB, x, xcPB[,1:48])
pi = apply(cc, 2, function(i) tapply(1:nrow(cc), ecot, function(s) wilcox.test(i[s], i[-s])$p.value))
pi[] = p.adjust(pi, method='fdr')
#w = rev(toDisp[!is.na(toDisp[,2]),2])
w = toDisp[,3];
z = cc[,w]; a = colnames(z) %in% geneList[wg,"List.of.interesting.genes"];
#colnames(z)[a] = paste(colnames(z)[a], "(gene)");
lbls = paste0("'", colnames(z), "'"); lbls[a] = paste0("italic(", lbls[a], ")"); 
#lbls=sub("^(.+) \\(gene\\)$", "paste(italic('\\1'), ' gene')", colnames(z));
#ww=grepl("^paste", lbls); lbls[!ww]=paste0("'", lbls[!ww], "'");
lbls = parse(text=lbls);
wxc = colnames(z) %in% colnames(xcPB);
colLbls = c('grey45', 'black', 'darkblue'); names(colLbls) = c("Cell types", "Signatures", "Single genes");
#colLbls=colLbls[c(2,1,3)];
col = rep(colLbls["Signatures"], length(a));
col[a] = colLbls["Single genes"]; col[wxc] = colLbls["Cell types"];
#o = order(factor(rev(cl), levels=unique(rev(cl))), factor(col,levels=colLbls[c(1,3,2)]));
#lbls=lbls[o]; z=z[,o]; col=col[o];

#cairo_pdf(paste0(figDir, "EcoSingle.horiz.pdf"), width=8.6, height=2.5) # Version with ET description
#par(mar=c(1,18,9,3), oma=c(0,0,0,4), mgp=c(1.5,.5,0), cex=.5, xpd=NA);
cairo_pdf(paste0(figDir, "EcoSingle.horiz.pdf"), width=7, height=2.5)
par(mar=c(1,3,9,3), oma=c(0,0,0,4), mgp=c(1.5,.5,0), cex=.5, xpd=NA);
ETnames = paste0("SA", 1:9);
dotPlot(z, factor(ETnames[ecot], levels=rev(ETnames)), toDisp=t(pi[nrow(pi):1,w])<.05, at=at,
  horizontal=TRUE, axPos=2, srt=30, inMa=1, lbls=lbls, col.lbl=col, legend=FALSE)
#wi = (which(!is.na(toDisp[,1]))-1)[-1];
wi = at[which(!duplicated(cl))[-1]]
rect(wi-1.9, rep(9.7,4), wi-.1, rep(10.5, 4), col='white', xpd=NA, border=NA);
wi = c(1, wi, at[nrow(toDisp)]+2);
mtext(txt<-unique(cl), side=1, line=-.3, at=(wi[-1]+wi[-length(wi)])/2-1, cex=.5)
for (i in 1:(length(wi)-1))
{ lines(wi[(0:1)+i]+c(.5,-.5)-1, c(.3,.3), xpd=NA, col=colDisp[txt[i]], lwd=2)
}
wi2 = setdiff(at[which(toDisp[-1,2] != toDisp[-nrow(toDisp),2])+1], wi)
rect(wi2-.8, rep(9.7,4), wi2-.2, rep(10.5, 4), col='white', xpd=NA, border=NA);

legendDotPlot(par('usr')[2]+2, par('usr')[4]+1, significantLabel=c("FDR < 5%", "FDR ≥ 5%"), interline=1.5)
#legend(par('usr')[2]+2, par('usr')[4]-6, names(colLbls)[c(2,3,1)], fill=colLbls[c(2,3,1)], bty='n', xjust=.5,
#  text.col=colLbls[c(2,3,1)], border=NA)
dev.off();

x = data.frame(ID=rownames(cli), `Spatial Archetype`=paste0("SA", ecot), cc);
write.xlsx2(x, file=paste0(figDataDir, "EcoSingle.horiz.xlsx"), row.names=FALSE)

## Targets
###############
tgt = c("NECTIN4", "ERBB2", "TACSTD2", "CD274", "AR")
names(tgt)=tgt; names(tgt)[3:4] = c("Trop-2", "PD-L1");

e = ecotRec
x = do.call(cbind, lapply(names(tgt), function(nm)
{ unlist(lapply(ds, function(i) scale(i$dn[tgt[nm],])));
})); colnames(x) = tgt;
hm = do.call(rbind, tapply(1:nrow(x), e, function(i) colMeans(x[i,])));
pi = do.call(rbind, tapply(1:nrow(x), e, function(i) apply(x,2,function(j) wilcox.test(j[i], j[-i], alt='g')$p.value)));
colnames(hm) = colnames(pi) = names(tgt);
qi = pi; qi[] = p.adjust(pi, method='fdr')

hm2 = (hm+1)/2; hm2[hm2<0] = 0; hm2[hm2>1] = 1;
hm2 = t(hm2); hm2 = hm2[5:1,9:1];
qi2 = t(qi)[5:1,9:1];
pdf(paste0(figDir, "hmTargets.pdf"), height=2.1, width=1.4)
par(mar=c(1+2,3,9-4,.5), mgp=c(1.5,.5,0), cex=.5)#, xpd=NA);
image(hm2, xaxt='n', yaxt='n', col=colorRampPalette(c('blue', 'yellow'))(100))
mtext(paste0("SA", colnames(hm2)), at=(0:8)/8, side=2, line=.5, las=2, cex=par('cex')) 
mtext(rownames(hm2), at=(0:4)/4, side=3, line=.2, las=2, font=3, cex=par('cex')) 

w = which(qi2<1e-5, arr.ind=TRUE);
points(((0:4)/4)[w[,1]], ((0:8)/8)[w[,2]], pch=16, col='black')
par(xpd=NA);
plotScale(colorRampPalette(c('blue', 'yellow'))(100), c("-1 s.d.", "0", "+1 s.d."), c(0, .5, 1),
  posx=par('usr')[1]-0*strwidth("M"), posy=par('usr')[3]-3*strheight("M"), width=13, height=1, title=NULL, title.font=2,
  horizontal=TRUE)
dev.off()

x = do.call(rbind, lapply(ds, function(i) scale(t(i$dn[tgt,]))))
x = data.frame(ID=unlist(lapply(ds, function(i) rownames(i$cli))),
  Study=rep(names(ds), sapply(ds, function(i) nrow(i$cli))), `Spatial Archetype`=paste0("SA", ecotRec),
  x, check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "hmTargets.xlsx"), row.names=FALSE);

## Alluvials ETs
#################
f = function(i) factor(i, levels=c("IM", "BL", "M", "MSL", "LAR"));
#coET = colEcot; names(coET)=paste0("SA", 1:9);
coET = ET.colors;
pdf(paste0(figDir, "alluvialET.pdf"));
par(mfrow=c(1,1), mar=c(3,3,0.5,.5), mgp=c(1.5,.5,0), xpd=NA)
t = plyr::count(data.frame(Subtype=f(cli$barPB), `Spatial Archetype`=factor(paste0("SA",ecot)),
  TIME=factor(cli$TIME_classes.bypathologist, levels=rev(c('ID', 'MR', 'SR', 'FI'))), check.names=FALSE));
t = t[complete.cases(t),];
colnames(t)[2]="Spatial Archetype"
alluvial(t[,1:3], freq=t$freq, col=colBar[as.character(t$Subtype)], blockCol=c(colBar,  colTIME, coET),
  draw_ticks=FALSE)
dev.off();

pdf(paste0(figDir, "alluvialET2.pdf"));
par(mfrow=c(1,1), mar=c(3,3,0.5,.5), mgp=c(1.5,.5,0), xpd=NA)
t = plyr::count(data.frame(Tumor=f(cli$barT.an), Ecotype=factor(paste0("ET",ecot)),
  Stroma=f(cli$barS.an)));
t = t[complete.cases(t),];
alluvial(t[,1:3], freq=t$freq, col=colBar[as.character(t$Tumor)], blockCol=c(colBar,  colTIME, coET),
  draw_ticks=FALSE)
dev.off();

n = do.call(cbind, tapply(seq_along(ecot), ecot, function(i) colMeans(fa[i,])))
#n = 100*n/rep(colSums(n), each=nrow(n));
colnames(n) = paste0("ET", colnames(n));
rownames(n) = paste0("MC", rownames(n));
n2 = data.frame(Ecotype=rep(colnames(n), each=nrow(n)), Megacluster=factor(rownames(n), levels=paste0("MC", 1:14)),
  freq=as.vector(n));
n2 = n2[n2$freq > 1e-5,];

pdf(paste0(figDir, "alluvialMCvsET.pdf"), height=5, width=6)
par(mar=c(2,3,1,5), mgp=c(1.5,.5,0));
alluvial(n2[,1:2], freq=n2$freq, col=coET[as.character(n2$Ecotype)], blockCol=c(MC.colors, coET),
  draw_ticks=FALSE)
dev.off();

## Source data for other panels
###################################
x = data.frame(ID=rownames(cli), `Spatial Archetype`=paste0("SA", ecot),
	`TLS ST signature`=csPB[,"TLS ST"], check.names=FALSE);
write.xlsx2(x, file=paste0(figDataDir, "TLSsigVsET.xlsx"), row.names=FALSE);


## Merge source data
######################
or = list(
  fig3=c(b="morphoTotal", c="morphoDistUp", e="morphoPatch.Tumor"),
  fig4=c(b="alluvialDec.2", c="tumorCompartment", d="stromaCompartment",
    e="TvsS.Mbp.hm", f="survTSdec.final"),
  fig5=c(b="tlsViolin", c="tlsGO.VSlympho.blueOnly", d="tlsGene", f="TLSsigVsBar", g="TLSsigVsTIME"),
  fig6=c(a="tlsAllCliSmallA", b="tlsAllCliSmallB", c="tlsAllCliSmallC", d="tlsAllCliSmallD",
    e="controlTLSby_immunoT.all", f="controlTLSby_immunoT.all"),
  fig7=c(c="barByClust", d="tblsMC2", e="forestMC.PB", f="forestMC.All2"),
  fig8=c(a="heatmapMC", b="barplotMCvsET", c="barplotETs", d="TLSsigVsET",
    e="EcoSingle.horiz", f="hmTargets", g="forestET.All3", h="survEcotypeNew2"),
  figS1=c(a="morphoDistUp",b="morphoPatch.Stroma"),
  figS2=c(a="aucs"),
  figS3=c(a1="TvsS.bar.dec.Tumor.Annot.horiz", a2="TvsS.bar.dec.Tumor.Genes.horiz", a3="TvsS.bar.dec.Tumor.Sigs.horiz",
    b1="TvsS.bar.dec.Stroma.Annot.horiz", b2="TvsS.bar.dec.Stroma.xCell.horiz",
    b3="TvsS.bar.dec.Stroma.Genes.horiz", b4="TvsS.bar.dec.Stroma.Sigs.horiz"),
  figS4=c(a="patchSizeVsSig.Tumor", b="patchSizeVsSig.Stroma"),
  figS5=c(a="tlsXcl", b="tlsGO.VSlympho", c="tlsGO.VSrest"),
  figS6=c(a="aucTLS2", b='aucTLS', c="compTLSsigs", d="TLSbyCompartment"),
  figS7=c(ab="tlsIspy", c='TLSimmunoOS', d='fpTLSimmuno'),
  figS8=c(ab="controlTLS_TNBC"),
  figS9=c(ab="controlTLSispy", cdef="controlTLSby_immunoT.all.SI"),
  figS11=c(a="meanAUC", bc="MCcomp.Sigs"),
  figS12=c(ab="MCcomp.Genes"),
  figS13=c(ab="MCcomp.xCell"),
  figS14=c(ab="forestMC.PB"),
  figS15=c(ab="forestMC.Bulk"),
  figS16=c(ab="forestMC.MB"),
  figS17=c(ab="forestMC.SCANB"),
  figS18=c(ab="forestMC.All3"),
  figS19=c(a="N_ET", b="barEcot", c="EcoComp.Annot", d="EcoComp.Cli"),
  figS20=c(ac="ET.HM.xc", bd="ET.HM.sig"),
  figS21=c(ab="ET.HM.genes"),
  figS22=c(a="EcoComp.Sigs.Bulk", b="EcoComp.Sigs.Bulk"),
  figS23=c(a="EcoComp.Sigs.ScanB"),
  figS24=c(a="EcoComp.xCell.Bulk"),
  figS25=c(a="EcoComp.xCell.MB"),
  figS26=c(a="EcoComp.xCell.ScanB"),
  figS27=c(a="EcoComp.Genes.Bulk", b="EcoComp.Genes.MB"),
  figS28=c(a="EcoComp.Genes.ScanB"),
  figS29=c(ab="forestET.PB"),
  figS30=c(ab="forestET.Bulk"),
  figS31=c(ab="forestET.MB"),
  figS32=c(ab="forestET.SCANB"),
  figS33=c(ab="forestET.All3")
);

names(or) = sub("figS", "Supplementary Fig. ", names(or))
names(or) = sub("^fig", "Fig. ", names(or))
for (f in names(or)[6])
{ wb = createWorkbook();
  csT = CellStyle(wb, font=Font(wb, isBold=TRUE))
  for (i in names(or[[f]]))
  { fi = paste0(figDataDir, or[[f]][[i]],".xlsx");
    if (!file.exists(fi)) { warning("No file", fi); next; } 
    x = read.xlsx2(fi, sheetIndex=1, check.names=FALSE);
    sheet = createSheet(wb, paste0(f, i));
    addDataFrame(x, sheet, row.names=FALSE, colnamesStyle=csT)
  }
  saveWorkbook(wb, file=paste0(figDataDir, "final/", f, ".xlsx"))
}