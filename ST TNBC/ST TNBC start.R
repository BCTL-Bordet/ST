source("~/prog/ST/ST TNBC/STscripts.R")
library(TeachingDemos)
library(xlsx)

#load("~/Data/Spatial/TNBC/res/start2.RData")
# source('~/prog/ngs/Studies/ST TNBC start.R')

###############################
# Generate starting data
##############################
  
# Signatures
sigs = readRDS(paste0(dataDir, "misc/signatures.RDS"))
sigs = sigs[c("ESR1", "CIN70", "GGI", "GENE70", "Immune1", "Immune2", "Stroma1", "Stroma2",
  "AR_gene", "TLS_Lundeberg", "VCpred_TN", "Parpi7", "Parpi7.norm")]
names(sigs)[names(sigs)=="ESR1"] = "ESR1 gene";
names(sigs) = sub("_", " ", names(sigs));

si = readRDS(paste0(dataDir, "misc/hallmarks.RDS"))
names(si) = tolower(sub("HALLMARK_", "h_", names(si)));
si$h_spermatogenesis = si$h_pancreas_beta_cells = si$h_uv_response = NULL;
nsi = strsplit(names(si), "_");
for (i in seq_along(nsi))
{ w= nsi[[i]] %in% c("e2f", "il2", "stat5", "il6", "jak", "stat3", "kras", "mtorc1", "myc", "notch", "pi3k",
    "akt",  "mtor", "tnfa", "nfkb", "wnt", "dna", "g2m", "tgf");
  nsi[[i]][w] = toupper(nsi[[i]][w]);
}
names(si) = sapply(nsi, function(i) paste("H", paste(i[-1], collapse=" ")))
substr(names(si),3,3) = toupper(substr(names(si),3,3))
names(si) = sub("MTOR", "mTOR", names(si)); names(si) = sub("P53", "p53", names(si));
names(si) = sub("NFKB", "NF-kB", names(si))
names(si) = sub("^H ", "", names(si))

x = readRDS(paste0(dataDir, "misc/signatures.RDS"))[c("CAF_29455927", "Trm_36026440", "TAM")];
names(x) = sub("_.+", "", names(x));

y = readRDS(paste0(dataDir, "misc/signature ST.RDS"))

sigH = c(si,x, y);
sigH = lapply(sigH, function(i) as.character(i$name));
  
load(paste0(dataDir, "misc/sigInfo.RData")) # sigInfo, colSig
colSig[5]="darkred";

# 1. Clinical etc
###################
ids = readRDS(paste0(dataDir, "Clinical/ids.RDS"))
cli = readRDS(paste0(dataDir, "Clinical/Clinical.RDS"))

cli$TIME_classes.bypathologist[cli$TIME_classes.bypathologist=="nd"] = NA;

ctrl = cli[,c("Age", "LN.status", "T_TNM")]
ctrl[,2]=ctrl[,2]>0; ctrl[,3]=ctrl[,3]>1;
colnames(ctrl) = c("Age", "LN", "Size");
ctrlForm = ~Age+LN+Size;

for (i in which(sapply(cli, is, "Surv")))
{ a = cli[[i]]; w = a[,1]>10; a[w,1]=10; a[w,2]=0; cli[[i]]=a;
}

# Patches
##########
pat = lapply(d<-dir(paste0(dataDir, "patches")),
  function(i) readRDS(paste0(dataDir, "patches/",i)))
names(pat) = sub("TNBC([0-9]+).RDS", "\\1", d);
pat = pat[rownames(cli)]

N = sapply(pat, function(i)
{ sapply(list(Tumor=c("Tumor","Tumor region"), Stroma=c("Low TIL stroma", "Vessels",
      "High TIL stroma", "Stroma cell"), Lymphocytes=c("High TIL stroma", "Lymphocyte"),
      Fat="Fat tissue", Necrosis="Necrosis", `in situ`="in situ", Ducts="Lactiferous duct"),
    function(w) sum(i$N[w]));
});

pp = lapply(rownames(N), function(w) sapply(pat, function(x) cPar(x$patches[[w]])))
names(pp) = rownames(N); for (i in 2:length(pp)) { rownames(pp[[i]]) = rownames(pp[[1]]); }
for (i in names(pp)) { cli[[paste0("patches.", i)]] = t(pp[[i]]); }

# Bulk RNA seq
################
geneMap = readRDS(paste0(dataDir, "misc/geneMap.RDS"))
bulk = readRDS(paste0(dataDir, "bulk_count.RDS"))
bulk = rowsum(bulk, geneMap$bulk[rownames(bulk)]); # Map to the genes of the ST
bulk = normRNAseq(bulk, lim=1e3);

cli = cli[colnames(bulk),]

PBraw = readRDS(paste0(dataDir, "PB_count.RDS"))
PBraw = rowsum(PBraw, geneMap$PB[rownames(PBraw)]);
#PBraw = t(rowsum(t(PBbsR), ids[colnames(PBbsR), "id"]))
PBraw = PBraw[,rownames(cli)];
PB = normRNAseq(PBraw, lim=1e3)
g = intersect(rownames(bulk), rownames(PB))

cli$bar = factBareche(TNBCclassif(bulk, version="bareche", shortName=TRUE, coef=FALSE))
cli$barco = TNBCclassif(bulk, version="bareche", shortName=TRUE, coef=TRUE)
cli$leh = TNBCclassif(bulk, version="lehmann", shortName=TRUE, coef=FALSE)
cli$lehco = TNBCclassif(bulk, version="lehmann", shortName=TRUE, coef=TRUE)
cli$barPB = factBareche(TNBCclassif(PB, version="bareche", shortName=TRUE, coef=FALSE))
cli$barPBco = TNBCclassif(PB, version="bareche", shortName=TRUE, coef=TRUE)
cli$lehPB = TNBCclassif(PB, version="lehmann", shortName=TRUE, coef=FALSE)
cli$lehPBco = TNBCclassif(PB, version="lehmann", shortName=TRUE, coef=TRUE)
cli$bur = TNBCclassif(bulk, version='burstein', shortName=TRUE)
cli$bur2 = TNBCclassif(bulk, version='burstein2', shortName=TRUE)
cli$burPB = TNBCclassif(PB, version='burstein', shortName=TRUE)
cli$bur2PB = TNBCclassif(PB, version='burstein2', shortName=TRUE)
cli$burPBco = TNBCclassif(PB, version='burstein', shortName=TRUE, coef=TRUE)
cli$bur2PBco = TNBCclassif(PB, version='burstein2', shortName=TRUE, coef=TRUE)

#cli$TIME = calcTIME(bulk); cli$TIMEpb = calcTIME(PB);

 
xcPB = t(microenvironmentScores(xCellAnalysis(PB, cell.types.use=xCellPoss)))
xcBulk = t(microenvironmentScores(xCellAnalysis(bulk, cell.types.use=xCellPoss)))

csPB = calcAllSig(PB); 
csBulk = calcAllSig(bulk);

N = cli$annotations;
N = N[,!(colnames(N) %in% c("Nothing", "Artefacts", "Hole (whitespace)"))]

n = 100*N/rowSums(N)
n = cbind(n, Stroma=NA, `Total stroma`=NA);
n[,"Tumor"] = n[,"Tumor"] + n[,"Tumor region"];
n[,"Stroma"] = n[,"Stroma cell"] + n[,"Low TIL stroma"]/2 + n[,"High TIL stroma"]/2
n[,"Acellular stroma"] = n[,"Acellular stroma"] + n[,"Low TIL stroma"]/2;
n[,"Lymphocyte"] = n[,"Lymphocyte"] + n[,"High TIL stroma"]/2
n[,"Total stroma"] = n[,"Stroma"] + n[,"Acellular stroma"];
n = n[, c("Tumor", "Stroma", "Acellular stroma", "Total stroma", "Lymphocyte", 'in situ',
  "Fat tissue", "Necrosis", "Lactiferous duct", "Vessels", "Lymphoid nodule", "Heterologous elements", "Nerve")];
cli$annotClean=n;


# Stuff for plotting
#######################

# Genes of interest
x = read.xlsx(file=paste0(dataDir, "misc/list of interesting genes.xlsx"), 1)
x$List.of.interesting.genes = x$New.list.of.interested.gene.names
x$List.of.interesting.genes = sub(" +$", "", x$List.of.interesting.genes)
x$List.of.interesting.genes[x$List.of.interesting.genes == "mTOR"] = "MTOR";
x$List.of.interesting.genes[x$List.of.interesting.genes == "TOP1A"] = "TOP1";
x = x[x$List.of.interesting.genes != "TOP1B",];
x$symbol = entrez2Symb(symb2Entrez(x$List.of.interesting.genes, tryHard=TRUE))
x$symbol[x$List.of.interesting.genes=="Trop-2"] = "TACSTD2";
geneList = x;
rownames(geneList) = geneList$symbol;
co = unique(geneList[, "Role.of.genes"]);
colGenesPal = palette.colors(palette="Polychrome 36")[-2][seq_along(co)]; names(colGenesPal) = co;
colGenes = colGenesPal[geneList[, "Role.of.genes"]]; names(colGenes) = geneList$symbol;

# Groups and colors for xCell
xcD=read.xlsx(paste0(dataDir, "misc/xCells - redefinition.xlsx"), 1, check.names=FALSE);
xcD = xcD[!grepl("Abbreviations", names(xcD))]
a = rep(names(xcD), each=nrow(xcD)); names(a) = unlist(xcD);
a = a[!is.na(names(a))];
colXc = c(`Lymphocytes B`="#FFB730" ,`Lymphocytes T CD4+`="#00A74E" ,`Lymphocytes T CD8+`="#00A7E9" ,
  `Other lymphocytes`="#C74E1F" ,`Myeloid cells`="#FF0017" ,`Stromal cells`="#BB8322" ,
  `Stem cells`="#6A2B91" ,`Other cells`="#A9BFE1")
colXct = colXc[a]; names(colXct) = names(a);

# Short clinical
ccn = cli[,c("Age", "Histology.invasive.tumor", "Number.of.tumor.foci", "Tumor.size_millimeters",
  "Grade", "N_TNM", "Ki67_percentage", "BRCA", "sTILs.percentage.bypathologist", "CD3.percentage.bypathologist",
  "CD20.percentage.by pathologist")];
ccn$medullary = ccn$Histology.invasive.tumor=="medullary"; ccn$Histology.invasive.tumor=NULL;
ccn$Number.of.tumor.foci = ccn$Number.of.tumor.foci==1
ccn$Grade[ccn$Grade=="nd"] = NA;
for (i in which(sapply(ccn, is.numeric))) { ccn[,i] = scale(ccn[,i]); }
ccn = sapply(ccn, as.numeric)
cliShort=ccn;

# x. Deconvolution wrt classification by spot
##############################################
allS = mclapply(d<-dir(paste0(dataDir, "deconvolution")), function(i)
{ load(paste0(dataDir, "deconvolution/", i))
  list(prot=rowsum(prot, geneMap$PB[rownames(prot)]), Ntot=Ntot);
}); names(allS) = sub("TNBC([0-9]+).RData", "\\1", d); allS=allS[rownames(cli)];
Ntot = unlist(lapply(allS, function(i) i$Ntot))
allS = lapply(allS, function(i) i$prot);

g = table(unlist(lapply(allS, function(i) rownames(i)))); g = names(g)[g>=max(g)/5]
xCanRaw = fullMat3(allS, g);
idCan = cbind(what=colnames(xCanRaw), id=rep(names(allS), sapply(allS, ncol)))
colnames(xCanRaw) = paste(idCan[,1], idCan[,2]);
Nt = matrix(NA, nrow=length(unique(idCan[,2])), ncol=length(unique(idCan[,1])));
rownames(Nt) = rownames(cli); colnames(Nt) = unique(idCan[,1]);
Nt[idCan[,c(2,1)]] = Ntot;

xCan = normIt(xCanRaw);
idAn = matrix(NA, nrow=length(unique(idCan[,2])), ncol=length(unique(idCan[,1])));
dimnames(idAn) = dimnames(Nt);
idAn[idCan[,c(2,1)]] = 1:nrow(idCan);
csAn = calcAllSig(xCan);
xcAn = t(microenvironmentScores(xCellAnalysis(xCan, cell.types.use=xCellPoss)))
goAn = gsva(xCan, lapply(readRDS(paste0(dataDir, "misc/MSigDB_c5_GO.RDS")),
  function(i) i$name), min.sz=5, parallel.sz=4L)
  
# Tumor / Stroma PB
#######################
PBn2 = normIt(PBraw)

w = setdiff(colnames(idAn), c("Tumor", "Necrosis", "in situ"))#, "Lactiferous duct"))
stromaAn = sapply(1:nrow(idAn), function(i)
{ rowSums(xCanRaw[,idAn[i,w]]*rep(Nt[i,w], each=nrow(xCanRaw)), na.rm=TRUE) })
stromaAn = normIt(stromaAn)
i = intersect(rownames(stromaAn), rownames(PBn2));
stromaAn2 = (stromaAn[i,]-rowMeans(PBn2[i,]))/rowSds(PBn2[i,]);
cli$barS.an = TNBCclassif(stromaAn2, version='bareche',shortName=TRUE, rescale=FALSE)

#w = c("Tumor", "in situ")
#tumorAn = sapply(1:nrow(idAn), function(i)
#{ rowSums(xCan[,idAn[i,w]]*rep(Nt[i,w], each=nrow(xCan)), na.rm=TRUE) })
tumorAn = xCanRaw[,idAn[,"Tumor"]]
tumorAn = normIt(tumorAn)
i = intersect(rownames(tumorAn), rownames(PBn2));
tumorAn2 = (tumorAn[i,]-rowMeans(PBn2[i,]))/rowSds(PBn2[i,]);
cli$barT.an = TNBCclassif(tumorAn2, version='bareche',shortName=TRUE, rescale=FALSE)

x = t(xCellAnalysis(cbind(tumorAn, stromaAn), cell.types.use=xCellPoss))
xcAnT = x[1:94,]; xcAnS = x[95:(2*94),]
csAnT = calcAllSig(tumorAn);
csAnS = calcAllSig(stromaAn);

########################
# Other datasets (METABRIC and SCAN-B)
#######################################
ds = readRDS(file=paste0(dataDir, "external datasets/otherTNBC.RDS"))

# Clean survival SCAN-B
x = ds$`SCAN-B`$cli;
ws = which(sapply(x, is, "Surv"));
for (i in ws) # Death is an event
{ w = which(x$OS[,2]==1 & x[[i]][,2]==0 & x$OS[,1]-1<=x[[i]][,1])
  x[[i]][w,2] = 1; x[[i]][w,1] = x$OS[w,1];
}
for (i in ws)
{ w = x[[i]][,1]>10; x[[i]][w,1]=10; x[[i]][w,2]=0; # Censor at 10 yr
}
ds$`SCAN-B`$cli = x;
# Keep only follow up cohort in ScanB
x = ds$`SCAN-B`; w = x$cli$Follow.up.cohort=="TRUE"; 
x$cli = x$cli[w,]; x$dn = x$dn[,w]; x$d = x$d[,w]; x$cs = x$cs[w,]; x$ecot=x$ecot[w,]; x$xcl = x$xcl[w,];
x[["15"]]$H = x[["15"]]$H[w,];
# Rename survival in SCAN-B
z = colnames(x$cli);
z[z=="DRFi"]="DRFS"; z[z=="BCFi"] = "iBCFS"; colnames(x$cli) = z;
ds$`SCAN-B`=x;

# Clean survival METABRIC
x = ds$METABRIC$cli;
x$OS = Surv(x$T/365.25, x$Death)
x$DRFS = Surv(x$TDR/365.25, x$DR | x$Death)
x$EFS = Surv(pmin(x$TDR, x$TLR)/365.25, x$DR | x$LR | x$Death)
for (i in c("OS", "DRFS", "EFS"))
{ su = x[[i]]; su[su[,1]==0,] = NA;
  w = su[,1]>10; su[w,1]=10; su[w,2]=0; # Censor at 10 yr
  x[[i]] = su;
}
w = colnames(x)[sapply(x, is, "Surv")];
w = setdiff(w, c("OS", "DRFS", "EFS"))
x[w]=NULL;
ds$METABRIC$cli = x;

# Survival all 3 dataasets together
sus = list(iBCFS = rbind(ds[[1]]$cli$EFS, ds[[2]]$cli$iBCFS, ds[[3]]$cli$iBCFS),
  DRFS = rbind(ds[[1]]$cli$DRFS, ds[[2]]$cli$DRFS, ds[[3]]$cli$DRFS),
  OS = rbind(ds[[1]]$cli$OS, ds[[2]]$cli$OS, ds[[3]]$cli$OS))
sus = lapply(sus, function(su) # METABRIC then Spatial
{ su = Surv(su[,1], su[,2]);
  su[su[,1]==0,] = NA;
  w = su[,1]>10; su[w,1]=10; su[w,2]=0; # Censor at 10 yr
  return(su);
})
susDs = sus;

# Order survival
o = c("RFS", "EFS", "iBCFS", "DRFS", "OS");
for (i in names(ds))
{ z = ds[[i]]$cli; s = z[o[o %in% names(z)]];
  z = cbind(z[!sapply(z, is, "Surv")], s);
  ds[[i]]$cli=z; 
}

stud = rep(c("M", "S", "SC"), sapply(ds, function(i) nrow(i$cli)))
ecotRec = unlist(lapply(ds, function(i) i$ecot[,1])); # Ecotypes on the other studies

# To control for MV analysis
ctrlMB = data.frame(sapply(ds$METABRIC$cli[,c("LYMPH_NODES_EXAMINED_POSITIVE", "TUMOR_SIZE", "AGE_AT_DIAGNOSIS")], as.numeric))
ctrlMB[,1] = ctrlMB[,1]>0; ctrlMB[,2] = ctrlMB[,2]>20;
colnames(ctrlMB) = c("LN", "Size", "Age"); ctrlMB = ctrlMB[,colnames(ctrl)];

ctrlSB = ds[["SCAN-B"]]$cli[,c("Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..", "LN", "T.size")];
colnames(ctrlSB)=c('Age', "LN", "Size");

## iSpy2
############
dIspy = readRDS(paste0(dataDir, "external datasets/iSpy2Expression.RDS"));
cliIspy = readRDS(paste0(dataDir, "external datasets/iSpy2Clinical.RDS"));
csI = calcAllSig(dIspy);
cliIspy$wI = cliIspy$Receptor.Subtype=="TN"; 
armI =  cliIspy$Arm..short.name.
armI[!(armI %in% c("Pembro", "Ctr"))] = NA;
cliIspy$arm = armI;
cliIspy$pcr = factor(cliIspy$pCR, levels=c(0,1), labels=c("no pCR", "pCR"));
cliIspy[cliIspy$wI,"bar"] = TNBCclassif(dIspy[,cliIspy$wI], "bareche", shortName=TRUE);

## Immunotherapies
######################
a = dir(bd<-paste0(dataDir, "external datasets/Immunotherapies/"));
dsIm = lapply(a, function(i)
{ cl = readRDS(paste0(bd, i, "/Clinical.RDS"));
  d = readRDS(paste0(bd, i, "/Expression.RDS"));
  z = symb2Entrez(rownames(d));
  d2 = d;
  d2 = (d2-rowMeans(d2))/rowSds(d2)
  g = calcAllSig(d2);
  s = names(cl)[sapply(cl, is, "Surv")]
  return(list(g=g, cl=cl, s=s));
}); names(dsIm) = a; 

s = unique(unlist(sapply(dsIm, function(i) colnames(i$g))))
resp = factor(unlist(lapply(dsIm, function(i) i$cl$response)));
d = do.call(rbind, lapply(dsIm, function(i) t(fullMat(t(i$g), s, NA))))
d = d[,!(colnames(d) %in% c("AR_gene", "ESR1"))];
d = scale(d);
dIm = d;
ty = unlist(lapply(dsIm, function(i) i$cl$cancer_type))

st = paste(rep(names(dsIm), sapply(dsIm, function(i) nrow(i$g))), ty)

s = unique(unlist(lapply(dsIm, function(i) i$s)))
suI = do.call(data.frame, lapply(s, function(w) 
{ r = do.call(rbind, lapply(dsIm, function(i)
  { if (!is.null(i$cl[[w]])) { return(i$cl[[w]]); } else { return(Surv(rep(NA, nrow(i$cl)))) }
  } ))
  Surv(r[,1], r[,2]);
})); colnames(suI)=s;

## Survival
ps = lapply(colnames(suI), function(what)
{ s = suI[,what];
  p = t(apply(d, 2, function(i) getCoxP(coxph(s~i+ strata(st)), HR=TRUE))) #, subset=cancer=="Bladder")))
  p2 = t(apply(d, 2, function(i) coef(summary(coxph(s~d[,"TLS ST"]+i+ strata(st))))[, c('Pr(>|z|)', "coef")]))
  return(list(p, p2));
}); names(ps) = colnames(suI);
psSurv=ps;

## Reponse
p = t(apply(d, 2, function(i) coef(summary(glmer(resp ~ i+(1|st), family=binomial)))["i", c("Pr(>|z|)", "Estimate")]))
#p2 = t(apply(d[,colnames(d)!="TLS ST"], 2, function(i) coef(summary(glmer(resp ~ d[,"TLS ST"]+i+(1|st), family=binomial)))[2:3, c("Pr(>|z|)", "Estimate")]))
pResp = p;

# x. Tumor / stroma PB
#########################
PBn2 = normIt(PBraw)

w = setdiff(colnames(idAn), c("Tumor", "Necrosis", "in situ"))#, "Lactiferous duct"))
stromaAn = sapply(1:nrow(idAn), function(i)
{ rowSums(xCanRaw[,idAn[i,w]]*rep(Nt[i,w], each=nrow(xCanRaw)), na.rm=TRUE) })
stromaAn = normIt(stromaAn)
i = intersect(rownames(stromaAn), rownames(PBn2));
stromaAn2 = (stromaAn[i,]-rowMeans(PBn2[i,]))/rowSds(PBn2[i,]);
cli$barS.an = factBareche(TNBCclassif(stromaAn2, version='bareche',shortName=TRUE, rescale=FALSE))

tumorAn = xCanRaw[,idAn[,"Tumor"]]
tumorAn = normIt(tumorAn)
i = intersect(rownames(tumorAn), rownames(PBn2));
tumorAn2 = (tumorAn[i,]-rowMeans(PBn2[i,]))/rowSds(PBn2[i,]);
cli$barT.an = factBareche(TNBCclassif(tumorAn2, version='bareche',shortName=TRUE, rescale=FALSE))
rm(PBn2);

# 3. TLS
############
# Select only the samples with actual TLS
hasTLS = c(12,15,16,19,2,20,27,28,29,30,31,32,33,38,39,42,43,44,45,46,48,5,50,52,53,54,56,57,6,61,62,63,64,65,67,69,70,71,72,
  73,77,78,80,82,85,86,87,88,89,9,92,94)
cli$hasTLS = rownames(cli) %in% hasTLS;
# Calculate genes upregulated in TLS
tt = setdiff(colnames(idAn), "Lymphoid nodule"); # Comparisons to do

pinp = lapply(tt, function(z) # Not paired
{ a = xCan[,idAn[cli$hasTLS,"Lymphoid nodule"]]; b = xCan[, idAn[,z]];
  a = a[,!is.na(a[1,])]; b = b[,!is.na(b[1,])]
  p = unlist(mclapply(1:nrow(a), function(i) wilcox.test(a[i,], b[i,], alt='g')$p.value))
  m = rowMedians(a)-rowMedians(b);
  cbind(m=m, p=p);
}); names(pinp) = tt;

## Sigs from Charoentong
g = lapply(readRDS(paste0(dataDir, "misc/signatures.RDS")), function(i) i$name);
g = g[grep("Charoentong", names(g))];
csX = gsva(xCan, g)

# Enrichment..
psGo = lapply(a<-c("xCell", "GO", "Sigs"), function(what)
{ x = list(xCell=xcAn, GO=t(goAn), Sigs=csAn)[[what]]
  p = lapply(tt, function(z) # Not paired
  { a = x[idAn[cli$hasTLS,"Lymphoid nodule"],]; b = x[idAn[,z], ];
    a = a[,!is.na(a[1,])]; b = b[,!is.na(b[1,])]
    p = do.call(rbind, mclapply(1:ncol(a), function(i)
      c(wilcox.test(a[,i], b[,i], alt='g')$p.value, wilcox.test(a[,i], b[,i], alt='l')$p.value)))
    rownames(p)=colnames(a)
    return(p)
  });
  p = do.call(cbind, p); colnames(p) = paste(rep(tt, each=2), c("Higher", "Lower"))
  p = cbind(`Max other higher`=rowMaxs(p[,seq(1, ncol(p), by=2)][,tt!="Lymphocyte"]),
    `Max other lower`=rowMaxs(p[,seq(2, ncol(p), by=2)][,tt!="Lymphocyte"]), p)
  return(p)
}); names(psGo) = a;

# 2. Annotations & classification per spots
###############################################
ist = readRDS(paste0(dataDir, "classification/classifAll.RDS"))
ist = lapply(ist, function(ii)
{ ii = ii$pr;
  ii = ii[,c('Fat tissue','in situ','Lactiferous duct','Lymphoid nodule','Necrosis',
    "Tumor", "Stroma", "Lymphocyte", 'Vessels')];
})
ist = ist[rownames(cli)] # Classification per spot

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

# 3. Megaclusters
#####################
# Load prototypes
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

xCn.xc = t(xCellAnalysis(xCn, cell.types.use=xCellPoss))
xCn.cs = calcAllSig(xCn);

## Bareche by cluster
y2 = rowsum(exp(xCn)-1, geneMap$PB[rownames(xCn)]);

g = intersect(rownames(y2), rownames(PBraw));

pb = log(1e4*PBraw[g,]/rep(colSums(PBraw[g,]), each=nrow(PBraw[g,]))+1)
y3 = log(1e4*y2[g,]/rep(colSums(y2[g,]), each=nrow(y2[g,]))+1)

y3n = (y3-rowMeans(pb))/rowSds(pb);
tn =  TNBCclassif(y3n, version='bar', shortName=TRUE, rescale=FALSE)
tnCo = TNBCclassif(y3n, version='bar', shortName=TRUE, rescale=FALSE, coef=TRUE)
idC$bar = factor(tn, levels=levels(cli$bar));
idC$barCo = tnCo;

# Kmeans
kms = lapply(d<-dir(fd<-paste0(dataDir, "clustering/Kmeans MC/")), function(n) readRDS(paste0(fd, n)))
names(kms) = sub("km([0-9]+).RDS", "\\1", d); kms = kms[order(as.integer(names(kms)))];

y = 1e4*xC/rep(colSums(xC), each=nrow(xC))
w = rowMeans(y>1)>.05;
si = sigma[w,];
y = y[w,];
y = log(y+1);

tmp = assignClust(kms[["15"]], t(y), idC, annot, cli, minPatient=1, clLim=0)

o = c(1, 2, 5, 6, 3, 4, 7, 8, 11, 12, 9, 10, 13, 14); oo = o; oo[o]=1:14 # Reorder
tmp$k = oo[tmp$k]; tmp$tp = tmp$tp[o,]; tmp$Np = tmp$Np[o]; tmp$desc = tmp$desc[o];
K = tmp$k;
pres = calcPres(K, idC);

# Map to old order...
#oOld = match(oo, c(5, 2, 1, 11, 6, 8, 7, 3, 9, 4, 12, 10, 14, 13));

# Deconvolution MC per spot
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
#fa = t(sapply(m2s, function(i) colMeans(i>.2)));colnames(fa)=1:ncol(fa)
fa[fa<.01]=0

faN = fa; faN[faN<.01] = .01;

hc = hclust(dist(log10(100*(faN))), method='ward.D2');
a = cutree(hc, 9); m = unique(a[hc$order]); n = rep(NA, length(m)); n[m]=seq_along(m); 
ecot = n[a];
