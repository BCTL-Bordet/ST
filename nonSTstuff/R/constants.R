colAnn2 = c(Nothing="#FFFFFF", Tumor="#017801", Necrosis="#000000", `Fat tissue`="#000080",
  `Low TIL stroma`="#ff904f", Vessels="#dc0000", Artefacts="#6e2400",
  `Lactiferous duct`="#9980e6", `High TIL stroma`="#e9e900", `in situ`="#ccffcc",
  `Lymphoid nodule`="#80801a", `Hole (whitespace)`="#40e5f6", Lymphocyte="#c4417f",
  `Stroma cell`="#ff9980", Nerve="#4d8080", `Heterologous elements`="#808080",
  `Acellular stroma`="#e9d1bb", `Tumor region`="#258a15")

colsBY=c('lightblue', 'yellow');

colTIME = colorRampPalette(c('grey', 'red'), space="Lab")(6)[-c(3,5)];
names(colTIME) = c('ID', 'MR', 'SR', 'FI')

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