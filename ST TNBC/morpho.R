source("~/prog/ST/ST TNBC/STscripts.R")

nAn = seq_along(colAnn2); names(nAn) = names(colAnn2);

fs = dir(paste0(dataDir, "annotationRecoded/"));

# Quantified patches
nbs = mclapply(fs, function(f)
{ message(f);
  if (file.exists(paste0(dataDir, "patches/", f))) { return("Already there"); }
  im = readRDS(paste0(dataDir, "annotationRecoded/", f));
  
  Np = tabulate(im, nbins=length(nAn)); names(Np) = names(nAn);
  
  patches = list();
  # patches...
  for (what in c("Tumor", "Stroma", "Lymphocytes", "Fat", "Necrosis", "in situ", "Ducts"))
  { a = Image( im %in% nAn[list(Tumor=c("Tumor","Tumor region"), Stroma=c("Low TIL stroma", "Vessels",
      "High TIL stroma", "Stroma cell"), Lymphocytes=c("High TIL stroma", "Lymphocyte"),
      Fat="Fat tissue", Necrosis="Necrosis", `in situ`="in situ", Ducts="Lactiferous duct")[[what]]], dim=dim(im));
    # Remove very small
    if (sum(a)==0) { next; }
    lbl = bwlabel(a);
    if (what %in% c("Ducts")) { lbl = fillHull(lbl); a = lbl>0; }
    t = tabulate(lbl); #table(lbl)[-1];
    a[lbl %in% as.integer(names(t)[t<=20])] = 0;
    b = a; if (any(t>1e4)) { b[lbl %in% as.integer(names(t)[t>=1e4])] = 0; } # Do not dilate large ones
    b = dilate(b, makeBrush(15, shape='disc'));
    a = a | b | (im==nAn["Tumor region"]);
    lbl = bwlabel(a); rm(a,b); gc();
    patches[[what]] = tabulate(lbl); # Size of each patch
    rm(lbl); gc();
  }
 
  saveRDS(list(Np = Np, patches=patches), file=paste0(dataDir, "patches/", f));
  return("ok");
}, mc.preschedule=FALSE);

# Count each annotation
pat = lapply(d<-dir(paste0(dataDir, "patches")),
  function(i) readRDS(paste0(dataDir, "patches/",i)))
names(pat) = sub("TNBC([0-9]+).RDS", "\\1", d);
pat = pat[rownames(cli)]

N = do.call(rbind, lapply(pat, function(i)
{ i$Np
}));

cli = readRDS(paste0(dataDir, "Clinical/Clinical.RDS"))
identical(N, cli$annotations) # Should be TRUE