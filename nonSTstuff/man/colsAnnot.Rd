\name{colsAnnot}
\alias{colsAnnot}
\title{Colors for typical annotations}
\description{
	Give the color scheme for biological annotations. Mostly those are breast-cancer related.
}
\usage{
colsAnnot(what=c("pam50", "IHC", "IntSub", "SCM", "leh", "all"))
}
\arguments{
	\item{what}{Gives the colors chosen for a set of biological annotations.
	  pam50 is the PAM50 classification of breast cancers.
	  IHC is the 3-classes classification based on HR and HER2.
	  IntSub is the same as IHC with different names.
	  SCM is the 3 genes classifier + prolif.
	  leh is the 5 classes triple negative breast cancer classification.
	   }
}

\seealso{\link{STstuff}
}