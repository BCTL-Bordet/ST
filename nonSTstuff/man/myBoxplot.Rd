\name{myBoxplot}
\title{Draw a nice boxplot}
\description{
	Draw a nice boxplot which display p-values and other bells and whistles.
	In particular, global p-values (Kruskal-Wallis by defaut) are displayed, as well
	as basic between group comparisons as "star" notation.
}
\usage{
myBoxplot(x, f, subset=NULL, subset2=NULL, baseLab=NULL, index=NULL, ylim, subLine=NA, subSide=1, subAt=NA,
  colPoints='blue', Plog=NULL, alternative="two.sided", leg45=FALSE, pch=1, normGroup=NULL, compToAll=FALSE, isComp=FALSE, alpha=1,
  xlab="", pvalLim=1e-20, xaxt, sepBaselab=" - ", boldBaselab=FALSE, paramTest=FALSE, pairsToShow=NULL,
  log=ifelse(is.null(Plog) || is.na(Plog), "", "y"), ...)
}
\arguments{
	\item{x}{The set of numbers for the boxplot. If a matrix, then the columns are considered as
	  the different categories for the boxplot and f is not taken into account.}
	\item{f}{The categories for the values of x.}
	\item{subset}{The subset of values to display.}
	\item{subset2}{The subset of values used to calculate the p-values.}
	\item{baseLab}{An optional additional text to put next to the p-value.}
	\item{index}{}
	\item{ylim}{The range of the y-axis.}
	\item{subLine}{On which line should the p-value be displayed.}
	\item{subSide}{On which side?}
	\item{subAt}{If not NA, where on the x-axis should the p-value be displayed, in user unit.}
	\item{colPoints}{The color of the points.}
	\item{Plog}{If not null, the values displayed will be log(x+Plog), although the axe labels will
	  be corrected to fit the original values. Useful for log-like values that contains 0 (e.g. counts). }
	\item{alternative}{For pairwise comparisons, what kind of test should be performed.}
	\item{leg45}{Should the categories be displayed on the x-axis with a 45 deg rotation?}
	\item{pch}{The plotting character}
	\item{normGroup}{For pairwise comparison, what is the base group.}
	\item{compToAll}{Should the comparison be done 1 vs. all?}
	\item{isComp}{Should points be linked between the groups, if x is a matrix?}
	\item{alpha}{Points transparency, used if <1.}
	\item{xlab}{Label of the x-axis.}
	\item{pValLim}{P-values below this are put to 0.}
	\item{xaxt}{Should the x-axis be displayed?}
	\item{sepBaselab}{The separation between baseLab and the p-value.}
	\item{boldBaselab}{The font of the baseLab.}
	\item{paramTest}{Should non-parametric tests (kruskal.test and wilcox.test) be used, or
	  the parametric equivalent (aov and t.test)?}
	\item{pairsToShow}{}
	\item{log}{Should the y-axis be in log?}
	\item{...}{Additional parameters passed to boxplot}
}
\value{ The global p-value comparing the groups (typically the kruskal.test). }

\seealso{\link{nonSTstuff}, \link{allForest},
}