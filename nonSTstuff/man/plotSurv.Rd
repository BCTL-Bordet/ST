\name{plotSurv}
\title{Plot a survival analysis.}
\description{
	Plot a survival analysis, sort of like survminer in base graphics.
}
\usage{
plotSurv(su, x, pval=TRUE, col=NULL, legend.title=NULL, subset=NULL, addMedian=FALSE, medianUnit="",
  addN=FALSE, hr=NULL, ppos=.5, pposx=par("usr")[2]*.95, legSep=" ", useName=TRUE,
  legendPos='topright', addTable=FALSE, marTable=0,
    leftTable=3, cexTable=par("cex"), test="survdiff", useCox=FALSE, parseCateg=FALSE, mark.time=TRUE,
    nperm=1e4, args.legend=NULL, moveTimeOne=0, ...}
\arguments{
	\item{su}{Survival object}
	\item{x}{Groups on which to the the analysis. If a matrix/data.frame, groups are created by
	  combining the features.}
	\item{col}{Colors of the survival curves}
	\item{legend.title}{Title of the legend}
	\item{subset}{Subset of cases to analyse.}
	\item{addMedian}{Show the median survival for each group in the legend.}
	\item{medianUnit}{Unit of the median}
	\item{addN}{Show the number of cases for each group in the legend.}
	\item{hr}{Should the HR be displayed next to the p-values?}
	\item{ppos}{y position of the p-value (0 to 1)}
	\item{pposx}{x position of the p-value (user unit)}
	\item{legSep}{Separator for combination of groups if x is a matrix.}
	\item{useName}{Use the names of the features or only their values if x is a matrix.}
	\item{legendPos}{Position of the legend, as in \link{legend}.}
	\item{addTable}{Display the table below the KM?}
	\item{marTable}{Add extra margin lines to the table? Useful if having multiple survival plots
	  with different number of categories and wanting to show all figs the same size.}
	\item{leftTable}{How much space between the left numbers and the labels for the table?}
	\item{cexTable}{Character size of the table.}
	\item{test}{Which test to use. Can be "survdiff", "cox" or "perm" (the permlogrank from clinfun,
	  useful for cases with few samples). Cox is the logrank test in coxph. }
	\item{useCox}{If TRUE, sets test at "cox".}
	\item{parseCateg}{If TRUE, parse the group names to expression for display.}
	\item{mark.time}{Should the ticks be displayed for censoring?}
	\item{nperm}{Number of permutations used for permlogrank.}
	\item{args.legend}{Further arguments to be passed to legend (list).}
	\item{moveTimeOne}{Time point on which to calculate the first column of the table.
	  Can be useful if strange things happen at t=0 (e.g. left truncated data). }
}

\seealso{\link{nonSTstuff}
}