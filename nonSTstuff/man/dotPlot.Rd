\name{dotPlot}
\alias{dotPlot}
\title{dotplot showing significant association results}
\description{
	dotplot showing significant association results.
}
\usage{
dotPlot(x, cl, efi=NULL, lbls=NULL, toDisp=NULL, maxRange=4, oma=c(0,0,0,5), inMa=.5,
  col.lbl=NULL, maxP=.1, horizontal=FALSE, pch=c("full", "left", "right"), add=FALSE,
  family="Arial Unicode MS", blackBorder=TRUE, at=NULL, axPos=1, srt=45, pt.lwd=.5,
  cex.pch=1, ...)
}
\arguments{
	\item{x}{Matrix / dataframe of features to compare.}
	\item{cl}{Factor with the classes to be compared.}
	\item{efi}{Matrix with pre-computed effect sizes. Otherwise will be estimated using logistic regressions.}
	\item{lbls}{Names of the features to be used for display. Default to colnames(x).}
	\item{toDisp}{Binary matrix of size ncol(x), length(levels(cl)) indicating which dots to plot.
	  If absent, determined using computed p-values and maxP.}
	\item{maxRange}{Effect sizes larger (in absolute value) that maxRange will be bound to maxRange.}
	\item{oma}{Outer margins, used to display the legend.}
	\item{inMa}{Innmer margin.}
	\item{col.lbl}{Color for the feature labels.}
	\item{maxP}{Cutoff on the spots FDR used to decide which dot to display.}
	\item{horizontal}{Whether the plot should be presented vertically (with the classes on the x-axis) or vertically
	  (with the classes on the y-axis).}
	\item{pch}{Wheteher the dots should be full dots, or only left/right half dots.}
	\item{add}{If TRUE, the plot is added to the previous plot. In this case, only the dots are displayed.}
	\item{family}{The font family used. For pch=left or right an Unicode family where the right characters are present is needed.}
	\item{blackBorder}{If TRUE, the dots are shown with a black outline.}
	\item{at}{The positions of the features on the figure. Can be used to add space between groups of related features.}
	\item{axPos}{On which side the axis (with the factor levels) should be drawn (1 or 3 if horizontal is FALSE, 2 or 4 if TRUE).}
	\item{srt}{Angle at which the feature names are displayed, if horizontal is TRUE.}
	\item{pt.lwd}{Width of the black outline of the spots.}
	\item{cex.pch}{General scaling factor for the spots.}
	\item{...}{Not used at this stage}
}

\seealso{\link{nonSTstuff}
}