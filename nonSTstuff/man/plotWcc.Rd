\name{plotWcc}
\alias{plotWcc}
\title{Scatter plot with correlation}
\description{
	Does a regular scatter plot, but add the correlation, CI on the correlation, p-values etc
	as requested. Also can be used for \verb{log(x+a)} type of axis.
}
\usage{
plotWcc(x, y, sub='', method='s', addP=FALSE, addCI=FALSE, subSide=1,
  subLine=par('mgp')[1]+1, pWhat="cor", atP=NA,
  cexP=par('cex'), Plog, digits=0, ...)
}
\arguments{
	\item{x, y}{The two variables to display on the scatter plot}
	\item{sub}{Some text to add to the left of the line with the correlation}
	\item{method}{The correlation method used (as in \verb{cor})}
	\item{addP, addCI}{Should the p-value (\verb{cor.test}) or confidence intervals be displayed?}
	\item{subSide}{On which side to display the correlation}
	\item{subLine}{Line on which to display the correlation (as in \verb{mtext})}
	\item{pWhat}{The meaning of the p-value displayed. If \verb{cor}, the p-value is related to the correlation.
	  If \verb{diff}, the p-value is a paired Wilcoxon test between x and y}
	\item{atP}{Where on the x-axis to display the correlation (default to centered)}
	\item{cexP}{Character expansion factor for the correlation information}
	\item{Plog}{If not \code{FALSE}, the value to add to the data prior to plotting in a log-like axis
	  with zeros (see \code{\link{axisP}}).
	  If a single scalar, then the same value is used for both axis, otherwise the 1st value is for
	  the x-axis and the 2nd for the y-axis. Implies that the axis will be log-like.}
	\item{digits}{Number of significant digits to display for the correlation}
	\item{...}{Further parameters passed to \verb{plot}}
}

\seealso{\link{nonSTstuff}, \link{axisP}
}