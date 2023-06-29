\name{basicForest}
\alias{basicForest}
\title{Forestplot in base graphics}
\description{
	Display a basic forestplot in base graphics
}
\usage{
basicForest(x, a, adj=NULL, xlim=NULL, xlog=FALSE, xlab="", col='black', cex.axis=.7,
  titles=names(x), annotDir=NULL, lineHeight=1.5, cexAnnot=0.7, colWidth=NULL)}
\arguments{
	\item{x}{The table to be displayed, as a matrix, dataframe or list of list (necessary if expressions are to be used).}
	\item{a}{The numbers for the forestplot, in the form of min / mean / max. }
	\item{adj}{Column adjustment (see par). Default is left adjustment for the first column and center adjustment for the others.}
	\item{xlim}{The limits on the graph part.}
	\item{xlog}{Should the values be displayed in logs?}
  \item{xlab}{Label of the x axis.}
  \item{col}{Colors of the dots (recycled).}
  \item{cex.axis}{Cex for the axis of the plot.}
  \item{titles}{Titles of the columns.}
  \item{annotDir}{Text of lenght 2 to annotate the meaning of the direction on the plot (e.g. worse outcome / better outcome). }
  \item{lineHeight}{Height of the lines (recycled).}
  \item{cexAnnot}{Cex for the annotDir}
  \item{colWidth}{Width of the columns. Default is based on the max width of each column.}
}

\seealso{\link{STstuff}, \link{plotTable}
}