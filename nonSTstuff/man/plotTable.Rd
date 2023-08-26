\name{plotTable}
\alias{plotTable}
\title{Display a table}
\description{
	Display a table (matrix) in a figure
}
\usage{
plotTable(x, colHeader='grey80', space=.01, colCell=c('grey96', 'grey90'),
  colSize=NULL, rowSize=NULL, tableLabels=NULL, showZeros=TRUE, parseCells=FALSE)
}
\arguments{
	\item{x}{The matrix / dataframe / table to be plotted.}
	\item{colHeader}{The background of the header cells.}
	\item{space}{Space around the cells.}
	\item{colCell}{Background color of the cell lines (recycled).}
	\item{colSize}{Relative sizes of the columns. Default to columns of identical sizes.}
	\item{rowSize}{Relative sizes of the rowns. Default to rows of identical sizes.}
	\item{tableLabels}{If \code{TRUE} and \code{x} is a \code{table},
	  annotate the two dimensions of the plot using the labels from \code{x}.
	  Can also be a vector of length two in which case it is directly using to annotate the dimensions.}
	\item{showZeros}{If \code{FALSE}, values that are 0 are shown as empty cells (and not 0).}
	\item{parseCells}{If \code{TRUE}, the cells are parsed with \code{parse(text=...)}, allowing to use
	  \code{plotmath}.}
}
\seealso{\link{nonSTstuff}, \link{allForest}
}