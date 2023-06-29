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
	\item{x}{The matrix / dataframe to be plotted.}
	\item{colHeader}{The background of the header cells.}
	\item{space}{Space around the table.}
	\item{colCell}{Background color of the cell lines (recycled).}
	\item{colSize}{Relative sizes of the columns. Default to columns of identical sizes.}
	\item{rowSize}{Relative sizes of the rowns. Default to rows of identical sizes.}
	\item{tableLabels}{If TRUE, annotate the two dimensions of the plot. Can also be a vector of length two.}
	\item{showZeros}{If FALSE, values that are 0 are shown as empty cells (and not 0).}
	\item{parseCells}{If TRUE, the cells are parsed with parse(text=...).}
}
\seealso{\link{STstuff}
}