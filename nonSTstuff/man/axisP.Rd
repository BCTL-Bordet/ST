\name{axisP}
\alias{axisP}
\title{Display a modified log axis}
\description{
	Often data are naturally on a log scale, but there are some zeros that prevent them to be displayed.
	In those cases, it is natural to display \verb{log(x+e)} where \verb{x} is the value to
	display and \verb{e} is a constant. 
	\verb{axisP} will display the axis corresponding to this case.
}
\usage{
axisP(Plog, side=2, nice=FALSE, maxVal=NULL)
}
\arguments{
	\item{Plog}{The value added.}
	\item{side}{On which side to display the axis.}
	\item{nice}{If \code{TRUE}, shows the axis label in scientific notation.}
	\item{maxVal}{Maximum value to display on the axis (to truncate the axis).}
}
\seealso{\link{nonSTstuff}, \code{\link{myBoxplot}}, \code{\link{plotWcc}}
}