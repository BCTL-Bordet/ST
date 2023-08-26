\name{legendDotPlot}
\alias{legendDotPlot}
\title{Legend for dotplots}
\description{
	Legend for dotplots.
}
\usage{
legendDotPlot(x, y, double=FALSE, doubleAnnot, interline=1,
  family="Arial Unicode MS", pt.lwd=.5, horizontal=FALSE, cex.pch=1)
  }
\arguments{
	\item{x, y}{Position of the legend}
	\item{double}{If \code{TRUE}, legend is for double analysis (left and right dots)}
	\item{doubleAnnot}{Names of the double annotations}
	\item{interline}{Interline between the dots}
	\item{family}{Font family}
	\item{pt.lwd}{Width of the black outline of the dots}
	\item{horizontal}{Whether the legend should be vertical or horizontal}
	\item{cex.pch}{General scaling factor the dots}
}

\seealso{\link{nonSTstuff}, \link{dotPlot}
}