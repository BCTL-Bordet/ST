\name{addAnnot}
\alias{addAnnot}
\title{Add an annotation line to an existing plot}
\description{
	Adds an annotation line made of colored rectangles to an existing plot, in the plot margins.
	Useful to annotate heatmaps or dendrograms for instance.
}
\usage{
addAnnot(name, cols, line, side=1, at=seq_along(cols), heights=1, pch=NA,
  point.col=NA, point.cex=1)
}
\arguments{
	\item{name}{Name of the annotation (displayed next to it)}
	\item{cols}{Colors of the annotation}
	\item{line}{Line on which to draw the annotation line (as in mtext). }
	\item{side}{Side on which to draw the annotation (1: bottom, 2: left, 3:top, 4:right).}
	\item{at}{Positions of the items to be annotated, that will be the center of the boxes.}
	\item{heights}{Height of the annotation line, in terms of lines (as in mtext).}
	\item{pch}{If not NA, will add points (with the specified pch) in the middle of the rectangles for
	  a potential 2nd level of annotation.}
	\item{point.col}{Colors of the points (if present). }
	\item{point.cex}{Size of the points (if present). }
}

\seealso{\link{nonSTstuff}
}