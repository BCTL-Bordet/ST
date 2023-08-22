\name{plotScale}
\title{Add a scale to an existing plot}
\description{
	Add a color scale to an existing plot.
}
\usage{
plotScale(cols, v, atV, posx, posy, width=1, height=10, title=NULL, title.font=2,
  horizontal=FALSE)}
\arguments{
	\item{cols}{The color scheme for the scale}
	\item{v}{The values to show on the side of the scale}
	\item{atV}{At which ticks should the values be represented? (between 0 and 1)}
	\item{posx, posy}{The position of bottom left corner the scale in user coordinates}
	\item{width}{The width of the scale in terms of character width.}
	\item{height}{The height of the scale in terms of character height.}
	\item{title}{The title of the scale (if not NULL)}
	\item{title.font}{The font used to display the title}
}

\seealso{\link{nonSTstuff}
}