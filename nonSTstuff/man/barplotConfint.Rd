\name{barplotConfint}
\alias{barplotConfint}
\title{barplot with automatic confidence intervals}
\description{
	Display a barplot with confidence intervals. Those are computed using \verb{binom.test}.
}
\usage{
barplotConfint(a, b, length=.1, addN=FALSE, addP=FALSE, lineP=3, sideP=1, p, ...)
}
\arguments{
	\item{a}{Either a table (with 2 rows or 2 columns) or the number of TRUE (if b is missing)}
	\item{b}{Missing if a is a table. Otherwise, the total number of cases.}
	\item{length}{Length of the bars at the top of the CIs.}
	\item{addN}{Should the total number of cases be displayed?}
	\item{addP}{Should the total p-value (Fisher test) be displayed?}
	\item{lineP}{If addP is TRUE, gives the line (\link{title}) on which the p-value is displayed. }
	\item{sideP}{If addP is TRUE, gives the side on which the p-value should be displayed.}
	\item{p}{The p-value to be displayed (default to a Fisher test).}
	\item{...}{Further parameters passed to \verb{barplot}.}					
}
\seealso{\link{nonSTstuff}
}