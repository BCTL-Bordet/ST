\name{formatNice}
\alias{formatNice}
\title{Format a number for display}
\description{
	Format a number for display, round to a given number of significant digits and using
	prettier notation for the exponent. Note that ending 0s will be kept to fit the number of 
	digits requested.
}
\usage{
formatNice(x, parse=TRUE, lim=1e-20, nDigits=2)
}
\arguments{
	\item{x}{The number to format}
	\item{parse}{If TRUE, output an expression if needed. If FALSE, output a character
	  string that can be parsed into an expression (may be useful for column names etc). }
	\item{lim}{Values below lim are put at 0.}
	\item{nDigits}{Number of significant digits.}
}
\value{ Either a character string (if parse is FALSE or there is no exponent), or an
  expression (if parse is TRUE and there is an exponent). }

\seealso{\link{STstuff}, \link{formatNiceP}
}