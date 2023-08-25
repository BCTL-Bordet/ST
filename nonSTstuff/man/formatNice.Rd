\name{formatNice}
\alias{formatNice}
\title{Format a number for display}
\description{
	Format a number for display, round to a given number of significant digits and use
	prettier notation for the exponent. Note that ending 0s will be kept to fit the number of 
	digits requested.
}
\usage{
formatNice(x, parse=TRUE, lim=1e-20, nDigits=2)
}
\arguments{
	\item{x}{The number to format}
	\item{parse}{If \code{TRUE}, output an expression if needed. If \code{FALSE}, output a character
	  string that can be parsed into an expression (may be useful for column names that can be
	  parsed later, for instance with \code{\link{allForest}}). }
	\item{lim}{Values below \code{lim} are put at \code{0}.}
	\item{nDigits}{Number of significant digits.}
}
\value{ Either a character string (if \code{parse==FALSE} or there is no exponent), or an
  expression (if \code{parse==TRUE} and there is an exponent). }

\seealso{\link{nonSTstuff}, \code{\link{formatNiceP}}
}