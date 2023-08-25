\name{formatNiceP}
\alias{formatNiceP}
\title{Format a p-value for display}
\description{
	Format a p-value for display, rounded to a given number of significant digits and using
	prettier notation for the exponent. Note that ending 0s will be kept to fit the number of 
	digits requested.
}
\usage{
formatNiceP(x, descr=NULL, parse=TRUE, lim=1e-20, pName="p", sepDescr=" - ",
  boldDescr=FALSE)
}
\arguments{
	\item{x}{The number to format}
	\item{descr}{An optional description text to be appended before the p-value.}
	\item{parse}{If \code{TRUE}, output an expression if needed. If \code{FALSE}, output a character
	  string that can be parsed into an expression (may be useful for column names etc). }
	\item{lim}{Values below \code{lim} are put at 0.}
	\item{pName}{The name of the p-value-like number (e.g. p, fdr, q).}
	\item{sepDescr}{Separation between \code{descr} and the p-value.}
	\item{boldDescr}{Whether the \code{descr} part should be in bold.}
}
\value{ Either a character string (if \code{parse==FALSE} or there is no exponent), or an
  expression (if \code{parse==TRUE} and there is an exponent). }

\seealso{\link{nonSTstuff}, \code{\link{formatNice}}
}