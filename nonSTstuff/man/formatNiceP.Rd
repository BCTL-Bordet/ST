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
	\item{parse}{If TRUE, output an expression if needed. If FALSE, output a character
	  string that can be parsed into an expression (may be useful for column names etc). }
	\item{lim}{Values below lim are put at 0.}
	\item{pName}{The name of the p-value-like number (e.g. p, fdr, q).}
	\item{sepDescr}{Separation between descr and the p-value.}
	\item{boldDescr}{Whether the descr part should be in bold.}
}
\value{ Either a character string (if parse is FALSE or there is no exponent), or an
  expression (if parse is TRUE and there is an exponent). }

\seealso{\link{nonSTstuff}, \link{formatNice}
}