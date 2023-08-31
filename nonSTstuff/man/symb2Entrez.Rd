\name{symb2Entrez}
\alias{symb2Entrez}
\title{Finds Entrez IDs corresponding to gene symbols}
\description{
	Finds Entrez IDs corresponding to gene symbols. It is based on the org.Hs.eg.db library.
	}
\usage{
symb2Entrez(d, tryHard=FALSE)
}
\arguments{
	\item{d}{Vector of gene symbols.}
	\item{tryHard}{If \code{FALSE}, will only use official gene symbols. If \code{TRUE},
	  will also search for alternate aliases for those that do not fit using the official symbols.}
}
\value{
A character vector of Entrez IDs. The symbols for which no match was found are set to \code{NA}.
}
\seealso{\link{nonSTstuff}, \code{\link{entrez2Symb}}, \code{\link{calcSig}}
}