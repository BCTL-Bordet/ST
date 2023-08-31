\name{entrez2Symb}
\alias{entrez2Symb}
\title{Finds the gene symbol corresponding to a set of Entrez IDs}
\description{
	Finds the gene symbol corresponding to a set of Entrez IDs. It is based on the org.Hs.eg.db library.
	}
\usage{
entrez2Symb(d)
}
\arguments{
	\item{d}{Vector of Entrez IDs.}
}
\value{
A character vector of gene symbols.
}
\seealso{\link{nonSTstuff}, \code{\link{symb2Entrez}}, \code{\link{calcSig}}
}