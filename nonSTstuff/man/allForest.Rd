\name{allForest}
\alias{allForest}
\title{Display all tests between dataframe of predictors and an output as a forest plot.}
\description{
	Displays all tests between dataframe of predictors and an output as a forest plot.
	The output can be either binary (in which case logistic regressions will be used for the
	determination of the confidence intervals and potentially \verb{wilcox.test} or \verb{fisher.test}
	for the p-values) or
	a survival object (in which case Cox regression will be used for the CI and potentially
	\verb{survdiff} for the p-value).
}
\usage{
allForest = function(x, y, subset=NULL, control=~1, sig=.05, sigOnFDR=FALSE, fdr=TRUE,
  new_page=FALSE, clip=c(.2,5), ticks=NULL, dispN=FALSE, useWilcox=FALSE, genesItal=FALSE,
  parseNames=FALSE, rmNAp=TRUE, normReal=TRUE,
  capHR=TRUE, p=NULL, weights=NULL, useNested=TRUE, columns=c("OR", "CI", "P", "FDR"),
  boxsize=.2, version=3, FirthIfNeeded=TRUE, graphwidth = unit(4, "cm"), plot=TRUE, ...)}
\arguments{
	\item{x}{The dataframe or matrix of predictors to test.}
	\item{y}{The outcome variable, as a 2-levels factor or a survival object.}
	\item{subset}{The subset of cases to take into account (if not \verb{NULL}).}
	\item{control}{Features that are used as control variables in a multivariable setting.}
	\item{sig}{Significance threshold for coloring the confidence intervals.}
  \item{sigOnFDR}{If \verb{TRUE}, the threshold \verb{sig} for coloring applies to FDR, otherwise to p-values.}
  \item{fdr}{Should the FDR be reported?}
  \item{new_page}{For version 1 to 2, should the plot be done on a new page?}
  \item{clip}{Limits on the OR/HR.}
  \item{ticks}{Where to put the ticks on the x-axis.}
  \item{dispN}{Should the number of cases for each feature (i.e. the number of non-NA values) be reported?}
  \item{useWilcox}{If \verb{TRUE}, use \verb{wilcox.test}/\verb{fisher.test} for categorical outcome,
    \verb{survdiff} for survival outcome.}
  \item{genesItal}{If \verb{TRUE}, the feature names that fit gene names are put in italics.}
  \item{parseNames}{If \verb{TRUE}, parse the column names from dataframe \verb{x} for display, allowing to use
    all of \verb{plotmath} possibilities.}
  \item{rmNAp}{If \verb{TRUE}, features on which p-values cannot be evaluated are not displayed.}
  \item{normReal}{If \verb{TRUE}, scale the numerical features in \verb{x} if \verb{x} is a dataframe.
    This will impact the OR/HR but not the p-values. }
  \item{capHR}{If \verb{TRUE}, shows the OR/HR that are outside of clip as triangles on the limits of the plot. }
  \item{weights}{Case weights to be used for the analysis. In the case of Cox regression, will lead to 
    use the robust se from \verb{coxph}, and is incompatible with Firth method and nested models.}
  \item{useNested}{If \verb{TRUE}, used nested models + log likelihood test from \verb{anova} to assess significance.}
  \item{columns}{The set of columns to display.}
  \item{boxsize}{The size of the boxes for the forest plot.}
  \item{version}{1, 2 or 3. 1 and 2 are different layout using \verb{forestplot} and so grid graphics.
      3 is similar to 2, but uses \verb{basicForest} which uses base graphics.}
  \item{FirthIfNeeded}{Use Firth penalized likelihood method in case the CI are too extreme. P-values
    are still from the base models.}
  \item{graphwidth}{The width of the plot part, in case of version 1 or 2.}
  \item{plot}{If \verb{TRUE}, make the plot (for version 1 or 2). }
  \item{...}{Additional parameters passed to \verb{forestplot} (version 1 or 2) or \code{\link{basicForest}} (version 3). }
}

\value{
   A list containing the following
    components:
    \item{tt}{The table to be displayed on the forest plot, as a list of lists.}
    \item{tn}{The confidence intervals and values to display in the forest plot, as a matrix with
      columns mean, lower and upper.}
    \item{p}{The p-values as numbers.}
    \item{grob}{If version is 1 or 2, the grob that can be plotted.}
}

\seealso{\link{nonSTstuff}, \code{\link{plotTable}}, \code{\link{basicForest}}
}