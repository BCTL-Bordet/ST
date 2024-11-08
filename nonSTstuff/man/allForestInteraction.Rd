\name{allForestInteraction}
\alias{allForestInteraction}
\title{Display all interaction tests between dataframe of predictors, a binary variable and an output as a forest plot.}
\description{
	Displays all interaction tests between dataframe of predictors, and a binary variable and an output as a forest plot.
	For instance, this can be useful to test whether a feature has significantly different hazard ratios in two
	treatment arms.
	The output can be either binary (in which case logistic regressions will be used for the
	determination of the confidence intervals and potentially \verb{wilcox.test} or \verb{fisher.test}
	for the p-values) or
	a survival object (in which case Cox regression will be used for the CI and potentially
	\verb{survdiff} for the p-value).
}
\usage{
allForestInteraction = function(x, y, inter = "", subset=NULL,
				control=~1, sig=.05,
				full=TRUE, fullBar=full,
        clip=c(.2,5), ticks=c(.2,.5,1,2,5), 
        lbl1="", lbl2="inter", genesItal=FALSE, 
        weights=NULL, boxsize= 0.3, names=NULL,
        graphwidth = unit(4, "cm"), plot=TRUE, new_page=TRUE, ...}
\arguments{
	\item{x}{The dataframe or matrix of predictors to test.}
	\item{y}{The outcome variable, as a 2-levels factor or a survival object.}
	\item{iter}{The column name from \code{x}  meant to be used as the interaction variable.}
	\item{subset}{The subset of cases to take into account (if not \verb{NULL}).}
	\item{control}{Features that are used as control variables in a multivariable setting, as a right-hand side
		formula (e.g. \code{~ x + y} if \code{x} and \code{y} are the covariate to control for). }
	\item{sig}{Significance threshold for coloring the confidence intervals.}
	\item{full}{If TRUE, columns shown include the global HR/OR for the feature of interest. If FALSE,
		only the interaction is shown.}
	\item{fullBar}{If TRUE, the CIs for both arms are plotted. Otherwise, the CI on the interaction estimate is plotted.}
	\item{clip}{Limits on the OR/HR x-axis.}
  \item{ticks}{Where to put the ticks on the x-axis.}
	\item{lbl1}{How the p-values, HR/OR and CI column names are prefixed for the global effect.}
	\item{lbl2}{How the p-values, HR/OR and CI column names are prefixed for the interaction.}
  \item{genesItal}{If \verb{TRUE}, the feature names that fit gene names are put in italics.}
  \item{weights}{Case weights to be used for the analysis. In the case of Cox regression, will lead to 
    use the robust se from \verb{coxph}, and is incompatible with Firth method and nested models.}
  \item{boxsize}{The size of the boxes for the forest plot.}
  \item{names}{Vector or list of names to display instead of using the columns from x. Can be used
  	for complicated characters, or to use mathematical notation.}
  \item{graphwidth}{The width of the plot part.}
  \item{new_page}{Should the plot be done on a new page?}
  \item{plot}{Whether the plot should be made. Otherwise the returned grob can be used.}
  \item{...}{Additional parameters passed to \verb{forestplot}. }
}

\value{
   A list containing the following
    components:
    \item{tt}{The table to be displayed on the forest plot, as a list of lists.}
    \item{mean, lower, upper}{The confidence intervals (lower to upper) and values (mean)
    	to display in the forest plot, as 2-column matrices.}
    \item{p}{The p-values as numbers.}
    \item{grob}{tTe grob that can be plotted.}
}

\seealso{\link{nonSTstuff}, \code{\link{plotTable}}, \code{\link{allForest}}
}