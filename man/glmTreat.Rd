\name{glmTreat}
\alias{glmTreat}

\title{Test for Differential Expression Relative to a Threshold}

\description{Conduct genewise statistical tests for a given coefficient or contrast relative to a specified fold-change threshold.}

\usage{
glmTreat(glmfit, coef = ncol(glmfit$design), contrast = NULL, lfc = 0, null = "interval")
}

\arguments{
\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmFit} or \code{glmQLFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.  Values must be columns or column names of \code{design}. Defaults to the last coefficient. Ignored if \code{contrast} is specified.}

\item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested against the log2-fold-change threshold. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}

\item{lfc}{numeric scalar specifying the absolute value of the log2-fold change threshold above which differential expression is to be considered.}

\item{null}{character string, choices are \code{"worst.case"} or \code{"interval"}.
If \code{"worst.case"}, then the null hypothesis asssumes that the true logFC is on the boundary of the possible values, either at \code{lfc} or \code{-lfc}, whichever gives the largest p-value.
This gives the most conservative results.
If \code{"interval"}, then the null hypotheses assumes the true logFC to belong to a bounded interval of possible values.}
}

\value{
\code{glmTreat} produces an object of class \code{DGELRT} with the same components as for \code{glmfit} plus the following:
	\item{lfc}{absolute value of the specified log2-fold-change threshold.}
	\item{table}{data frame with the same rows as \code{glmfit} containing the log2-fold-changes, average log2-counts per million and p-values, ready to be displayed by \code{topTags}.}
	\item{comparison}{character string describing the coefficient or the contrast being tested.}

The data frame \code{table} contains the following columns:
	\item{logFC}{shrunk log2-fold-change of expression between conditions being tested.}
	\item{unshrunk.logFC}{unshrunk log2-fold-change of expression between conditions being tested. Exists only when \code{prior.count} is not equal to 0 for \code{glmfit}.}
	\item{logCPM}{average log2-counts per million, the average taken over all libraries.}
	\item{PValue}{p-values.}
}

\details{
\code{glmTreat} implements a test for differential expression relative to a minimum required fold-change threshold.
Instead of testing for genes which have log-fold-changes different from zero, it tests whether the log2-fold-change is greater than \code{lfc} in absolute value.
\code{glmTreat} is analogous to the TREAT approach developed by McCarthy and Smyth (2009) for microarrays.

\code{glmTreat} detects whether \code{glmfit} was produced by \code{glmFit} or \code{glmQLFit}.
In the former case, it conducts a modified likelihood ratio test (LRT) against the fold-change threshold.
In the latter case, it conducts a quasi-likelihood (QL) F-test against the threshold.

If \code{lfc=0}, then \code{glmTreat} is equivalent to \code{glmLRT} or \code{glmQLFTest}, depending on whether likelihood or quasi-likelihood is being used.

If there is no shrinkage on log-fold-changes, i.e., fitting glms with \code{prior.count=0}, then \code{unshrunk.logFC} and \code{logFC} are essentially the same. Hence they are merged into one column of \code{logFC} in \code{table}.
Note that \code{glmTreat} constructs test statistics using \code{unshrunk.logFC} rather than \code{logFC}.

\code{glmTreat} with positive \code{lfc} gives larger p-values than would be obtained with \code{lfc=0}.
If \code{null="worst.case"}, then \code{glmTreat} conducts a test closely analogous to the \code{treat} function in the limma package.
This conducts a test if which the null hypothesis puts the true logFC on the boundary of the \code{[-lfc,lfc]} interval closest to the observed logFC.
If \code{null="interval"}, then the null hypotheses assumes an interval of possible values for the true logFC.
This approach is somewhat less conservative.
}

\note{
\code{glmTreat} was previously called \code{treatDGE}.
}

\author{Yunshun Chen and Gordon Smyth}

\seealso{
\code{\link{topTags}} displays results from \code{glmTreat}.

\code{\link[limma]{treat}} is the corresponding function in the limma package, designed on normally distributed expression data rather than negative binomial counts.
}

\references{
McCarthy, D. J., and Smyth, G. K. (2009).
Testing significance relative to a fold-change threshold is a TREAT.
\emph{Bioinformatics} 25, 765-771.
\url{http://bioinformatics.oxfordjournals.org/content/25/6/765}
}

\examples{
ngenes <- 100
n1 <- 3
n2 <- 3
nlibs <- n1+n2
mu <- 100
phi <- 0.1
group <- c(rep(1,n1), rep(2,n2))
design <- model.matrix(~as.factor(group))

### 4-fold change for the first 5 genes
i <- 1:5
fc <- 4
mu <- matrix(mu, ngenes, nlibs)
mu[i, 1:n1] <- mu[i, 1:n1]*fc

counts <- matrix(rnbinom(ngenes*nlibs, mu=mu, size=1/phi), ngenes, nlibs)
d <- DGEList(counts=counts,lib.size=rep(1e6, nlibs), group=group)

gfit <- glmFit(d, design, dispersion=phi)
tr <- glmTreat(gfit, coef=2, lfc=1)
topTags(tr)
}
