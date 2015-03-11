\name{glmTreat}
\alias{glmTreat}

\title{Testing for Differential Expression Relative to a Threshold}

\description{Conduct genewise statistical tests for a given coefficient or coefficient contrast relative to a specified threshold.}

\usage{
glmTreat(glmfit, coef=ncol(glmfit$design), contrast=NULL, lfc=0)
}

\arguments{
\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmFit} or \code{glmQLFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.  Values must be columns or column names of \code{design}. Defaults to the last coefficient.  Ignored if \code{contrast} is specified.}

\item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested against the log2-fold change threshold. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}

\item{lfc}{numeric scalar specifying the absolute value of the log2-fold change threshold above which differential expression is to be considered.}
}

\value{
\code{glmTreat} produces an object of class \code{DGELRT} with the same components as for \code{glmfit} plus the following:
	\item{lfc}{absolute value of the specified log2-fold change threshold.}
	\item{table}{data frame with the same rows as \code{glmfit} containing the log2-fold changes, average log2-counts per million and p-values, ready to be displayed by \code{topTags.}.}
	\item{comparison}{character string describing the coefficient or the contrast being tested.}

The data frame \code{table} contains the following columns:
	\item{logFC}{log2-fold change of expression between conditions being tested.}
	\item{logCPM}{average log2-counts per million, the average taken over all libraries.}
	\item{PValue}{p-values.}
}

\details{
\code{glmTreat} implements a two-sided modified likelihood ratio test.
}

\author{Yunshun Chen and Gordon Smyth}

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
