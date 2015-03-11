\name{glmQLFit}
\alias{glmQLFit}
\alias{glmQLFit.DGEList}
\alias{glmQLFit.default}
\alias{glmQLFTest}

\title{Quasi-likelihood methods with empirical Bayes shrinkage}

\description{Fit a quasi-likelihood negative binomial generalized log-linear model to count data.
Conduct genewise statistical tests for a given coefficient or coefficient contrast.}

\usage{
\method{glmQLFit}{DGEList}(y, design=NULL, dispersion=NULL, offset=NULL, abundance.trend=TRUE, 
        robust=FALSE, winsor.tail.p=c(0.05, 0.1), \dots)
\method{glmQLFit}{default}(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, 
        abundance.trend=TRUE, AveLogCPM=NULL, robust=FALSE, 
        winsor.tail.p=c(0.05, 0.1), \dots)
glmQLFTest(glmfit, coef=ncol(glmfit$design), contrast=NULL) 
}

\arguments{
\item{y}{a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the tagwise linear models.}

\item{dispersion}{numeric scalar or vector of negative binomial dispersions. Must provide if \code{y} is a matrix of counts or a \code{DGEList} object without any dispersion component. If \code{NULL} will be extracted from the \code{DGEList} object \code{y}, with order of precedence: trended dispersions, common dispersion, a constant value of 0.05.}

\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalor or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix. If \code{NULL} will be computed by \code{getOffset(y)}.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size)}. Defaults to \code{colSums(y)}.}

\item{abundance.trend}{logical, whether to allow an abundance-dependent trend when estimating the prior values for the quasi-likelihood multiplicative dispersion parameter.}

\item{AveLogCPM}{average log2-counts per million, the average taken over all libraries in \code{y}. If \code{NULL} will be computed by \code{aveLogCPM(y)}.}

\item{robust}{logical, whether to estimate the prior degrees of freedom robustly.}

\item{winsor.tail.p}{numeric vector of length 2 giving proportion to trim (Winsorize) from lower and upper tail of the distribution of genewise deviances when estimating the hyperparameters. Positive values produce robust empirical Bayes ignoring outlier small or large deviances. Only used when \code{robust=TRUE}.} 

\item{\dots}{other arguments are passed to \code{\link{glmFit}}.}

\item{glmfit}{a \code{DGEGLM} object, usually output from \code{qlmQLFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.}

\item{contrast}{numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero.}
}

\value{
\code{glmQLFit} produces an object of class \code{DGEGLM} with the same components as that produced by \code{\link{glmFit}}, plus:
	\item{df.residual.zeros}{a numeric vector containing the number of residual degrees of freedom for the GLM fit of each gene.}	
	\item{df.prior}{a numeric vector or scalar, containing the prior degrees of freedom.}
	\item{var.prior}{a numeric vector of scalar, containing the location of the prior distribution.}
	\item{var.post}{a numeric vector containing the shrunk quasi-likelihood dispersion for each gene.} 
If \code{robust=TRUE}, \code{df.prior} is a vector whereby the prior d.f. is reduced for any outlier genes. Otherwise, a scalar is returned with the same prior d.f. for all genes.
If \code{abundance.trend=FALSE}, \code{var.prior} is a vector specifying a abundance-dependent trend for the location. Otherwise, a single location is returned for allg eens.

\code{glmQFTest} produces objects of class \code{DGELRT} with the same components as that produced by \code{\link{glmLRT}}.
It also stores \code{df.total}, a numeric vector containing the total degrees of freedom for the F-test.
The data frame of results is stored as \code{table} and can be visualized with \code{\link{topTags}}. 
It contains the following columns:
	\item{logFC}{log2-fold change of expression between conditions being tested.}
	\item{logCPM}{average log2-counts per million, the average taken over all libraries in \code{y}.}
	\item{F}{F-statistics.}
	\item{PValue}{p-values.}
}

\details{
\code{glmQLFTest} implements the quasi-likelihood method of Lund et al (2012).
It behaves the same as \code{glmLRT} except that it replaces likelihood ratio tests with quasi-likelihood F-tests for coefficients in the linear model.
This function calls the limma function \code{\link{squeezeVar}} to conduct empirical Bayes smoothing of the genewise multiplicative dispersions.
Note that the \code{QuasiSeq} package provides a alternative implementation of Lund et al (2012), with slightly different glm, trend and FDR methods.

There are a number of subtleties involved in the use of QL models.
The first is that the negative binomial dispersions \emph{must} be trended or common values.
This is because the function assumes that the supplied values are the true values.
For the trended/common values, the assumption is reasonable as information from many genes improves precision.
This is not the case for the tagwise dispersions due to the limited information for each gene.

Another subtlety involves the handling of zero counts.
Observations with fitted values of zero provide no residual degrees of freedom.
This must be considered when computing the value of the quasi-likelihood dispersion for genes with many zeros.
Finally, a lower bound is defined for the p-value of each gene, based on the likelihood ratio test.
This avoids spurious results involving weak shrinkage with very low quasi-likelihood dispersions.
}

\references{
Lund, SP, Nettleton, D, McCarthy, DJ, and Smyth, GK (2012).
Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates.
\emph{Statistical Applications in Genetics and Molecular Biology} Volume 11, Issue 5, Article 8.
\url{http://www.statsci.org/smyth/pubs/QuasiSeqPreprint.pdf}
}

\author{Davis McCarthy and Gordon Smyth, with modifications by Aaron Lun}

\examples{
nlibs <- 4
ntags <- 1000
dispersion.true <- 1/rchisq(ntags, df=10)
design <- model.matrix(~factor(c(1,1,2,2)))

# Generate count data
y <- rnbinom(ntags*nlibs,mu=20,size=1/dispersion.true)
y <- matrix(y,ntags,nlibs)
d <- DGEList(y)
d <- calcNormFactors(d)

# Fit the NB GLMs with QL methods
d <- estimateDisp(d, design)
fit <- glmQLFit(d, design)
results <- glmQLFTest(fit)
topTags(results)
fit <- glmQLFit(d, design, robust=TRUE)
results <- glmQLFTest(fit)
topTags(results)
fit <- glmQLFit(d, design, abundance.trend=FALSE)
results <- glmQLFTest(fit)
topTags(results)
}

\seealso{
\code{\link{topTags}} displays results from \code{glmQLFTest}.

\code{\link{plotQLDisp}} can be used to visualize the distribution of QL dispersions after EB shrinkage from \code{glmQLFit}.

The \code{QuasiSeq} package gives an alternative implementation of \code{glmQLFTest} based on the same statistical ideas.
}

\keyword{models}
