\name{glmFit}
\alias{glmFit}
\alias{glmFit.DGEList}
\alias{glmFit.default}
\alias{glmLRT}
\alias{glmQLFTest}

\title{Genewise Negative Binomial Generalized Linear Mdels}

\description{Fit a negative binomial generalized log-linear model to the read counts for each gene or transcript.
Conduct statistical tests based on the fitted models.}

\usage{
\method{glmFit}{DGEList}(y, design=NULL, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL,
       prior.count.total=0.5, start=NULL, method="auto", ...)
glmLRT(glmfit, coef=ncol(glmfit$design), contrast=NULL)
glmQLFTest(glmfit, coef=ncol(glmfit$design), contrast=NULL, abundance.trend=TRUE)
}

\arguments{
\item{y}{an object that contains the raw counts for each library (the measure of expression level); alternatively, a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the tagwise linear models.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{dispersion}{numeric scalar or vector providing the value for the dispersion parameter that is used in fitting the GLM for each transcript. Can be a common value for all tags, or a vector of values can provide a unique dispersion value for each tag. If \code{NULL} will be extracted from \code{y}, with order of precedence: tagwise dispersion, trended dispersions, common dispersion.}

\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalor or a vector of length \code{ncol{y}}, in which case it is expanded out to a matrix.}

\item{weights}{optional numeric matrix giving prior weights for the observations (for each library and transcript) to be used in the GLM calculations.  Not supported by methods \code{"linesearch"} or \code{"levenberg"}.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size}). Defaults to \code{colSums(y)}.}

\item{prior.count.total}{the total number of prior counts to be added to each row of data to shrink the estimated log-fold-changes towards zero.}

\item{start}{optional numeric matrix of initial estimates for the linear model coefficients.}

\item{method}{which fitting algorithm to use.  Possible values are \code{"auto"}, \code{"linesearch"}, \code{"levenberg"} or \code{"simple"}.}

\item{...}{other arguments are passed to lower-level functions, for example to \code{mglmLS}.}

\item{glmfit}{a \code{DGEGLM} object, usually as output from \code{glmFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.  Values must be columns or column names of \code{design}. Defaults to the last coefficient.  Ignored if \code{contrast} is specified.}

\item{contrast}{numeric vector specifying a contrast of the linear model coefficients to be tested equal to zero. Must be of length equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}

\item{abundance.trend}{logical, whether to allow an abundance-dependent trend when estimating the prior values for the quasi-likelihood multiplicative dispersion parameter.}
}

\value{
\code{glmFit} produces an object of class \code{DGEGLM} containing components \code{counts}, \code{samples}, \code{genes} and \code{abundance} from \code{y} plus the following new components:
	\item{design}{design matrix as input.}
	\item{weights}{matrix of weights as input.}
	\item{df.residual}{numeric vector of residual degrees of freedom, one for each tag.}
	\item{offset}{numeric matrix of linear model offsets.}
	\item{dispersion}{vector of dispersions used for the fit.}
	\item{coefficients}{numeric matrix of estimated coefficients from the glm fits, on the natural log scale, of size \code{nrow(y)} by \code{ncol(design)}.}
	\item{fitted.values}{matrix of fitted values from glm fits, same number of rows and columns as \code{y}.}
	\item{deviance}{numeric vector of deviances, one for each tag.}

\code{glmLRT} and \code{glmQFTest} produce objects of class \code{DGELRT} contain the same components as \code{glmfit} plus the following:
	\item{table}{data frame with same rows as \code{y} and columns containing the following:
	log2-fold change of expression between conditions being tested (\code{logFC}),
	average log2-counts per million (\code{logCPM}),
	likelihood ratio statistics (\code{LR}) for \code{glmLRT} or F-statistics (\code{F}) for \code{glmQFTest},
	p-values (\code{PValue}).}
	\item{comparison}{character vector of names of for the coefficients or the contrast tested.}
}

\details{
\code{glmFit} and \code{glmLRT} implement generalized linear model (glm) methods developed by McCarthy et al (2012).

\code{glmFit} fits genewise negative binomial glms, all with the same design matrix but possibly different dispersions, offsets and weights.
When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the glms are fitting very quickly using \code{\link{mglmOneGroup}}.
Otherwise the default fitting method, implemented in \code{\link{mglmLS}}, is a parallelized line search algorithm described by McCarthy et al (2012).
Other possible fitting methods are \code{\link{mglmLevenberg}} and \code{\link{mglmSimple}}.

Positive \code{prior.count.total} cause the returned coefficients to be shrunk in such a way that fold-changes between the treatment conditions are decreased.
In particular, infinite fold-changes are avoided.
Larger values cause more shrinkage.
The returned coefficients are affected but not the likelihood ratio tests or p-values.

\code{glmLRT} conducts likelihood ratio tests for one or more coefficients in the linear model.
If \code{coef} is used, the null hypothesis is that all the coefficients indicated by \code{coef} are equal to zero.
If \code{contrast} is non-null, then the null hypothesis is that the specified contrast of the coefficients is equal to zero.
For example, a contrast of \code{c(0,1,-1)}, assuming there are three coefficients, would test the hypothesis that the second and third coefficients are equal.

\code{glmQLFTest} implements the quasi-likelihood method of Lund et al (2012).
It behaves the same as \code{glmLRT} except that it replaces likelihood ratio tests with
quasi-likelihood F-tests for coefficients in the linear model.
This function calls the limma function \code{\link{squeezeVar}} to conduct empirical Bayes smoothing of the genewise multiplicative dispersions.
Note that the \code{QuasiSeq} package provides a alternative implementation of Lund et al (2012),
with slightly different glm, trend and FDR methods.
}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}

Lund, SP, Nettleton, D, McCarthy, DJ, Smyth, GK (2012). Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates.
\emph{Statistical Applications in Genetics and Molecular Biology}. (Accepted 31 July 2012)
}

\author{Davis McCarthy and Gordon Smyth}

\examples{
nlibs <- 3
ntags <- 100
dispersion.true <- 0.1

# Make first transcript respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ntags-1)))
mu.true <- 2^(beta.true \%*\% t(design))

# Generate count data
y <- rnbinom(ntags*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ntags,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("Gene",1:ntags,sep="")
d <- DGEList(y)

# Normalize
d <- calcNormFactors(d)

# Fit the NB GLMs
fit <- glmFit(d, design, dispersion=dispersion.true)

# Likelihood ratio tests for trend
results <- glmLRT(fit, coef=2)
topTags(results)

# Estimate the dispersion (may be unreliable with so few tags)
d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
}

\seealso{
Low-level computations are done by \code{\link{mglmOneGroup}}, \code{\link{mglmLS}}, \code{\link{mglmLevenberg}} or \code{\link{mglmSimple}}.

\code{\link{topTags}} displays results from \code{glmLRT} or \code{glmQLFTest}.

The \code{QuasiSeq} package gives an alternative implementation of \code{glmQLFTest} based on the same statistical ideas.
}

\keyword{models}
