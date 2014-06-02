\name{mglm}
\alias{mglm}
\alias{mglmOneGroup}
\alias{mglmOneWay}
\alias{mglmLevenberg}
\alias{designAsFactor}

\title{Fit Negative Binomial Generalized Linear Model to Multiple Response Vectors: Low 
Level Functions}

\description{
Fit the same log-link negative binomial or Poisson generalized linear model (GLM) to each 
row of a matrix of counts.
}

\usage{
mglmOneGroup(y, dispersion=0, offset=0, weights=NULL, maxit=50, tol=1e-10, verbose=FALSE)
mglmOneWay(y, design=NULL, dispersion=0, offset=0, weights=NULL, maxit=50, tol=1e-10)
mglmLevenberg(y, design, dispersion=0, offset=0, weights=NULL,
              coef.start=NULL, start.method="null", maxit=200, tol=1e-06)
designAsFactor(design)
}

\arguments{
\item{y}{numeric matrix containing the negative binomial counts.  Rows for tags and 
columns for libraries.}

\item{design}{numeric matrix giving the design matrix of the GLM.
Assumed to be full column rank.}

\item{dispersion}{numeric scalar or vector giving the dispersion parameter for each GLM.
Can be a scalar giving one value for all tags, or a vector of length equal to the number 
of tags giving tag-wise dispersions.}

\item{offset}{numeric vector or matrix giving the offset that is to be included in the 
log-linear model predictor.  Can be a scalar, a vector of length equal to the number of 
libraries, or a matrix of the same size as \code{y}.}

\item{weights}{numeric vector or matrix of non-negative quantitative weights.
Can be a vector of length equal to the number of libraries, or a matrix of the same size 
as \code{y}.}
\item{coef.start}{numeric matrix of starting values for the linear model coefficients.
Number of rows should agree with \code{y} and number of columns should agree with 
\code{design}.}

\item{start.method}{method used to generate starting values when \code{coef.stat=NULL}. 
Possible values are \code{"null"} to start from the null model of equal expression levels 
or \code{"y"} to use the data as starting value for the mean.}

\item{tol}{numeric scalar giving the convergence tolerance. For \code{mglmOneGroup}, 
convergence is judged successful when the step size falls below \code{tol} in absolute 
size.}

\item{maxit}{scalar giving the maximum number of iterations for the Fisher scoring 
algorithm.}

\item{verbose}{logical. If \code{TRUE}, warnings will be issued when \code{maxit} 
iterations are exceeded before convergence is achieved.}
}

\details{
The functions \code{mglmOneGroup}, \code{mglmOneWay} and \code{mglmLevenberg} all fit 
negative binomial generalized linear models, with the same design matrix but possibly 
different dispersions, offsets and weights, to a series of response vectors.
The functions are all low-level functions in that they operate on atomic objects such as 
matrices.
They are used as work-horses by higher-level functions in the edgeR package, especially 
by \code{glmFit}.

\code{mglmOneGroup} fit the null model, with intercept term only, to each response vector.
In other words, it treats the libraries as belonging to one group.
It implements Fisher scoring with a score-statistic stopping criterion for each tag.
Excellent starting values are available for the null model, so this function seldom has 
any problems with convergence.
It is used by other edgeR functions to compute the overall abundance for each tag.

\code{mglmLevenberg} fits an arbitrary log-linear model to each response vector.
It implements a Levenberg-Marquardt modification of the glm scoring algorithm to prevent 
divergence.
The main computation is implemented in C++.

All these functions treat the dispersion parameter of the negative binomial distribution 
as a known input.

\code{deviances.function} chooses the appropriate deviance function to use given a scalar 
or vector of dispersion parameters.
If the dispersion values are zero, then the Poisson deviance function is returned; if the 
dispersion values are positive, then the negative binomial deviance function is returned.
}

\value{
\code{mglmOneGroup} produces a vector of length equal to the number of tags/genes (number 
of rows of \code{y}) providing the single coefficent from the GLM fit for each tag/gene. 
This can be interpreted as a measure of the 'average expression' level of the tag/gene.

\code{mglmLevenberg} produces a list with the following components: 
	\item{coefficients}{matrix of estimated coefficients for the linear models}
	\item{fitted.values}{matrix of fitted values}
	\item{deviance}{residual deviances}
	\item{iter}{number of iterations used}
	\item{fail}{logical vector indicating tags for which the maximum damping was 
exceeded before convergence was achieved}

\code{deviances.function} returns a function to calculate the deviance as appropriate for 
the given values of the dispersion.

\code{designAsFactor} returns a factor of length equal to \code{nrow(design)}.
}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor 
RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\author{Gordon Smyth, Yunshun Chen, Davis McCarthy, Aaron Lun.  C++ code by Aaron Lun.}

\examples{
y <- matrix(rnbinom(1000,mu=10,size=2),ncol=4)
lib.size <- colSums(y)
dispersion <- 0.1

abundance <- mglmOneGroup(y, dispersion=dispersion, offset=log(lib.size))
AveLogCPM <- log1p(exp(1e6*abundance))/log(2)
summary(AveLogCPM)

## Same as above:
AveLogCPM <- aveLogCPM(y, dispersion, offset=log(lib.size))

## Fit the NB GLM to the counts with a given design matrix
f1 <- factor(c(1,1,2,2))
f2 <- factor(c(1,2,1,2))
x <- model.matrix(~f1+f2)
fit <- mglmLevenberg(y, x, dispersion=dispersion, offset=log(lib.size))
head(fit$coefficients)
}

\seealso{
\code{\link{glmFit}}, for more object-orientated GLM modelling for DGE data.
}
