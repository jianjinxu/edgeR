\name{nbinomDeviance}
\alias{nbinomDeviance}
\alias{nbinomUnitDeviance}

\title{Negative Binomial Deviance}

\description{
Fit the same log-link negative binomial or Poisson generalized linear model (GLM) to each 
row of a matrix of counts.
}

\usage{
nbinomUnitDeviance(y, mean, dispersion=0)
nbinomDeviance(y, mean, dispersion=0, weights=NULL)
}

\arguments{
\item{y}{numeric vector or matrix containing the negative binomial counts.  If a matrix, 
then rows for tags and columns for libraries.  \code{nbinomDeviance} treats a vector as a 
matrix with one row.}

\item{mean}{numeric vector matrix of expected values, of same dimension as \code{y}.}

\item{dispersion}{numeric vector or matrix of negative binomial dispersions.
Can be a scalar, or a vector of length equal to the number of tags, or a matrix of same 
dimensions as \code{y}.}

\item{weights}{numeric vector or matrix of non-negative weights, as for \code{glmFit}.}
}

\details{
\code{nbinomUnitDeviance} computes the unit deviance for each \code{y} observation.
\code{nbinomDeviance} computes the total residual deviance for each row of \code{y} 
observation, i.e., weighted row sums of the unit deviances.

Care is taken to ensure accurate computation for small dispersion values.
}

\value{
\code{nbinomUnitDeviance} returns a numeric vector or matrix of the same size as \code{y}.

\code{nbinomDeviance} returns a numeric vector of length equal to the number of rows of 
\code{y}.
}

\references{
Jorgensen, B. (2006).
Generalized linear models. Encyclopedia of Environmetrics, Wiley.
\url{http://onlinelibrary.wiley.com/doi/10.1002/9780470057339.vag010/full}.

McCarthy, DJ, Chen, Y, Smyth, GK (2012).
Differential expression analysis of multifactor RNA-Seq experiments with respect to 
biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\author{Gordon Smyth, Yunshun Chen, Aaron Lun.  C++ code by Aaron Lun.}

\examples{
y <- matrix(1:6,3,2)
mu <- matrix(3,3,2)
nbinomUnitDeviance(y,mu,dispersion=0.2)
nbinomDeviance(y,mu,dispersion=0.2)
}
