\name{adjustedProfileLik}
\alias{adjustedProfileLik}

\title{Adjusted Profile Likelihood for the Negative Binomial Dispersion Parameter}

\description{
Compute adjusted profile-likelihoods for estimating the dispersion parameters of genewise negative binomial glms.
}

\usage{
adjustedProfileLik(dispersion, y, design, offset, adjust=TRUE)
}

\arguments{
\item{dispersion}{numeric scalar or vector of dispersions.}
\item{y}{numeric matrix of counts.}
\item{design}{numeric matrix giving the design matrix.}
\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalor or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix.}
\item{adjust}{logical, if \code{TRUE} then Cox-Reid adjustment is made to the log-likelihood, if \code{FALSE} then the log-likelihood is returned without adjustment.}
}

\value{
vector of adjusted profile log-likelihood values, one for each row of \code{y}.
}

\details{
For each row of data, compute the adjusted profile-likelihood for estimating the dispersion parameter of the negative binomial glm.
The adjusted profile likelihood is described by McCarthy et al (2012), and is based on the method of Cox and Reid (1987).

The adjusted profile likelihood is an approximate log-likelihood for the dispersion parameter, conditional on the estimated values of the coefficients in the NB log-linear models.
The conditional likelihood approach is a technique for adjusting the likelihood function to allow for the fact that nuisance parameters have to be estimated in order to evaluate the likelihood.
When estimating the dispersion, the nuisance parameters are the coefficients in the linear model.

This implementation calls the LAPACK library to perform the Cholesky decomposition during adjustment estimation.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\author{Yunshun Chen, Gordon Smyth, Aaron Lun}
\examples{
y <- matrix(rnbinom(1000, mu=10, size=2), ncol=4)
design <- matrix(1, 4, 1)
dispersion <- 0.5
apl <- adjustedProfileLik(dispersion, y, design, offset=0)
apl
}

\seealso{
\code{\link{glmFit}}
}
