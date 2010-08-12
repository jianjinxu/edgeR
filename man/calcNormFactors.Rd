\name{calcNormFactors}
\Rdversion{1.1}
\alias{calcNormFactors}
\title{Calculates Normalization Factors for a Matrix of Count Data}
\description{
Using a reference sample, calculate the normalization factors, over and above accounting for library size. 
}
\usage{
calcNormFactors(object, refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{either a \code{matrix} of raw (read) counts or a \code{DGEList} object}
  \item{refColumn}{column to use as reference}
  \item{logratioTrim}{amount of trim to use on log-ratios ("M" values)}
  \item{sumTrim}{amount of trim to use on the combined absolute levels ("A" values)}
  \item{doWeighting}{logical, whether to compute (asymptotic binomial precision) weights}
  \item{Acutoff}{cutoff on "A" values to use before trimming}
}
\details{
The weighted trimmed mean of M values (to the reference) is used as the normalization factor, where the weights are from the delta method on Binomial data (more details to come).

The normalization factor for the reference sample will always be 1.
}
\value{
If a \code{matrix} is given for \code{object}, the output is a vector with length \code{ncol(object)} giving the relative normalization factors.
If a \code{DGEList} object is given for \code{object}, the output is a \code{DGEList} object containing the normalization factors in the \code{samples$norm.factors} element.
}
\author{ Mark Robinson }
\examples{
  d <- matrix( rpois(1000, lambda=5), nrow=200 )
  f <- calcNormFactors(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
