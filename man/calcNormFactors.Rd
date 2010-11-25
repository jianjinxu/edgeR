\name{calcNormFactors}
\Rdversion{1.1}
\alias{calcNormFactors}
\title{Calculates Normalization Factors for a Matrix of Count Data}
\description{
Using a reference sample, calculate the normalization factors, over and above accounting for library size. 
}
\usage{
calcNormFactors(object, method=c("TMM","RLE","quantile"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, quantile=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{either a \code{matrix} of raw (read) counts or a \code{DGEList} object}
  \item{method}{method to use to calculate the scale factors}
  \item{refColumn}{column to use as reference, only used when \code{method="TMM"}}
  \item{logratioTrim}{amount of trim to use on log-ratios ("M" values), only used when \code{method="TMM"}}
  \item{sumTrim}{amount of trim to use on the combined absolute levels ("A" values), only used when \code{method="TMM"}}
  \item{doWeighting}{logical, whether to compute (asymptotic binomial precision) weights, only used when \code{method="TMM"}}
  \item{Acutoff}{cutoff on "A" values to use before trimming, only used when \code{method="TMM"}}
  \item{quantile}{quantile used to compute scale factors from, only used when \code{method="Quantile"}}
}
\details{

When \code{method="TMM"}, the weighted trimmed mean of M values (to the reference) is used as the normalization factor, where the weights are from the delta method on Binomial data.  If \code{refColumn} is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.  When \code{method="RLE"} (which stands for relative log expression), a median library is calculated from the geometric mean of all columns and the median ratio of each sample to the median library is taken as the scale factor (this is the implementation proposed by the DESeq package).  When \code{method="Quantile"}, the scale factors are calculated from the quantiles (default=75%) of each library.

For symmetry, normalization factors are adjusted to multiply to 1.

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
