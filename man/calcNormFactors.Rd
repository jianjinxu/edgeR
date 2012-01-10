\name{calcNormFactors}
\Rdversion{1.1}
\alias{calcNormFactors}
\title{Calculates Normalization Factors for a Matrix of Count Data}
\description{
Using a reference sample, calculate the normalization factors, over and above accounting for library size. 
}
\usage{
calcNormFactors(object, method=c("TMM","RLE","upperquartile"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
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
  \item{p}{percentile (between 0 and 1) used to compute scale factors from, only used when \code{method="upperquartile"}}
}
\details{
\code{method="TMM"} is the weighted trimmed mean of M-values (to the reference) proposed by Robinson and Oshlack (2010), where the weights are from the delta method on Binomial data.  If \code{refColumn} is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.

\code{method="RLE"} is the scaling factor method proposed by Anders and Huber (2010).
We call it "relative log expression", as median library is calculated from the geometric mean of all columns and the median ratio of each sample to the median library is taken as the scale factor.

\code{method="upperquartile"} is the upper-quartile normalization method of Bullard et al (2010), in which the scale factors are calculated from the 75\% quantile of the counts for each library,
after removing transcripts which are zero in all libraries.
We generalize it to allow scaling by any quantile of the distributions.

For symmetry, normalization factors are adjusted to multiply to 1.
}

\value{
If a \code{matrix} is given for \code{object}, the output is a vector with length \code{ncol(object)} giving the relative normalization factors.
If a \code{DGEList} object is given for \code{object}, the output is a \code{DGEList} object containing the normalization factors in the \code{samples$norm.factors} element.
}

\author{ Mark Robinson }

\references{
Anders, S, Huber, W (2010).
Differential expression analysis for sequence count data
\emph{Genome Biology} 11, R106.

Bullard JH, Purdom E, Hansen KD, Dudoit S. (2010)
Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments.
\emph{BMC Bioinformatics} 11, 94.
A scaling normalization method for differential expression analysis of RNA-seq data.

Robinson MD, Oshlack A (2010).
\emph{Genome Biology} 11, R25.
}

\examples{
  d <- matrix( rpois(1000, lambda=5), nrow=200 )
  f <- calcNormFactors(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
