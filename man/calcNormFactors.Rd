\name{calcNormFactors}
\alias{calcNormFactors}
\title{Calculate Normalization Factors to Align Columns of a Count Matrix}
\description{
Calculate normalization factors to scale the raw library sizes.
}
\usage{
calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"), refColumn = NULL,
      logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
}
\arguments{
  \item{object}{either a \code{matrix} of raw (read) counts or a \code{DGEList} object}
  \item{method}{normalization method to be used}
  \item{refColumn}{column to use as reference for \code{method="TMM"}. Can be a column number or a numeric vector of length \code{nrow(object))}.}
  \item{logratioTrim}{amount of trim to use on log-ratios ("M" values) for \code{method="TMM"}}
  \item{sumTrim}{amount of trim to use on the combined absolute levels ("A" values) for \code{method="TMM"}}
  \item{doWeighting}{logical, whether to compute (asymptotic binomial precision) weights for \code{method="TMM"}}
  \item{Acutoff}{cutoff on "A" values to use before trimming for \code{method="TMM"}}
  \item{p}{percentile (between 0 and 1) of the counts that is aligned when \code{method="upperquartile"}}
}

\details{
\code{method="TMM"} is the weighted trimmed mean of M-values (to the reference) proposed by Robinson and Oshlack (2010), where the weights are from the delta method on Binomial data.
If \code{refColumn} is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.

\code{method="RLE"} is the scaling factor method proposed by Anders and Huber (2010).
We call it "relative log expression", as median library is calculated from the geometric mean of all columns and the median ratio of each sample to the median library is taken as the scale factor.

\code{method="upperquartile"} is the upper-quartile normalization method of Bullard et al (2010), in which the scale factors are calculated from the 75\% quantile of the counts for each library,
after removing transcripts which are zero in all libraries.
This idea is generalized here to allow scaling by any quantile of the distributions.

If \code{method="none"}, then the normalization factors are set to 1.

For symmetry, normalization factors are adjusted to multiply to 1.
The effective library size is then the original library size multiplied by the scaling factor.

Note that rows that have zero counts for all columns are trimmed before normalization factors are computed.
Therefore rows with all zero counts do not affect the estimated factors.
}

\value{
If \code{object} is a \code{matrix}, the output is a vector with length \code{ncol(object)} giving the relative normalization factors.
If \code{object} is a \code{DGEList}, then it is returned as output with the relative normalization factors in \code{object$samples$norm.factors}.
}

\author{Mark Robinson, Gordon Smyth}

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
y <- matrix( rpois(1000, lambda=5), nrow=200 )
calcNormFactors(y)
}
