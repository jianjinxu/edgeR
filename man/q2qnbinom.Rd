\name{q2qnbinom}
\alias{q2qpois}
\alias{q2qnbinom}
\title{Quantile to Quantile Mapping between Negative-Binomial Distributions}

\description{Approximate quantile to quantile mapping between negative-binomial distributions with the same dispersion but different means.
The Poisson distribution is a special case.}

\usage{
q2qpois(x, input.mean, output.mean)
q2qnbinom(x, input.mean, output.mean, dispersion=0)
}

\arguments{ 

\item{x}{numeric matrix of unadjusted count data from a \code{DGEList} object}

\item{input.mean}{numeric matrix of estimated mean counts for tags/genes in unadjusted libraries }

\item{output.mean}{numeric matrix of estimated mean counts for tags/genes in adjusted (equalized) libraries, the same for all tags/genes in a particular group, different between groups}

\item{dispersion}{numeric scalar, vector or matrix of \code{dispersion} parameters}
}

\details{
This function finds the quantile with the same left and right tail probabilities relative to the output mean as \code{x} has relative to the input mean.
\code{q2qpois} is equivalent to \code{q2qnbinom} with \code{dispersion=0}.

This is the function that generates the pseudodata for \code{equalizeLibSizes}, and is required by \code{estimateCommonDisp} to adjust the library sizes to be equal.
The function takes values of the estimated mean for the unadjusted libraries (input.mean) and the estimated mean for the equalized libraries (output.mean) for each tag, as well as values (tagwise or common) for the dispersion parameter (\code{dispersion}).

The function calculates the percentiles that the counts in the unadjusted library represent for the normal and gamma distributions with mean and variance defined by the negative binomial rules: mean=\code{input.mean} and variance=\code{input.mean*(1+dispersion*input.mean)}.
The percentiles are then used to obtain quantiles from the normal and gamma distributions respectively, with mean and variance now defined as above but using \code{output.mean} instead of \code{input.mean}.
The function then returns as the pseudodata, i.e., equalized libraries, the arithmetic mean of the quantiles for the normal and the gamma distributions.

The mapping performed by this function is approximate but reasonably accurate.
Unlike the \code{qnbinom} function, it does not produce \code{Inf} when the right tail probabilities become small.
}

\value{numeric matrix of the same size as \code{x} with quantile-adjusted pseudodata }

\author{Gordon Smyth}

\examples{
y<-matrix(rnbinom(10000,size=2,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000,1010),2))
conc<-estimatePs(d,r=2)
N<-exp(mean(log(d$samples$lib.size)))
in.mean<-matrix(0,nrow=nrow(d$counts),ncol=ncol(d$counts))
out.mean<-matrix(0,nrow=nrow(d$counts),ncol=ncol(d$counts))
for(i in 1:2) {
	in.mean[,d$samples$group==i]<-outer(conc$conc.group[,i],d$samples$lib.size[d$samples$group==i])
	out.mean[,d$samples$group==i]<-outer(conc$conc.group[,i],rep(N,sum(d$samples$group==i)))
}
pseudo<-q2qnbinom(d$counts, input.mean=in.mean, output.mean=out.mean, dispersion=0.5)
}
\keyword{algebra}
