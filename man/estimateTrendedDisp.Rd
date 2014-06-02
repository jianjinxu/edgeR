\name{estimateTrendedDisp}
\alias{estimateTrendedDisp}

\title{Estimate Empirical Bayes Trended Dispersion Values}

\description{
Estimates trended dispersion values by an empirical Bayes method.
}

\usage{
estimateTrendedDisp(object, method="bin.spline", df=5, span=2/3)
}

\arguments{ 

\item{object}{object of class \code{DGEList} containing (at least) the elements 
\code{counts} (table of raw counts), \code{group} (factor indicating group), 
\code{lib.size} (numeric vector of library sizes) and \code{pseudo.alt} (numeric matrix 
of quantile-adjusted pseudocounts calculated under the alternative hypothesis of a true 
difference between groups; recommended to use the \code{DGEList} object provided as the 
output of \code{estimateCommonDisp}}

\item{method}{method used to estimated the trended dispersions. Possible values are 
\code{"spline"}, and \code{"loess"}.}

\item{df}{integer giving the degrees of freedom of the spline function if \code{"spline"} 
method is used, see \code{ns} in the splines package. Default is 5.}

\item{span}{scalar, passed to \code{loess} to determine the amount of smoothing for the 
loess fit when \code{"loess"} method is used. Default is \code{2/3}.}
}

\details{
This function takes the binned common dispersion and abundance, and fits a smooth curve 
through these binned values using either natural cubic splines or loess. From this smooth 
curve it predicts the dispersion value for each gene based on the gene's overall 
abundance. This results in estimates for the NB dispersion parameter which have a 
dependence on the overall expression level of the gene, and thus have an 
abundance-dependent trend.
}

\value{
An object of class \code{DGEList} with the same components as for 
\code{\link{estimateCommonDisp}} plus the trended dispersion estimates for each gene or 
tag.
}

\author{Yunshun Chen and Gordon Smyth}
\examples{
ngenes <- 1000
nlib <- 4
log2cpm <- seq(from=0,to=16,length=ngenes)
lib.size <- 1e7
mu <- 2^log2cpm * lib.size * 1e-6
dispersion <- 1/sqrt(mu) + 0.1
counts <- rnbinom(ngenes*nlib, mu=mu, size=1/dispersion)
counts <- matrix(counts,ngenes,nlib)
y <- DGEList(counts,lib.size=rep(lib.size,nlib))
y <- estimateCommonDisp(y)
y <- estimateTrendedDisp(y)
y <- estimateTagwiseDisp(y)
plotBCV(y)
}

\seealso{
\code{\link{estimateCommonDisp}} estimates a common value for the dispersion parameter 
for all tags/genes - should generally be run before \code{estimateTrendedDisp}.
}

