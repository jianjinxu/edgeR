\name{camera.DGEList}
\alias{camera.DGEList}
\title{Competitive Gene Set Test for Digital Gene Expression Data Accounting for 
Inter-gene Correlation}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of 
differential expression, accounting for inter-gene correlation.
}
\usage{
\method{camera}{DGEList}(y, index, design=NULL, contrast=ncol(design), \dots)
}

\arguments{
  \item{y}{a \code{DGEList} object containing dispersion estimates.}
  \item{index}{an index vector or a list of index vectors.  Can be any vector such that 
\code{y[indices,]} selects the rows corresponding to the test set.}
  \item{design}{the design matrix.}
  \item{contrast}{the contrast of the linear model coefficients for which the test is 
required. Can be an integer specifying a column of \code{design}, or else a numeric 
vector of same length as the number of columns of \code{design}.}
  \item{\dots}{other arguments are passed to \code{\link{camera.default}}.}
}

\details{
The camera gene set test was proposed by Wu and Smyth (2012) for microarray data.
This function makes the camera test available for digital gene expression data.
The negative binomial count data is converted to approximate normal deviates by computing 
mid-p quantile residuals (Dunn and Smyth, 1996; Routledge, 1994) under the null 
hypothesis that the contrast is zero.
See \code{\link{camera}} for more description of the test and for a complete list of 
possible arguments.

The design matrix defaults to the \code{model.matrix(~y$samples$group)}.
}

\value{
A data.frame. See \code{\link{camera}} for details.
}

\author{Yunshun Chen, Gordon Smyth}

\references{
Dunn, K. P., and Smyth, G. K. (1996).
Randomized quantile residuals.
\emph{J. Comput. Graph. Statist.}, 5, 236-244. 
\url{http://www.statsci.org/smyth/pubs/residual.html}

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.

Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for 
inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\url{http://nar.oxfordjournals.org/content/40/17/e133}
}

\seealso{
\code{\link{roast.DGEList}},
\code{\link{camera}}.
}

\examples{
mu <- matrix(10, 100, 4)
group <- factor(c(0,0,1,1))
design <- model.matrix(~group)

# First set of 10 genes that are genuinely differentially expressed
iset1 <- 1:10
mu[iset1,3:4] <- mu[iset1,3:4]+10

# Second set of 10 genes are not DE
iset2 <- 11:20

# Generate counts and create a DGEList object
y <- matrix(rnbinom(100*4, mu=mu, size=10),100,4)
y <- DGEList(counts=y, group=group)

# Estimate dispersions
y <- estimateDisp(y, design)

camera(y, iset1, design)
camera(y, iset2, design)

camera(y, list(set1=iset1,set2=iset2), design)
}
