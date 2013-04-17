\name{camera.DGEList}
\alias{camera.DGEList}
\title{Competitive Gene Set Test for Digital Gene Expression Data Accounting for Inter-gene Correlation}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, accounting for inter-gene correlation.
}
\usage{
\method{camera}{DGEList}(y, index, design, contrast=ncol(design), weights=NULL, use.ranks=FALSE, allow.neg.cor=TRUE, trend.var=FALSE, sort=TRUE)
}

\arguments{
  \item{y}{\code{DGEList} object.}
  \item{index}{an index vector or a list of index vectors.  Can be any vector such that \code{y[indices,]} selects the rows corresponding to the test set.}
  \item{design}{design matrix.}
  \item{contrast}{contrast of the linear model coefficients for which the test is required. Can be an integer specifying a column of \code{design}, or else a numeric vector of same length as the number of columns of \code{design}.}
  \item{weights}{can be a numeric matrix of individual weights, of same size as \code{y}, or a numeric vector of array weights with length equal to \code{ncol(y)}.}
  \item{use.ranks}{do a rank-based test (\code{TRUE}) or a parametric test (\code{FALSE}?}
  \item{allow.neg.cor}{should reduced variance inflation factors be allowed for negative correlations?}
  \item{trend.var}{logical, should an empirical Bayes trend be estimated?  See \code{\link{eBayes}} for details.}
  \item{sort}{logical, should the results be sorted by p-value?}
}

\details{
This function implements a method proposed by Wu and Smyth (2012) for the digital gene expression data, eg. RNA-Seq data.
\code{camera} performs a \emph{competitive} test in the sense defined by Goeman and Buhlmann (2007).
It tests whether the genes in the set are highly ranked in terms of differential expression relative to genes not in the set.
It has similar aims to \code{geneSetTest} but accounts for inter-gene correlation.
See \code{\link{roast.DGEList}} for an analogous \emph{self-contained} gene set test.

The function can be used for any sequencing experiment which can be represented by a Negative Binomial generalized linear model.
The design matrix for the experiment is specified as for the \code{\link{glmFit}} function, and the contrast of interest is specified as for the \code{\link{glmLRT}} function.
This allows users to focus on differential expression for any coefficient or contrast in a model by giving the vector of test statistic values.

\code{camera} estimates p-values after adjusting the variance of test statistics by an estimated variance inflation factor.
The inflation factor depends on estimated genewise correlation and the number of genes in the gene set.
}

\value{
A data.frame. See \code{\link{camera}} for details.
}

\author{Yunshun Chen, Gordon Smyth}

\references{
Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\url{http://nar.oxfordjournals.org/content/40/17/e133}

Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 
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
\keyword{htest}

