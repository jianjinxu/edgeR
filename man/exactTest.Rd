\name{exactTest}
\alias{exactTest}
\alias{exactTestDoubleTail}
\alias{exactTestBetaApprox}
\alias{exactTestBySmallP}
\alias{exactTestByDeviance}

\title{Exact Tests for Differences between Two Groups of Negative-Binomial Counts}

\description{Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts.}

\usage{
exactTest(object, pair=NULL, dispersion="auto", rejection.region="doubletail", big.count=900)
exactTestDoubleTail(y1, y2, dispersion=0, big.count=900)
exactTestBySmallP(y1, y2, dispersion=0, big.count=900)
exactTestByDeviance(y1, y2, dispersion=0, big.count=900)
exactTestBetaApprox(y1, y2, dispersion=0)
}

\arguments{
\item{object}{an object of class \code{\link[edgeR:DGEList-class]{DGEList}}.}

\item{pair}{vector of length two, either numeric or character, providing the pair of groups to be compared; if a character vector, then should be the names of two groups (e.g. two levels of \code{object$samples$group}); if numeric, then groups to be compared are chosen by finding the levels of \code{object$samples$group} corresponding to those numeric values and using those levels as the groups to be compared; if \code{NULL}, then first two levels of \code{object$samples$group} (a factor) are used. Note that the first group listed in the pair is the baseline for the comparison---so if the pair is \code{c("A","B")} then the comparison is \code{B - A}, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change).}

\item{dispersion}{
either a numeric vector of dispersions or a character string indicating that dispersions should be taken from the data object.
If a numeric vector, then can be either of length one or of length equal to the number of tags.
Allowable character values are \code{"common"}, \code{"tagwise"} or \code{"auto"}.
\code{"common"} uses \code{object$common.dispersion}, \code{"tagwise"} uses \code{object$tagwise.dispersion} and \code{"auto"} uses tagwise when available and otherwise common.
}

\item{rejection.region}{type of rejection region for two-sided exact test.  Possible values are \code{"doubletail"}, \code{"smallp"} or \code{"deviance"}.}

\item{big.count}{count size above which asymptotic beta approximation will be used.}

\item{y1}{numeric matrix of counts for the first the two experimental groups to be tested for differences.
Rows correspond to genes or transcripts and columns to libraries.
Libraries are assumed to be equal in size - e.g. adjusted pseudocounts from the output of \code{\link{equalizeLibSizes}}.}

\item{y2}{numeric matrix of counts for the second of the two experimental groups to be tested for differences.
Rows correspond to genes or transcripts and columns to libraries.
Libraries are assumed to be equal in size - e.g. adjusted pseudocounts from the output of \code{\link{equalizeLibSizes}}. Must have the same number of rows as \code{y1}.}
}

\value{
\code{exactTestDoubleTail} and friends produce a numeric vector of genewise p-values, one for each row of \code{y1} and \code{y2}.

\code{exactTest} produces an object of class \code{DGEExact} containing the following components:
\item{table}{a data frame containing the elements \code{logConc}, the log-average concentration/abundance for each tag in the two groups being compared, \code{logFC}, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, \code{p.value}, exact p-value for differential expression using the NB model}
\item{comparison}{a vector giving the names of the two groups being compared} 
\item{genes}{a data frame containing information about each transcript; taken from \code{object} and can be \code{NULL}}
}

\details{
The functions test for differential expression between two groups of count libraries.
They implement the exact test proposed by Robinson and Smyth (2008) for a difference in mean between two groups of negative binomial random variables.
The functions accept two groups of count libraries, and a test is performed for each row of data.
For each row, the test is conditional on the sum of counts for that row.
The test can be viewed as a generalization of the well-known exact binomial test, implemented in the function \code{binom.test} in the stats package, but generalized to overdispersed counts.

The low level functions \code{exactTestDoubleTail}, \code{exactTestBetaApprox}, \code{exactTestBySmallP} and \code{exactTestByDeviance} all assume
that the libraries have been normalized to have the same size (expected column sum under the null hypothesis).
The higher level function \code{exactTest} is intended to be called by users.
This has a more object-orientated flavor and produces an object containing all the necessary components for downstream analysis.
\code{exactTest} equalizes the library sizes using \code{\link{equalizeLibSizes}} before calling one of the low level functions.

The functions \code{exactTestDoubleTail}, \code{exactTestBySmallP} and \code{exactTestByDeviance} correspond to different ways to define the two-sided rejection region when the two groups have different numbers of samples.
\code{exactTestBySmallP} implements the method of small probabilities as proposed by Robinson and Smyth (2008).
This method corresponds to \code{binom.test} when the dispersion is near zero, but gives poor results when the dispersion is very large.
\code{exactTestDoubleTail} computes two-sided p-values by doubling the smaller tail probability.
\code{exactTestByDeviance} uses the deviance goodness of fit statistics to define the rejection region, and is therefore equivalent to a conditional likelihood ratio test.
This has good statistical properties but is relatively slow to compute.
For general remarks on different types of rejection regions for exact tests see Gibbons and Pratt (1975).

\code{exactTestBetaApprox} implements an asymptotic beta distribution approximation to the conditional count distribution.
}

\references{
Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data.
\emph{Biostatistics}, 9, 321-332.

Gibbons, JD and Pratt, JW (1975).  P-values: interpretation and methodology.
\emph{The American Statistician} 29, 20-25.
}

\author{Mark Robinson, Davis McCarthy, Gordon Smyth}

\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(80,size=1/0.2,mu=10),nrow=20,ncol=4)
rownames(y) <- paste("Gene",1:nrow(y),sep=".")
group <- factor(c(1,1,2,2))
d <- DGEList(counts=y,group=group,lib.size=rep(1000,4))

# estimate dispersions and find differences in expression
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
de <- exactTest(d)
topTags(de)

# same example using low level exactTest function directly
p.value <- exactTestDoubleTail(y[,1:2],y[,3:4],dispersion=0.2)
}

\seealso{
\code{\link{equalizeLibSizes}}, \code{\link{binomTest}}
}

\keyword{algebra}

