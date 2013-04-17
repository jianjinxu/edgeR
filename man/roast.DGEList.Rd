\name{roast.DGEList}
\alias{roast.DGEList}
\alias{mroast.DGEList}
\title{Rotation Gene Set Tests for Digital Gene Expression Data}
\description{
Rotation gene set testing for Negative Binomial generalized linear models.
}

\usage{
\method{roast}{DGEList}(y, index=NULL, design=NULL, contrast=ncol(design), set.statistic="mean",
     gene.weights=NULL, array.weights=NULL, weights=NULL, block=NULL, correlation,
     var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=999)
\method{mroast}{DGEList}(y, index=NULL, design=NULL, contrast=ncol(design), set.statistic="mean",
     gene.weights=NULL, array.weights=NULL, weights=NULL, block=NULL, correlation,
     var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=999, adjust.method="BH", midp=TRUE, sort="directional")
}

\arguments{
  \item{y}{\code{DGEList} object.}
  \item{index}{index vector specifying which rows (genes) of \code{y} are in the test set.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{y[iset,]} contains the values for the gene set to be tested.}
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{set.statistic}{summary set statistic. Possibilities are \code{"mean"},\code{"floormean"},\code{"mean50"} or \code{"msq"}.}
  \item{gene.weights}{optional numeric vector of weights for genes in the set. Can be positive or negative.  For \code{mroast.DGEList} this vector must have length equal to \code{nrow(y)}.  For \code{roast.DGEList}, can be of length \code{nrow(y)} or of length equal to the number of genes in the test set.} 
  \item{array.weights}{optional numeric vector of array weights.}
  \item{weights}{optional matrix of observation weights. If supplied, should be of same dimensions as \code{y} and all values should be positive.}
  \item{block}{optional vector of blocks.}
  \item{correlation}{correlation between blocks.}
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{trend.var}{logical, should a trend be estimated for \code{var.prior}?  See \code{eBayes} for details.  Only used if \code{var.prior} or \code{df.prior} are \code{NULL}.}
  \item{nrot}{number of rotations used to estimate the p-values.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for possible values.}
  \item{midp}{logical, should mid-p-values be used in instead of ordinary p-values when adjusting for multiple testing?}
  \item{sort}{character, whether to sort output table by directional p-values (\code{"directional"}), non-directional p-value (\code{"mixed"}), or not at all (\code{"none"}).}
}

\value{
\code{roast} produces an object of class \code{\link[limma:roast]{Roast}}. See \code{\link{roast}} for details.

\code{mroast} produces a data.frame. See \code{\link{mroast}} for details.
}

\details{
This function implements a method for the ROAST gene set test from Wu et al (2010) for the digital gene expression data, eg. RNA-Seq data.
Basically, the Negative Binomial generalized linear models are fitted for count data. The fitted values are converted into z-scores, and then it calls the \code{roast} function in \code{limma} package to conduct the gene set test.
It tests whether any of the genes in the set are differentially expressed.
This allows users to focus on differential expression for any coefficient or contrast in a generalized linear model.
If \code{contrast} is not specified, the last coefficient in the model will be tested.
The arguments \code{array.weights}, \code{block} and \code{correlation} have the same meaning as they for for the \code{\link{lmFit}} function.

The arguments \code{df.prior} and \code{var.prior} have the same meaning as in the output of the \code{\link{eBayes}} function.
If these arguments are not supplied, they are estimated exactly as is done by \code{eBayes}.

The argument \code{gene.weights} allows directions or weights to be set for individual genes in the set.

The gene set statistics \code{"mean"}, \code{"floormean"}, \code{"mean50"} and \code{msq} are defined by Wu et al (2010).
The different gene set statistics have different sensitivities to small number of genes.
If \code{set.statistic="mean"} then the set will be statistically significantly only when the majority of the genes are differentially expressed.
\code{"floormean"} and \code{"mean50"} will detect as few as 25\% differentially expressed.
\code{"msq"} is sensitive to even smaller proportions of differentially expressed genes, if the effects are reasonably large.

The output gives p-values three possible alternative hypotheses, 
\code{"Up"} to test whether the genes in the set tend to be up-regulated, with positive t-statistics,
\code{"Down"} to test whether the genes in the set tend to be down-regulated, with negative t-statistics,
and \code{"Mixed"} to test whether the genes in the set tend to be differentially expressed, without regard for direction.

\code{roast} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals (Langsrud, 2005), so p-values will vary slightly from run to run.
To get more precise p-values, increase the number of rotations \code{nrot}.
The p-value is computed as \code{(b+1)/(nrot+1)} where \code{b} is the number of rotations giving a more extreme statistic than that observed (Phipson and Smyth, 2010).
This means that the smallest possible p-value is \code{1/(nrot+1)}.

\code{mroast} does roast tests for multiple sets, including adjustment for multiple testing.
By default, \code{mroast} reports ordinary p-values but uses mid-p-values (Routledge, 1994) at the multiple testing stage.
Mid-p-values are probably a good choice when using false discovery rates (\code{adjust.method="BH"}) but not when controlling the family-wise type I error rate (\code{adjust.method="holm"}).

\code{roast} performs a \emph{self-contained} test in the sense defined by Goeman and Buhlmann (2007).
For a \emph{competitive} gene set test, see \code{\link{camera.DGEList}}.
}

\seealso{
\code{\link{roast}}, \code{\link{camera.DGEList}}
}

\author{Yunshun Chen and Gordon Smyth}

\references{
Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 

Langsrud, O (2005).
Rotation tests.
\emph{Statistics and Computing} 15, 53-60.

Phipson B, and Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Article 39.

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments.
\emph{Bioinformatics} 26, 2176-2182.
\url{http://bioinformatics.oxfordjournals.org/content/26/17/2176}
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

roast(y, iset1, design, contrast=2)
mroast(y, iset1, design, contrast=2)
mroast(y, list(set1=iset1, set2=iset2), design, contrast=2)
}
\keyword{htest}
