\name{roast.DGEList}
\alias{roast.DGEList}
\alias{mroast.DGEList}
\alias{fry.DGEList}
\title{Rotation Gene Set Tests for Digital Gene Expression Data}
\description{
Rotation gene set testing for Negative Binomial generalized linear models.
}

\usage{
\method{roast}{DGEList}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
      set.statistic = "mean", gene.weights = NULL, \dots)

\method{mroast}{DGEList}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
       set.statistic = "mean", gene.weights = NULL,
       adjust.method = "BH", midp = TRUE, sort = "directional", ...)

\method{fry}{DGEList}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
      sort = "directional", ...)
}

\arguments{
  \item{y}{\code{DGEList} object.}
  \item{index}{index vector specifying which rows (probes) of \code{y} are in the test set.
        Can be a vector of integer indices, or a logical vector of length \code{nrow(y)}, or a vector of gene IDs corresponding to entries in \code{geneid}.
        Alternatively it can be a data.frame with the first column containing the index vector and the second column containing directional gene weights.
        For \code{mroast} or \code{fry}, \code{index} is a list of index vectors or a list of data.frames. }
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required.
        Can be an integer specifying a column of \code{design}, or the name of a column of \code{design}, or a numeric contrast vector of length equal to the number of columns of \code{design}.}
  \item{geneid}{gene identifiers corresponding to the rows of \code{y}.
        Can be either a vector of length \code{nrow(y)} or the name of the column of \code{y$genes} containing the gene identifiers.
        Defaults to \code{rownames(y)}.}
  \item{set.statistic}{summary set statistic. Possibilities are \code{"mean"},\code{"floormean"},\code{"mean50"} or \code{"msq"}.}
  \item{gene.weights}{numeric vector of directional (positive or negative) genewise weights.
        For \code{mroast} or \code{fry}, this vector must have length equal to \code{nrow(y)}.
        For \code{roast}, can be of length \code{nrow(y)} or of length equal to the number of genes in the test set.} 
  \item{adjust.method}{method used to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for possible values.}
  \item{midp}{logical, should mid-p-values be used in instead of ordinary p-values when adjusting for multiple testing?}
  \item{sort}{character, whether to sort output table by directional p-value (\code{"directional"}), non-directional p-value (\code{"mixed"}), or not at all (\code{"none"}).}
  \item{\dots}{other arguments are currently ignored.}
}

\value{
\code{roast} produces an object of class \code{\link[limma:roast]{Roast}}. See \code{\link{roast}} for details.

\code{mroast} and \code{fry} produce a data.frame. See \code{\link{mroast}} for details.
}

\details{
The roast gene set test was proposed by Wu et al (2010) for microarray data.
This function makes the roast test available for digital gene expression data.
The negative binomial count data is converted to approximate normal deviates by computing mid-p quantile residuals (Dunn and Smyth, 1996; Routledge, 1994) under the null hypothesis that the contrast is zero.
See \code{\link{roast}} for more description of the test and for a complete list of possible arguments.

The design matrix defaults to the \code{model.matrix(~y$samples$group)}.

\code{mroast} performs \code{roast} tests for a multiple of gene sets.
}

\seealso{
\code{\link{roast}}, \code{\link{camera.DGEList}}
}

\author{Yunshun Chen and Gordon Smyth}

\references{
Dunn, PK, and Smyth, GK (1996).
Randomized quantile residuals.
\emph{J. Comput. Graph. Statist.}, 5, 236-244. 
\url{http://www.statsci.org/smyth/pubs/residual.html}

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
