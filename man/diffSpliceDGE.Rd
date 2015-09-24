\name{diffSpliceDGE}
\alias{diffSpliceDGE}

\title{Test for Differential Exon Usage}
\description{Given a negative binomial generalized log-linear model fit at the exon level, test for differential exon usage between experimental conditions.}
\usage{
diffSpliceDGE(fit.exon, coef=ncol(fit.exon$design), contrast=NULL, geneid, exonid=NULL, verbose=TRUE)
}

\arguments{
  \item{fit.exon}{an \code{DGEGLM} fitted model object produced by \code{glmFit}. Rows should correspond to exons.}
  \item{coef}{integer indicating which coefficient of the generalized linear model is to be tested for differential exon usage. Defaults to the last coefficient.}
  \item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested for differential exon usage. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(fit.exon)} or the name of the column of \code{fit.exon$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{exon identifiers. Either a vector of length \code{nrow(fit.exon)} or the name of the column of \code{fit.exon$genes} containing the exon identifiers.}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
}

\value{
\code{diffSpliceDGE} produces an object of class \code{DGELRT} containing the component \code{design} from \code{fit.exon} plus the following new components:
  \item{comparison}{character string describing the coefficient being tested.}
  \item{coefficients}{numeric vector of coefficients on the natural log scale. Each coefficient is the difference between the log-fold-change for that exon versus the average log-fold-change for the rest exons within the same gene.}
  \item{exon.F}{numeric vector of F-statistics for exons.}
  \item{exon.df.test}{numeric vector of testing degrees of freedom for exons.}
  \item{exon.df.prior}{numeric vector of prior degrees of freedom for exons.}
  \item{exon.df.residual}{numeric vector of residual degrees of freedom for exons.}
  \item{exon.p.value}{numeric vector of p-values for exons.}
  \item{genes}{data.frame of exon annotation}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs.}
  \item{exoncolname}{character string giving the name of the column of \code{genes} containing exon IDs.}
  \item{gene.df.test}{numeric vector of testing degrees of freedom for genes.}
  \item{gene.df.prior}{numeric vector of prior degrees of freedom for genes.}
  \item{gene.df.residual}{numeric vector of residual degrees of freedom for genes.}
  \item{gene.Simes.p.value}{numeric vector of Simes' p-values for genes.}
  \item{gene.F}{numeric vector of F-statistics for gene-level test.}
  \item{gene.F.p.value}{numeric vector of F-test p-values for genes.}
  \item{gene.genes}{data.frame of gene annotation.}
The information and testing results for both exons and genes are sorted by geneid and by exonid within gene.
}

\details{
This function tests for differential exon usage for each gene for a given coefficient of the generalized linear model.

Testing for differential exon usage is equivalent to testing whether the exons in each gene have the same log-fold changes as the other exons in the same gene. 
At exon-level, each exon is compared to the average of all other exons for the same gene using quasi-likelihood F-tests. 
At gene-level, two different tests are provided. The first is converting exon-level p-values to gene-level p-values by Simes method.
The other is an F-test for differences between the exon log-fold-changes within each gene.
}

\author{Yunshun Chen and Gordon Smyth}

\examples{
# Gene exon annotation
Gene <- paste("G", 1:10, sep="")
Gene <- rep(Gene, each=10)
Exon <- paste("Ex", 1:10, sep="")
Gene.Exon <- paste(Gene, Exon, sep=".")
genes <- data.frame(GeneID=Gene, Gene.Exon=Gene.Exon)

design <- model.matrix(~c(0,0,0,1,1,1))
mu <- matrix(20, 100, 6)
mu[1,4:6] <- 200
counts <- matrix(rnbinom(600,mu=mu,size=20),100,6)

y <- DGEList(counts=counts, lib.size=rep(1e6,6), genes=genes)
gfit <- glmFit(y, design, dispersion=0.05)

ds <- diffSpliceDGE(gfit, geneid="GeneID")
topSpliceDGE(ds)
plotSpliceDGE(ds)
}

