\name{goana.DGELRT}
\alias{goana.DGELRT}
\title{Gene Ontology Analysis of Differentially Expressed Genes}
\description{
Fisher's exact test or Wallenius' noncentral hypergeometric test on differentially expressed genes for over-representation of gene ontology (GO) terms.
}
\usage{
\method{goana}{DGELRT}(de, geneid = rownames(de), FDR = 0.05, 
      species = "Hs", trend = FALSE, \dots)
}
\arguments{
  \item{de}{\code{DGELRT} object.}
  \item{geneid}{Entrez Gene identifiers. Either a vector of length \code{nrow(de)} or the name of the column of \code{de$genes} containing the Entrez Gene IDs.}
  \item{FDR}{numeric. False discovery rate of differentially expressed genes less than this cutoff.}
  \item{species}{character string specifying the species. Possible values are \code{"Hs"}, \code{"Mm"}, \code{"Rn"} or \code{"Dm"}.}
  \item{trend}{adjust analysis for gene length or abundance?
  Can be logical, or a numeric vector of covariate values, or the name of the column of \code{de$genes} containing the covariate values.
  If \code{TRUE}, then \code{de$AveLogCPM} is used as the covariate.}
  \item{\dots}{further arguments.}
}
\details{
Performs a Gene Ontology enrichment analysis for one for more gene lists using the appropriate Bioconductor organism package.
The gene lists must be supplied as Entrez Gene IDs.

If applied to a linear model fit, then the differentially expressed genes are extracted automatically from the fit object, and analyses are done separately for the up and down significant genes.

If \code{trend=FALSE}, the function computes one-sided hypergeometric tests.
If \code{trend=TRUE} or a covariate is supplied, then a trend is fitted to the differential expression results and the method of Young et al (2010) is used to adjust for this trend.
The adjusted test uses Wallenius' noncentral hypergeometric distribution.
}
\value{
  A dataframe with rows for GO IDs and the following possible columns:
  \item{Term}{GO terms.}
  \item{Ont}{GO term ontology from BP, CC and MF.}
  \item{N}{Number of genes in the GO term.}
  \item{Up}{Number of up regulated differentially expressed genes for testing.}
  \item{Down}{Number of down regulated differentially expressed genes for testing.}
  \item{DE1}{Number of differentially expressed genes for testing.}
  \item{P.Up}{P value from hypergeometric test or Wallenius' noncentral hypergeometric test of up regulated differentially expressed genes for over representation of GO terms.}
  \item{P.Down}{P value from hypergeometric test or Wallenius' noncentral hypergeometric test of down regulated differentially expressed genes for over representation of GO terms.}
  \item{P.DE1}{P value from hypergeometric test or Wallenius' noncentral hypergeometric test of differentially expressed genes for over representation of GO terms.}
}

\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{goana}}
}
\author{Gordon Smyth and Yifang Hu}
\examples{

\dontrun{
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
go <- goana(lrt)
}

}

\keyword{gene set test}
