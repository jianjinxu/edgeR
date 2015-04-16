\title{Top table of differentially spliced genes or exons}
\name{topSpliceDGE}
\alias{topSpliceDGE}
\description{
Top table ranking the most differentially spliced genes or exons.
}
\usage{
topSpliceDGE(lrt, level="gene", gene.test="Simes", number=10, FDR=1)
}
\arguments{
  \item{lrt}{\code{DGELRT} object produced by \code{diffSpliceDGE}.}
  \item{level}{character string, should the table be by \code{"exon"} or by \code{"gene"}.}
  \item{gene.test}{character string, choice for the gene-level p-values. Possible values are "Simes" and "F".}
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this cutoff.}
}

\details{
Ranks exons or genes by p-values.
}

\value{A data.frame with any annotation columns found in \code{fit} plus the following columns
  \item{NExons}{number of exons if \code{level="gene"}}
  \item{Gene.Exon}{exon annotation if \code{level="exon"}}
  \item{logFC}{log-fold change of one exon vs all the exons for the same gene (if \code{level="exon"})}
  \item{F}{F-statistics for exons if \code{level="exon"}}
  \item{P.Value}{p-value}
  \item{FDR}{false discovery rate}
}

\author{Yunshun Chen and Gordon Smyth}

\examples{# See \code{\link{diffSpliceDGE}}}
