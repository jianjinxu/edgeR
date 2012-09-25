\name{edgeR-package}
\docType{package}
\alias{edgeR}
\alias{edgeR-package}
\title{Empirical analysis of digital gene expression data in R}
\description{
edgeR is a package for the analysis of digital gene expression data arising from RNA sequencing technologies such as SAGE, CAGE, Tag-seq or RNA-seq, with emphasis on testing for differential expression.

Particular strengths of the package include the ability to estimate biological variation between replicate libraries, and to conduct exact tests of significance which are suitable for small counts.
The package is able to make use of even minimal numbers of replicates.

An extensive User's Guide is available, and can be opened by typing \code{edgeRUsersGuide()} at the R prompt.
Detailed help pages are also provided for each individual function.

The edgeR package implements original statistical methodology described in the publications below.
}

\author{
Mark Robinson <mrobinson@wehi.edu.au>, Davis McCarthy <dmccarthy@wehi.edu.au>, Yunshun Chen <yuchen@wehi.edu.au>, Aaron Lun <alun@wehi.edu.au>, Gordon Smyth
}

\references{
  Robinson MD and Smyth GK (2007). Moderated statistical tests for
  assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887

  Robinson MD and Smyth GK (2008). Small-sample estimation of negative
  binomial dispersion, with applications to SAGE data. \emph{Biostatistics},
  9, 321-332

  Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data. \emph{Bioinformatics} 26, 139-140

  McCarthy, DJ, Chen, Y, Smyth, GK (2012).
  Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
  \emph{Nucleic Acids Research} 40, 4288-4297.

  Lund, SP, Nettleton, D, McCarthy, DJ, Smyth, GK (2012).
  Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates.
  \emph{Statistical Applications in Genetics and Molecular Biology}. (Accepted 31 July 2012)
}
\keyword{package}
