\name{decideTestsDGE}
\alias{decideTestsDGE}
\title{Multiple Testing Across Genes and Contrasts}
\description{
Classify a series of related differential expression statistics as up, down or not significant.
A number of different multiple testing schemes are offered which adjust for multiple testing down the genes as well as across contrasts for each gene.
}
\usage{
decideTestsDGE(object, adjust.method="BH", p.value=0.05)
}
\arguments{
  \item{object}{\code{deDGElist} object, output from \code{exactTest},
  or \code{DGELRT} object, output from \code{DGELRT}, from which
  p-values for differential expression and log-fold change values may be extracted.}

 \item{adjust.method}{character string specifying p-value adjustment method.  Possible values are \code{"none"}, \code{"BH"}, \code{"fdr"} (equivalent to \code{"BH"}), \code{"BY"} and \code{"holm"}. See \code{\link[stats]{p.adjust}} for details.}

  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
}
\value{
An object of class \code{TestResults} (see \code{\link[limma:TestResults]{TestResults}}).
This is essentially a numeric matrix with elements \code{-1}, \code{0}
or \code{1} depending on whether each DE p-value is classified as
significant with negative log-fold change, not significant or
significant with positive log-fold change, respectively.
}
\details{
These functions implement multiple testing procedures for determining
whether each log-fold change in a matrix of log-fold changes should be considered significantly different from zero.

}
\seealso{
Adapted from \code{\link[limma:decideTests]{decideTests}} in the limma package.
}
\author{Davis McCarthy, Gordon Smyth}
\keyword{htest}
