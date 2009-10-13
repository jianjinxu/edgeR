\name{edgeR-package}
\docType{package}
\alias{edgeR}
\alias{edgeR-package}
\title{Empirical analysis of digital gene expression data in R}
\description{
edgeR is a library for the analysis of digital gene expression data arising from RNA sequencing technologies such as SAGE, CAGE or RNA-seq, with emphasis on testing for differential expression.

Particular strengths of the package include the ability to estimate biological variation between replicate libraries, and to conduct exact tests of significance which are suitable for small counts.
The package is able to make use of even minimal numbers of replicates.

A User's Guide is available as well as the usual help page documentation for each of the individual functions.

The library implements statistical methodology developed by Robinson and Smyth (2007, 2008).
}

\details{
\tabular{ll}{
Package: \tab edgeR\cr
Version: \tab 1.3.7\cr
Date: \tab 2009/10/13\cr
Depends: \tab R (>= 2.3.0), methods\cr
Imports: \tab limma\cr
Suggests: \tab MASS\cr
biocViews: \tab Bioinformatics, DifferentialExpression, Statistics\cr
License: \tab LGPL version 2 or newer\cr
}

Index:
\preformatted{
DGEList-class           Digital Gene Expression data - class
EBList-class            differential expression of Digital Gene
                        Expression data - class
SmoothList-class        differential expression of Digital Gene
                        Expression data - class
Tu102                   Raw data for several SAGE libraries from the
                        Zhang 1997 Science paper.
alpha.approxeb          Estimate the prior weight, alpha
approx.expected.info    Approximate of expected information (Fisher
                        information)
commonCondLogLikDerDelta
                        Common conditional log-likelihood in terms of
                        delta
condLogLikDerDelta      Conditional log-likelihood in terms of delta
condLogLikDerSize       Log-likelihood of the common dispersion for a
                        single equalized group.
de4DGE                  Compute moderated differential expression
                        scores for digital gene expression (DGE) data
de4DGEList-class        differential expression for Digital Gene
                        Expression data - class
deDGE                   Compute moderated differential expression
                        scores for digital gene expression (DGE) data
deDGEList-class         differential expression of Digital Gene
                        Expression data - class
edgeR-package           Introduction to the edgeR Package
equalizeLibSizes        Quantile adjustment to equalize library sizes
                        for a fixed value of the dispersion parameter
estimateCommonDisp      Estimates the negative binomial common
                        dispersion by maximizing the negative binomial
                        conditional common likelihood
estimateDispIter        Normalizes a dataset using quantile adjustment
                        and iteratively estimates the dispersion
                        parameter
estimatePs              Estimate expression proportions
estimateSmoothing       Estimate the prior weight and common dispersion
estimateTagwiseDisp     Maximizes the negative binomial common
                        likelihood
exactTest               An exact test for differences between two
                        negative binomial groups
exactTestNB             An exact test for differences between two
                        negative binomial groups
findMaxD2               Maximizes the negative binomial likelihood
getCounts               Extract table of counts from DGEList object
interpolateHelper       Quantile Adjustment interpolator
logLikDerP              Log-likelihood for proportion
maPlot                  plots log fold change versus log concentration
                        (or, M versus A) for count data
plotFC                  Plot of fold change against abundance for DGE
                        data
plotSmear               plots log fold change versus log concentration
                        (or, M versus A) for count data
q2qnbinom               Quantile to quantile mapping between
                        negative-binomial distributions
quantileAdjust          Normalizes a dataset by using a quantile
                        adjustment
readDGE                 Read and merge a set of files containing DGE
                        data
splitIntoGroups         Split the counts from a DGEList object
                        according to group
splitIntoGroupsPseudo   Split pseudocounts according to group
subsetting              Subset DGEist Objects
tau2.0.objective        Objective function for tau2
topTags                 Table of the Top Differentially Expressed Tags
weightedCondLogLikDerDelta
                        Weighted conditional log-likelihood in terms of
                        delta
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{edgeR} \tab edgeR User's Guide (source, pdf)\cr
}
}
\author{
Mark Robinson <mrobinson@wehi.edu.au>, Davis McCarthy <dmccarthy@wehi.edu.au>, Gordon Smyth
}

\references{
Robinson, MD, and Smyth, GK (2007). Moderated statistical tests for assessing differences in tag abundance.
\emph{Bioinformatics} 23, 2881-2887.

Robinson, MD and Smyth, GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data.
\emph{Biostatistics} 9, 321-332.
}
\keyword{package}
