\name{DGEGLM-class}
\docType{class}
\alias{DGEGLM-class}
\alias{show,DGEGLM-method}
\title{Digital Gene Expression Generalized Linear Model results - class}

\description{
A simple list-based class for storing results of a GLM fit
to each tag/gene in a DGE dataset.
}

\section{Slots/List Components}{
Objects of this class contain the following list components:

 \code{coefficients}: matrix containing the coefficients
  computed from fitting the model defined by the design matrix to each gene/tag in the dataset.

  \code{df.residual}: vector containing the residual degrees of freedom for the model fit to each tag/gene in the dataset.

  \code{deviance}: vector giving the deviance from the model fit to each tag/gene.

 \code{design}: design matrix for the full model from the likelihood
  ratio test.

  \code{offset}: scalar, vector or matrix of offset values to be included in the GLMs for each tag/gene.

  \code{samples}: data frame containing information about the samples comprising the dataset.

  \code{genes}: data frame containing information about the genes or tags for which we have DGE data (can be \code{NULL} if there is no information available).

  \code{dispersion}: scalar or vector providing the value of the dispersion parameter used in the negative binomial GLM for each tag/gene.

  \code{lib.size}: vector providing the effective library size for each sample in the dataset.

  \code{weights}: matrix of weights used in the GLM fitting for each tag/gene.

  \code{fitted.values}: the fitted (expected) values--here they are counts--from the GLM for each tag/gene.

  \code{abundance}: vector of gene/tag abundances (expression level), on the log2 scale, computed from the mean count for each gene/tag after scaling count by normalized library size.

}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGEGLM} objects also have a \code{show} method.
}

\author{Davis McCarthy}

\keyword{classes}

