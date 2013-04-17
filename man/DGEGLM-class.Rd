\name{DGEGLM-class}
\docType{class}
\alias{DGEGLM-class}
\alias{show,DGEGLM-method}
\title{Digital Gene Expression Generalized Linear Model results - class}

\description{
A list-based S4 class for storing results of a GLM fit to each gene in a DGE dataset.
}

\section{List Components}{
For objects of this class, rows correspond to genomic features and columns to coefficients in the linear model.
The genomic features are called genes, but in reality might correspond to transcripts, tags, exons etc.

Objects of this class contain the following list components:
\tabular{ll}{
 \code{coefficients } \tab matrix containing the coefficients computed from fitting the model defined by the design matrix to each gene in the dataset.\cr

  \code{df.residual } \tab vector containing the residual degrees of freedom for the model fit to each gene in the dataset.\cr

  \code{deviance } \tab vector giving the deviance from the model fit to each gene.\cr

  \code{design } \tab design matrix for the full model from the likelihood ratio test.\cr

  \code{offset } \tab scalar, vector or matrix of offset values to be included in the GLMs for each gene.\cr

  \code{samples } \tab data frame containing information about the samples comprising the dataset.\cr

  \code{genes } \tab data frame containing information about the genes or tags for which we have DGE data (can be \code{NULL} if there is no information available).\cr

  \code{dispersion } \tab scalar or vector providing the value of the dispersion parameter used in the negative binomial GLM for each gene.\cr

  \code{lib.size } \tab vector providing the effective library size for each sample in the dataset.\cr

  \code{weights } \tab matrix of weights used in the GLM fitting for each gene.\cr

  \code{fitted.values } \tab the fitted (expected) values from the GLM for each gene.\cr

  \code{AveLogCPM } \tab numeric vector giving average log2 counts per million for each gene.
}
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.

The dimensions, row names and column names of a \code{DGEGLM} object are defined by those of the dataset, see \code{\link{dim.DGEGLM}} or \code{\link{dimnames.DGEGLM}}.

\code{DGEGLM} objects can be subsetted, see \code{\link{subsetting}}.

\code{DGEGLM} objects also have a \code{show} method so that printing produces a compact summary of their contents.

}
\author{edgeR team. First created by Davis McCarthy.}

\seealso{
Other classes defined in edgeR are \code{\link{DGEList-class}}, \code{\link{DGEExact-class}}, \code{\link{DGELRT-class}}, \code{\link{TopTags-class}}
}

\keyword{classes}

