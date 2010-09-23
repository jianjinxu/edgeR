\name{DGELRT-class}
\docType{class}
\alias{DGELRT-class}
\alias{show,DGELRT-method}
\title{Digital Gene Expression Likelihood Ratio Test data and results - class}

\description{
A simple list-based class for storing results of a GLM-based
differential expression analysis for DGE data, with evidence for
differential expression assessed using a likelihood ratio test.
}

\section{Slots/List Components}{
Objects of this class contain the following list components:

  \code{table}: data frame containing the log-concentration (i.e. expression level), the log-fold change in expression between the two groups/conditions and the exact p-value for differential expression, for each tag.
  
  \code{coefficients.full}: matrix containing the coefficients
  computed from fitting the full model (fit using \code{glmFit} and a
  given design matrix) to each gene/tag in the dataset.

  \code{coefficients.null}: matrix containing the coefficients
  computed from fitting the null model to each gene/tag in the
  dataset.  The null model is the model to which the full model is
  compared, and is fit using \code{glmFit} and dropping selected
  column(s) (i.e. coefficient(s)) from the design matrix for the full model.

  \code{design}: design matrix for the full model from the likelihood
  ratio test.

  \code{...}: if the argument \code{y} to \code{glmLRT} (which
  produces the \code{DGELRT} object) was itself a \code{DGEList} object, then
  the \code{DGELRT} will contain all of the elements of \code{y},
  except for the table of counts and the table of pseudocounts.

}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGELRT} objects also have a \code{show} method.
}

\author{Davis McCarthy}

\keyword{classes}

