\name{SmoothList-class}
\docType{class}
\alias{SmoothList-class}
\alias{show,SmoothList-method}

\title{differential expression of Digital Gene Expression data - class}

\description{
A simple list-based class for storing results of the approximate empirical Bayes rule parameters, which determine the extent to which tagwise estimates of the dispersion are squeezed towards a common value.
}

\section{Slots/List Components}{
Objects of this class contain the following list components:
  \code{prior.n}: numeric scalar weight given to the common likelihood when estimating tag/genewise dispersions
  \code{common.dispersion}: numeric scalar estimate of common dispersion on the phi scale.
  \code{common.dispersion.delta}: numeric scalar estimate of common dispersion on the delta scale (delta = \code{phi/(1+phi)}.
  \code{sigma2.0.est}: numeric scalar \code{sigma_0^2} estimate.
  \code{alpha}: numeric scalar alpha estimate.
  \code{scores}: numeric scalar (likelihood) score.
  \code{infos}: numeric vector containing the (likelihood) information for each tag.
  \code{quantileAdjusted}: list from output of \code{equalizeLibSizes}).
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{SmoothList} objects also have a \code{show} method.
}

\author{Mark Robinson, Davis McCarthy}

\keyword{classes}
