\name{weightedComLik}
\alias{weightedComLik}
\alias{weightedComLikMA}

\title{Weighted Common Log-Likelihood}

\description{Allow a flexible approach to accounting for a potential dependence of the dispersion on the abundance (expression level) of tags/genes by calculating a weighted 'common' log-likelihood for each gene.}

\usage{
weightedComLik(object,l0,prop.used=0.25)
weightedComLikMA(object,l0,prop.used=0.05)
}

\arguments{ 

\item{object}{\code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{l0}{matrix of the conditional log-likelihood evaluated at a variety of values for the dispersion (on the delta scale, \code{phi/(1 + phi)}) for each tag/gene. The matrix has number of rows equal to the number of tags/genes and number of columns equal to the number of grid values (between 0 and 1) for the dispersion at which the conditional log-likelihood is evaluated.}

\item{prop.used}{scalar giving the proportion of tags/genes in the whole dataset to use in computing the weighted common log-likelihood for each tag/gene. Default value is \code{0.25}, i.e. a quarter of the tags/genes in the dataset, for \code{weightedComLik} and \code{0.05} for \code{weightedComLikMA}.}

}

\value{matrix of weighted common log-likelihood values computed for each gene at each grid value for the dispersion. The matrix returned has the same dimensions as l0.}

\details{
Genes are ordered based on abundance (expression level) and for a given gene, a proportion of the genes close to it are used to compute the common log-likelihood with decreasing weight given to the genes further from the given gene. Weighting is done using the tricube weighting function for \code{weightedComLik}. Computation can be slow relative to other functions in \code{edgeR}, especially if the number of genes or the number of grid values (i.e. the dimensions of l0) are large. \code{weightedComLikMA} uses a moving average to do the weighting (using \code{\link{movingAverageByCol}}) and so is much faster than \code{weightedComLik}.
}


\author{Davis McCarthy}
