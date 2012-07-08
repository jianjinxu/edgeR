\name{predFC}
\alias{predFC}
\alias{predFC.DGEList}
\alias{predFC.default}

\title{Predictive log fold changes}

\description{Computes estimated coefficients for a generalised linear model with log-fold-changes shrunk towards zero.}

\usage{
\S3method{predFC}{DGEList}(y, design=NULL, prior.count.total=0.5, offset=NULL, dispersion=NULL)
\S3method{predFC}{default}(y, design=NULL, prior.count.total=0.5, offset=log(colSums(y)), dispersion=0)
}

\arguments{
\item{y}{a matrix of counts or a \code{DGEList} object} 

\item{design}{the design matrix for the experiment}

\item{prior.count.total}{the total prior number of counts to be added to the data. Larger values produce more shrinkage.}

\item{offset}{numeric vector or matrix giving the offset in the log-linear model predictor, as for \code{\link{glmFit}}. Usually equal to log library sizes.}

\item{dispersion}{the negative binomial dispersion}
}

\details{
This function adds \code{prior.count.total} counts to each row of \code{y}.
The counts are added is such a way that any log-fold-change that was zero prior to augmentation remains zero and non-zero log-fold-changes are shrunk towards zero.

The prior counts can be viewed as equivalent to a prior belief that the log-fold changes are small, and the output can be viewed as posterior log-fold-changes from this Bayesian viewpoint.
The output coefficients are called \emph{predictive} log fold-changes because, depending on the prior, they may be a better prediction of the true log fold-changes than the raw estimates.

This is done small count is added to each library in proportion to the library sizesA larger amount is added to counts from larger libraries, so that any log-fold-change that was zero prior to augmentation remains zero.
count is added to The specific 
stimates the predictive or posterior log-fold-changes for count data.
. If there are 2 groups in the experiment, n=2 for each group, the total prior count is 1, and the library sizes are equal, then in effect 0.5 of a count is added to each group, or 0.25 to each library. This prior count is the same for all genes or tags in the data, with the result that genes with low counts will be dampened more severely and genes with a large number of counts in each library will hardly be affected by the addition of a small count to each group. 

In order to get the predictive log-fold-changes, a generalised linear model is fitted to the augmented data, and the coefficients outputted in the form of a matrix. 

If \code{offset=NULL}, the offset will be taken from the \code{DGEList} object or computed from the column sums.

If \code{dispersion=NULL}, the dispersion used for the glm will be dependent on what is in the DGEList object; it is prioritised in the following manner: tagwise, trended, common and finally if no dispersion estimate is found it will set the dispersion to 0.

If \code{design=NULL}, then the function returns a matrix of the same size as \code{y} containing log2 counts-per-million, with zero values for the counts avoided.
This equivalent to choosing \code{design} to be the identity matrix with the same number of columns as \code{y}.
}

\value{
Numeric matrix of linear model coefficients (if \code{design} is given) or logCPM (if \code{design=NULL}) on the log2 scale.
}

\author{Belinda Phipson, Gordon Smyth}

\examples{
# generate counts for a two group experiment with n=2 in each group and 100 genes
dispersion <- 0.1
y <- matrix(rnbinom(400,size=1/dispersion,mu=4),nrow=100)
y <- DGEList(y,group=c(1,1,2,2))
design <- model.matrix(~group, data=y$samples)

#estimate the predictive log fold changes
predlfc<-predFC(y,design,dispersion=dispersion,prior.count=4)
logfc <- predFC(y,design,dispersion=dispersion,prior.count=0)
logfc.truncated <- pmax(pmin(logfc,100),-100)

#plot predFC's vs logFC's
plot(predlfc[,2],logfc.truncated[,2],xlab="Predictive log fold changes",ylab="Raw log fold changes")
abline(a=0,b=1)
}

\seealso{
\code{\link{glmFit}}, \code{\link{exactTest}}
}



