\name{predFC}
\alias{predFC}
\alias{predFC.DGEList}
\alias{predFC.default}

\title{Predictive log fold changes for RNASeq data}

\description{Estimates the predictive log fold changes for a given prior weight using generalised linear models.}

\usage{
\S3method{predFC}{DGEList}(y, design, prior.count.total=0.5, offset=NULL, dispersion=NULL)
\S3method{predFC}{default}(y, design, prior.count.total=0.5, offset=log(colSums(y)), dispersion=0)
}

\arguments{ 
\item{y}{a \code{DGEList} object} 

\item{design}{the design matrix for the experiment}

\item{prior.count.total}{the total number of counts to be augmented to the data}

\item{offset}{numeric vector or matrix giving the offset in the log-linear model predictor. Usually equal to log library sizes}

\item{dispersion}{the dispersion estimate for the count data}

}

\details{This function estimates the predictive or posterior log-fold-changes for RNASeq or any count-based data. A small count is added to each library in proportion to the library sizes. If there are 2 groups in the experiment, n=2 for each group, the total prior count is 1, and the library sizes are equal, then in effect 0.5 of a count is added to each group, or 0.25 to each library. This prior count is the same for all genes or tags in the data, with the result that genes with low counts will be dampened more severely and genes with a large number of counts in each library will hardly be affected by the addition of a small count to each group. 

In order to get the predictive log-fold-changes, a generalised linear model is fitted to the augmented data, and the coefficients outputted in the form of a matrix. 

If \code{offset=NULL}, the offset used in the glm will be the log of the library sizes.

If \code{dispersion=NULL}, the dispersion used for the glm will be dependent on what is in the DGEList object; it is prioritised in the following manner: tagwise, trended, common and finally if no dispersion estimate is found it will set the dispersion to 0.

}

\value{
Numeric matrix of log-fold-changes (natural log).  Each column corresponds to the relevant column of the design matrix.
}

\author{Belinda Phipson, Gordon Smyth}

\examples{
# generate counts from a negative binomial distribution for a two group experiment with n=2 in each group
y<-matrix(rnbinom(400,size=1,mu=10),nrow=100)
y<-DGEList(y,group=c(1,1,2,2))
design<-model.matrix(~y$samples$group)

# apply TMM normalisation
y<-calcNormFactors(y)

# estimate the common dispersion 
y <- estimateGLMCommonDisp(y,design)

# fit a glm to find differentially expressed genes
glm<-glmFit(y,design)
results<-glmLRT(y,glm,coef=2)

#estimate the predictive log fold changes
pfc<-predFC(y,design)

#plot predFC's vs logFC's
plot(pfc[,2],results$table$logFC,xlab="Predictive log fold changes",ylab="Raw log fold changes",pch=16,cex=0.8)
abline(a=0,b=1)

}

\seealso{
\code{\link{glmFit}}, \code{\link{glmLRT}} for generalised linear model fitting

\code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTrendedDisp}}, \code{\link{estimateGLMTagwiseDisp}} for estimating dispersions in the context of generalised linear models

\code{\link{estimateCommonDisp}}, \code{\link{estimateTagwiseDisp}} for estimating dispersions when the design of the experiment is a simple one-way layout.

\code{\link{calcNormFactors}} for TMM normalisation
}



