\name{deDGE}

\alias{deDGE}

\title{Compute moderated differential expression scores for digital gene expression (DGE) data}

\description{Runs weighted likelihood calculation for moderated estimates of dispersion, and tests for differences in 'tag' abundance between groups}

\usage{ 
deDGE(object,alpha=500,doPoisson=FALSE,verbose=TRUE)
}

\arguments{ 
\item{object}{\code{DGEList} containing elements \code{data} (matrix: rows-tags, columns-libraries), \code{lib.size} \code{group} indicating class} 

\item{alpha}{weight to put on the individual tag's likelihood}

\item{doPoisson}{logical, whether to fit Poisson model instead of Negative Binomial, default \code{FALSE}}

\item{verbose}{logical, whether to write comments, default \code{TRUE}}
}

\value{ \code{deDGEList} with elements \code{lr} (likelihood ratio test), \code{r} (estimates of 1/overdispersion), \code{ps} (list containing proportion estimates)}
\author{Mark Robinson}
\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))

# find alpha and call main procedure to find differences
alpha<-alpha.approxeb(d)
ms<-deDGE(d,alpha=alpha$alpha)
}

\references{
Robinson MD, Smyth GK. 'Small-sample estimation of negative binomial dispersion, with applications to SAGE data.' Biostatistics. 2008 Apr;9(2):321-32.

Robinson MD, Smyth GK. 'Moderated statistical tests for assessing differences in tag abundance.' Bioinformatics. 2007 Nov 1;23(21):2881-7.
}

\keyword{algebra}
