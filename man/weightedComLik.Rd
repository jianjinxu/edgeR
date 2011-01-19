\name{weightedComLik}
\alias{weightedComLik}

\title{Weighted Common Log-Likelihood}

\description{Allow a flexible approach to accounting for a potential dependence of the dispersion on the abundance (expression level) of tags/genes by calculating a weighted 'common' log-likelihood for each gene.}

\usage{
weightedComLik(object,l0,prop.used=0.25)
}

\arguments{ 

\item{object}{\code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{l0}{matrix of the conditional log-likelihood evaluated at a variety of values for the dispersion (on the delta scale, \code{phi/(1 + phi)}) for each tag/gene. The matrix has number of rows equal to the number of tags/genes and number of columns equal to the number of grid values (between 0 and 1) for the dispersion at which the conditional log-likelihood is evaluated.}

\item{prop.used}{scalar giving the proportion of tags/genes in the whole dataset to use in computing the weighted common log-likelihood for each tag/gene. Default value is 0.25, i.e. a quarter of the tags/genes in the dataset.}

}

\value{matrix of weighted common log-likelihood values computed for each gene at each grid value for the dispersion. The matrix returned has the same dimensions as l0.}

\details{
Genes are ordered based on abundance (expression level) and for a given gene, a proportion of the genes close to it are used to compute the common log-likelihood with decreasing weight given to the genes further from the given gene. Weighting is done using the tricube weighting function. Computation can be slow relative to other functions in \code{edgeR}, especially if the number of genes or the number of grid values (i.e. the dimensions of l0) are large.
}


\author{Davis McCarthy}

\examples{
counts<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=counts,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
d<-estimateCommonDisp(d)
ntags<-nrow(d$counts)
y<-splitIntoGroups(new("DGEList",list(counts=d$pseudo.alt,samples=d$samples)))
grid.vals<-seq(0.001,0.999,length.out=10)
l0<-0
for(i in 1:length(y)) {
  	 l0<-condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
}
m0 <- ntags*weightedComLik(d,l0,prop.used=0.25) # Weights sum to 1, so need to multiply by number of tags to give this the same weight overall as the regular common likelihood
}
\keyword{file}

