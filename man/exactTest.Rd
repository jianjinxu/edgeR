\name{exactTest}
\alias{exactTest}

\title{An exact test for differences between two negative binomial groups}

\description{Carries out an exact test for differences between two negative binomial groups; calculations performed by \code{exactTestNB}}

\usage{ 
exactTest(object,pair=NULL, common.disp=TRUE)
}

\arguments{ 

\item{object}{a \code{DGEList} object, output of \code{estimateCommonDisp}, on which to compute Fisher-like exact statistics for the pair of groups specified}

\item{pair}{pair of groups to be compared; if \code{pair=NULL}, then the first two levels of \code{object$group} are compared}

\item{common.disp}{logical, if \code{TRUE}, then testing carried out using common dispersion for each tag/gene, if \code{FALSE} then tag-wise estimates of the dispersion parameter are used; default \code{TRUE}}
}

\value{a \code{deDGEList} object containing the following elements.
	\item{table}{a data frame containing the elements \code{logConc}, the log-average concentration/abundance for each tag in the two groups being compared, \code{logFC}, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, \code{p.value}, exact p-value for differential expression using the NB model}
	\item{comparison}{a vector giving the names of the two groups being compared} 
}


\author{Mark Robinson, Davis McCarthy}
\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(80,size=1,mu=10),nrow=20)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$counts)<-paste("tagno",1:nrow(d$counts),sep=".")

# estimate common dispersion and find differences in expression
d<-estimateCommonDisp(d)
de<-exactTest(d)
}
\keyword{algebra}
