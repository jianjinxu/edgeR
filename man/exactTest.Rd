\name{exactTest}
\alias{exactTest}
\alias{exactTest.matrix}

\title{An Exact Test for Differences between Two Negative Binomial Groups}

\description{Carry out an exact test for differences between two negative binomial groups, based on conditioning on sums of (quantile-adjusted pseudo-)counts; calculations performed by \code{exactTest.matrix}}

\usage{ 
exactTest(object,pair=NULL, common.disp=TRUE)
exactTest.matrix(pseudo,group,pair=NULL,mus,r) 

}

\arguments{ 

\item{object}{a \code{DGEList} object, output of \code{estimateCommonDisp}, on which to compute Fisher-like exact statistics for the pair of groups specified}

\item{pair}{vector of length two providing the pair of groups to be compared; if \code{NULL}, then first two levels of \code{group} (as a factor) are used}

\item{common.disp}{logical, if \code{TRUE}, then testing carried out using common dispersion for each tag/gene, if \code{FALSE} then tag-wise estimates of the dispersion parameter are used; default \code{TRUE}}

\item{pseudo}{numeric matrix of counts (e.g. quantile adjusted pseudocounts) on which to compute exact statistics; rows are transcripts and columns are samples/libraries}

\item{group}{vector or factor indicating to which group the samples/libraries belong; must be same length as \code{ncol(pseudo)}}

\item{mus}{vector of count means for each transcript under the null hypothesis (of no difference between groups)}

\item{r}{vector of negative binomial \code{size} parameter values (\code{size} = \code{1/phi} where \code{phi} is the dispersion parameter in the NB model); if \code{r} is of length 1, then a common value of the dispersion is used for all transcripts, otherwise, must be a vector with length equal to the number of rows of \code{pseudo}.  If you want to run a Poisson test, set r very large (e.g. 1000)}

}

\value{\code{exactTest} produces an object of class \code{deDGEList} containing the following elements.
	\item{table}{a data frame containing the elements \code{logConc}, the log-average concentration/abundance for each tag in the two groups being compared, \code{logFC}, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, \code{p.value}, exact p-value for differential expression using the NB model}
	\item{comparison}{a vector giving the names of the two groups being compared} 
	\item{genes}{a data frame containing information about each transcript; taken from \code{object} and can be \code{NULL}}
\code{exactTest.matrix} produces a numeric vector of exact p-values with length equal to the number of transcripts, taken to be the number of rows of \code{pseudo}.
}

\details{
For each transcript, conditioning on the total sum of counts within each group and the total sum of counts across all groups allows us to construct an exact test for differences between two group. The conditional distribution for the sum of counts in a group is known (given the values for the mean counts, \code{mus}, and the dispersion parameter, 1/\code{r}), exact p-values can be computed by summing over all sums of counts that have a probability less than the probability under the null hypothesis of the observed sum of counts.

\code{exactTest.matrix} is the function that actually computes the exact p-values.
\code{exactTest} is intended to have a more object-orientated flavor as it produces objects containing all the necessary components for downstream analysis.

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

# example using exactTest.matrix directly
y<-matrix(rnbinom(20,mu=10,size=1.5),nrow=5)
group<-factor(c(1,1,2,2))
mus<-rep(10,5)
f<-exactTest.matrix(y,group,pair=c(1,2),mus,r=1.5)
}

\seealso{
Computing p-values for differential expression for each transcript between two (only) digital gene expression libraries can also be done using the \code{sage.test} function in the \code{statmod} package.
}

\keyword{algebra}
