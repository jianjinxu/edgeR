\name{exactTest}
\alias{exactTest}
\alias{exactTest.matrix}

\title{Exact Tests for Differences between Two Groups of Negative-Binomial Counts}

\description{Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts.}

\usage{
exactTest(object,pair=NULL,dispersion=NULL,common.disp=TRUE)
exactTest.matrix(y1,y2,r=NULL,dispersion=0,all.zeros=NULL)
}

\arguments{ 

\item{object}{a \code{DGEList} object containing pseudo-counts and group IDs.}

\item{pair}{vector of length two, either numeric or character, providing the pair of groups to be compared; if a character vector, then should be the names of two groups (e.g. two levels of \code{object$samples$group}); if numeric, then groups to be compared are chosen by finding the levels of \code{object$samples$group} corresponding to those numeric values and using those levels as the groups to be compared; if \code{NULL}, then first two levels of \code{object$samples$group} (a factor) are used. Note that the first group listed in the pair is the baseline for the comparison---so if the pair is \code{c("A","B")} then the comparison is \code{B - A}, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change).}

\item{dispersion}{optional vector either of length 1 or the same length as the number of tags. If not \code{NULL} (default), then the supplied value(s) will be used as the dispersion parameter for calculating p-values for differential expression. If \code{NULL}, then either the common or tagwise dispersion estimates from the \code{DGEList} object will be used, according to the value of \code{common.disp}. If \code{dispersion} is zero, then p-values are equivalent to exact Poisson rather than NB p-values.
For \code{exactTest.matrix} equal to \code{1/r} and ignored if \code{r} is supplied.}

\item{common.disp}{logical, if \code{TRUE}, then testing carried out using common dispersion for each tag/gene, if \code{FALSE} then tag-wise estimates of the dispersion parameter are used; default \code{TRUE}.}

\item{y1}{numeric matrix of counts for one of the two given experimental groups to be tested for differences. Libraries are assumed to be equal in size - e.g. adjusted pseudocounts from the output of \code{\link{equalizeLibSizes}}.}

\item{y2}{numeric matrix of counts for one of the two given experimental groups to be tested for differences. Libraries are assumed to be equal in size - e.g. adjusted pseudocounts from the output of \code{\link{equalizeLibSizes}}. Must have the same number of rows as \code{y1}.}

\item{r}{vector of negative binomial \code{size} parameter values (\code{size} = \code{1/phi} where \code{phi} is the dispersion parameter in the NB model); if \code{r} is of length 1, then a common value of the dispersion is used for all transcripts, otherwise, must be a vector with length equal to the number of rows of \code{y1} and \code{y2}.}

\item{all.zeros}{logical vector indicating which tags have zero total counts over all libraries.   The output p-value is set to \code{1} for these tags.  If not provided, is determined from the data.}

}

\value{
\code{exactTest.matrix} produces a numeric vector of genewise p-values with length equal to the number of transcripts, taken to be the number of rows of \code{y1}.

\code{exactTest} produces an object of class \code{DGEExact} containing the following components:
\item{table}{a data frame containing the elements \code{logConc}, the log-average concentration/abundance for each tag in the two groups being compared, \code{logFC}, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, \code{p.value}, exact p-value for differential expression using the NB model}
\item{comparison}{a vector giving the names of the two groups being compared} 
\item{genes}{a data frame containing information about each transcript; taken from \code{object} and can be \code{NULL}}
}

\details{
This function can be used to test for differential expression between groups when the library sizes are equal.
In the context of the edgeR package, it is applied after equalizing the library sizes using \code{equalizeLibSizes}, in other words, it is applied to quantile-adjusted pseudo-counts.
The dispersion may be different for genes, but is assumed equal for all the counts for a given gene.

The exact test uses the likelihood ratio test to generate a test statistics, then computes the exact p-value for this statistic from the conditional distribution of the group means conditional on the grand mean.
This test can be viewed as a generalization of the well-known exact binomial test, but generalized to overdispersed counts.

For each transcript, conditioning on the total sum of counts within each group and the total sum of counts across all groups allows us to construct an exact test for differences between two group. The conditional distribution for the sum of counts in a group is known (given the values for the mean counts, \code{mus}, and the dispersion parameter, 1/\code{r}), exact p-values can be computed by summing over all sums of counts that have a probability less than the probability under the null hypothesis of the observed sum of counts.

\code{exactTest.matrix} is the function that actually computes the exact p-values.
\code{exactTest} is intended to have a more object-orientated flavor as it produces objects containing all the necessary components for downstream analysis.

}

\references{
Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. \emph{Biostatistics}, 9, 321-332
}

\author{Mark Robinson, Davis McCarthy}

\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(80,size=1,mu=10),nrow=20)
d <- DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$counts) <- paste("tagno",1:nrow(d$counts),sep=".")

# estimate common dispersion and find differences in expression
d <- estimateCommonDisp(d)
de <- exactTest(d)
topTags(de)

# example using exactTest.matrix directly
y <- matrix(rnbinom(15,mu=10,size=1.5),5,3)
group <- factor(c(1,1,2))
p.value <- exactTest.matrix(y[,1:2],y[,3,drop=FALSE],r=1.5)
}

\seealso{
\code{\link{equalizeLibSizes}, \code{\link{binomTest}}
}

\keyword{algebra}

