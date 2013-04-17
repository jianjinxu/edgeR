\name{topTags}
\alias{topTags}
\alias{TopTags-class}
\alias{show,TopTags-method}
\alias{[.TopTags}

\title{Table of the Top Differentially Expressed Tags}

\description{Extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change.}

\usage{
topTags(object, n=10, adjust.method="BH", sort.by="PValue")
}

\arguments{ 
\item{object}{a \code{\link[edgeR:DGEList-class]{DGEExact}} object (output from \code{exactTest}) or a \code{\link[edgeR:DGELRT-class]{DGELRT}} object (output from \code{glmLRT}), containing the (at least) the elements \code{table}: a data frame containing the log-concentration (i.e. expression level), the log-fold change in expression between the two groups/conditions and the p-value for differential expression, for each tag. If it is a \code{DGEExact} object, then \code{topTags} will also use the \code{comparison} element, which is a vector giving the two experimental groups/conditions being compared. The object may contain other elements that are not used by \code{topTags}.}

\item{n}{scalar, number of tags to display/return}

\item{adjust.method}{character string stating the method used to adjust p-values for multiple testing, passed on to \code{p.adjust}}

\item{sort.by}{character string, should the top tags be sorted by p-value (\code{"PValue"}), by absolute log-fold change (\code{"logFC"}), or not sorted (\code{"none"}).}
}

\value{

an object of class \code{TopTags} containing the following elements for the top \code{n} most differentially expressed tags as determined by \code{sort.by}:

\item{table}{a data frame containing the elements \code{logFC}, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, \code{logCPM}, the log-average concentration/abundance for each tag in the two groups being compared, \code{PValue}, exact p-value for differential expression using the NB model, \code{FDR}, the p-value adjusted for multiple testing as found using \code{p.adjust} using the method specified.}
\item{adjust.method}{character string stating the method used to adjust p-values for multiple testing.}
\item{comparison}{a vector giving the names of the two groups being compared.}
\item{test}{character string stating the name of the test.}

The dimensions, row names and column names of a \code{TopTags} object are defined by those of \code{table}, see \code{\link{dim.TopTags}} or \code{\link{dimnames.TopTags}}.

\code{TopTags} objects also have a \code{show} method so that printing produces a compact summary of their contents.

}

\author{Mark Robinson, Davis McCarthy, Gordon Smyth}
\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(80,size=1,mu=10),nrow=20)
d <- DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$counts) <- paste("tag",1:nrow(d$counts),sep=".")

# estimate common dispersion and find differences in expression
# here we demonstrate the 'exact' methods, but the use of topTags is
# the same for a GLM analysis
d <- estimateCommonDisp(d)
de <- exactTest(d)

# look at top 10
topTags(de)
# Can specify how many tags to view
tp <- topTags(de, n=15)
# Here we view top 15
tp
# Or order by fold change instead
topTags(de,sort.by="logFC")
}

\references{
Robinson MD, Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. \emph{Biostatistics} 9, 321-332.

Robinson MD, Smyth GK (2007). Moderated statistical tests for assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887.
}

\seealso{
\code{\link{exactTest}}, \code{\link{glmLRT}}, \code{\link{p.adjust}}.

Analogous to \code{\link[limma:toptable]{topTable}} in the limma package.
}

\keyword{algebra}

