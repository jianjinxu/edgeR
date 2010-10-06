\name{meanvar}
\alias{binMeanVar}
\alias{pooledVar}
\alias{plotMeanVar}

\title{Explore the mean-variance relationship for DGE data}

\description{Appropriate modelling of the mean-variance relationship in DGE data is important for making inferences about differential expression. Here are functions to compute tag/gene means and variances, as well at looking at these quantities when data is binned based on overall expression level.}


\usage{ 
plotMeanVar(object, meanvar=NULL, show.raw.vars=FALSE, show.tagwise.vars=FALSE, dispersion.method="coxreid", scalar=NULL, NBline=FALSE, nbins=100, ...) 
binMeanVar(x, conc=NULL, group, nbins=100)
pooledVar(y,group)
}

\arguments{ 
\item{object}{\code{DGEList} object containing the raw data and dispersion value. According the method desired for computing the dispersion, either \code{CRDisp} or \code{estimateCommonDisp} and (possibly) \code{estimateTagwiseDisp} should be run on the \code{DGEList} object before using \code{plotMeanVar}.}

\item{meanvar}{list (optional) containing the output from \code{binMeanVar} or the returned value of \code{plotMeanVar}. Providing this object as an argument will save time in computing the tag/gene means and variances when producing a mean-variance plot. }

\item{show.raw.vars}{logical, whether or not to display the raw (pooled) gene/tag variances on the mean-variance plot. Default is \code{FALSE}.}

\item{show.tagwise.vars}{logical, whether or not to display the estimated genewise/tagwise variances on the mean-variance plot. Default is \code{FALSE}.}

\item{dispersion.method}{character string giving the method that has been used to estimate the common and tagwise dispersion values used to calculate the estimated variances. Default is \code{"coxreid"} indicating that the Cox-Reid method for GLMs has been used to compute the dispersions; other option is \code{"qcml"} to indicate that conditional inference methods (e.g. \code{estimateCommonDisp} and \code{estimateTagwiseDisp} were used.}

\item{scalar}{vector (optional) of scaling values to divide counts by. Would expect to have this the same length as the number of columns in the count matrix (i.e. the number of libraries).}

\item{NBline}{logical, whether or not to add a line on the graph showing the mean-variance relationship for a NB model with common dispersion.}

\item{nbins}{scalar giving the number of bins (formed by using the quantiles of the genewise mean expression levels) for  which to compute average means and variances for exploring the mean-variance relationship. Default is \code{100} bins}

\item{\dots}{further arguments passed on to \code{plot}}

\item{x}{matrix of count data, with rows representing tags/genes and columns representing samples}

\item{conc}{vector (optional) of values for the concentration (i.e. abundance) of each tag}

\item{group}{factor giving the experimental group or condition to which each sample (i.e. column of \code{x} or element of {y}) belongs}

\item{y}{vector of count data}

}

\value{ \code{plotMeanVar} produces a mean-variance plot for the DGE data using the options described above. \code{plotMeanVar} and \code{binMeanVar} both return a list with the following components:
	\item{avemeans}{vector of the average expression level within each bin of genes}
	\item{avevars}{vector of the average raw pooled gene-wise variance within each bin of genes}
	\item{bin.means}{list containing the average (mean) expression level for genes divided into bins based on amount of expression}
	\item{bin.vars}{list containing the pooled variance for genes divided into bins based on amount of expression}
	\item{means}{vector giving the mean expression level for each gene}
	\item{vars}{vectore giving the pooled variance for each gene}

\code{pooledVar} returns a scalar for the pooled variance of the given data vector.

}

\details{
This function is useful for exploring the mean-variance relationship in the data. Raw variances are, for each gene, the pooled variance of the counts from each sample, divided by a scaling factor (by default the effective library size). The function will plot the average raw variance for tags split into \code{nbins} bins by overall expression level. This provides a useful summary of how the variance of the gene counts change with respect to average expression level (abundance). A line showing the Poisson mean-variance relationship (mean equals variance) is always shown to illustrate how the genewise variances may differ from a Poisson mean-variance relationship. Optionally, the raw variances and estimated tagwise variances can also be plotted. Estimated tagwise variances can be calculated using either qCML estimates of the tagwise dispersions (\code{estimateTagwiseDisp}) or Cox-Reid conditional inference estimates (\code{CRDisp}). A log-log scale is used for the plot.
}



\author{Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
plotMeanVar(d) # Produce a straight-forward mean-variance plot
meanvar <- plotMeanVar(d, show.raw.vars=TRUE) # Produce a mean-variance plot with the raw variances shown and save the means and variances for later use

## If we want to show estimated tagwise variances on the plot, we must first estimate them!
d <- estimateCommonDisp(d) # Obtain an estimate of the dispersion parameter
d <- estimateTagwiseDisp(d)  # Obtain tagwise dispersion estimates
plotMeanVar(d, meanvar=meanvar, show.tagwise.vars=TRUE, NBline=TRUE, dispersion.method="qcml") # Use previously saved object to speed up plotting; set dispersion.method to 'qcml' instead of default 'coxreid'
## We could also estimate common/tagwise dispersions using the Cox-Reid methods using CRDisp() with an appropriate design matrix
}

\seealso{
\code{\link{plotMDS.dge}}, \code{\link{plotSmear}} and
\code{\link{maPlot}} provide more ways of visualizing DGE data.

}

\keyword{algebra}




