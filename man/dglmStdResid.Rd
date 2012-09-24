\name{dglmStdResid}
\alias{dglmStdResid}
\alias{getDispersions}

\title{Visualize the mean-variance relationship in DGE data using standardized residuals}

\description{Appropriate modelling of the mean-variance relationship in DGE data is important for making inferences about differential expression. However, the standard approach to visualizing the mean-variance relationship is not appropriate for general, complicated experimental designs that require generalized linear models (GLMs) for analysis. Here are functions to compute standardized residuals from a Poisson GLM and plot them for bins based on overall expression level of tags as a way to visualize the mean-variance relationship. A rough estimate of the dispersion parameter can also be obtained from the standardized residuals.}


\usage{
dglmStdResid(y, design, dispersion=0, offset=0, nbins=100, make.plot=TRUE,
          xlab="Mean", ylab="Ave. binned standardized residual", ...)
getDispersions(binned.object)
}

\arguments{ 
\item{y}{numeric matrix of counts, each row represents one tag, each column represents one DGE library.}

\item{design}{numeric matrix giving the design matrix of the GLM. Assumed to be full column rank.}

\item{dispersion}{numeric scalar or vector giving the dispersion parameter for each GLM. Can be a scalar giving one value for all tags, or a vector of length equal to the number of tags giving tag-wise dispersions.}

\item{offset}{numeric vector or matrix giving the offset that is to be included in teh log-linear model predictor. Can be a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}

\item{nbins}{scalar giving the number of bins (formed by using the quantiles of the genewise mean expression levels) for  which to compute average means and variances for exploring the mean-variance relationship. Default is \code{100} bins}

\item{make.plot}{logical, whether or not to plot the mean standardized residual for binned data (binned on expression level). Provides a visualization of the mean-variance relationship. Default is \code{TRUE}.}

\item{xlab}{character string giving the label for the x-axis. Standard graphical parameter. If left as the default, then the x-axis label will be set to "Mean".}

\item{ylab}{character string giving the label for the y-axis. Standard graphical parameter. If left as the default, then the y-axis label will be set to "Ave. binned standardized residual".}

\item{\dots}{further arguments passed on to \code{plot}}

\item{binned.object}{list object, which is the output of \code{dglmStdResid}.}

}

\value{ \code{dglmStdResid} produces a mean-variance plot based on standardized residuals from a Poisson model fitfor each tag for the DGE data. \code{dglmStdResid} returns a list with the following elements:
	\item{ave.means}{vector of the average expression level within each bin of observations}
	\item{ave.std.resid}{vector of the average standardized Poisson residual within each bin of tags}
	\item{bin.means}{list containing the average (mean) expression level (given by the fitted value from the given Poisson model) for observations divided into bins based on amount of expression}
	\item{bin.std.resid}{list containing the standardized residual from the given Poisson model for observations divided into bins based on amount of expression}
	\item{means}{vector giving the fitted value for each observed count}
	\item{standardized.residuals}{vector giving approximate standardized residual for each observed count}
	\item{bins}{list containing the indices for the observations, assigning them to bins}
	\item{nbins}{scalar giving the number of bins used to split up the observed counts}
	\item{ngenes}{scalar giving the number of genes/tags in the dataset}
	\item{nlibs}{scalar giving the number of libraries in the dataset}

\code{getDispersions} computes the dispersion from the standardized residuals and returns a list with the following components:
	\item{bin.dispersion}{vector giving the estimated dispersion value for each bin of observed counts, computed using the average standardized residual for the bin}
	\item{bin.dispersion.used}{vector giving the actual estimated dispersion value to be used. Some computed dispersions using the method in this function can be negative, which is not allowed. We use the dispersion value from the nearest bin of higher expression level with positive dispersion value in place of any negative dispersions.}
	\item{dispersion}{vector giving the estimated dispersion for each observation, using the binned dispersion estimates from above, so that all of the observations in a given bin get the same dispersion value.}

}

\details{
This function is useful for exploring the mean-variance relationship in the data. Raw or pooled variances cannot be used for complex experimental designs, so instead we can fit a Poisson model using the appropriate design matrix to each tag and use the standardized residuals in place of the pooled variance (as in \code{plotMeanVar}) to visualize the mean-variance relationship in the data. The function will plot the average standardized residual for observations split into \code{nbins} bins by overall expression level. This provides a useful summary of how the variance of the counts change with respect to average expression level (abundance). A line showing the Poisson mean-variance relationship (mean equals variance) is always shown to illustrate how the genewise variances may differ from a Poisson mean-variance relationship. A log-log scale is used for the plot.

The function \code{mglmLS} is used to fit the Poisson models to the data. This code is fast for fitting models, but does not compute the value for the leverage, technically required to compute the standardized residuals. Here, we approximate the standardized residuals by replacing the usual denominator of \code{ ( 1 - leverage )} by \code{ ( 1 - p/n ) }, where n is the number of observations per tag (i.e. number of libraries) and p is the number of parameters in the model (i.e. number of columns in the full-rank design matrix.

}



\author{Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=2),ncol=4)
design <- model.matrix(~c(0,0,1,1)+c(0,1,0,1))
binned <- dglmStdResid(y, design, dispersion=0.5)

getDispersions(binned)$bin.dispersion.used # Look at the estimated dispersions for the bins

}

\seealso{
\code{\link{plotMeanVar}}, \code{\link{plotMDS.DGEList}}, \code{\link{plotSmear}} and \code{\link{maPlot}} provide more ways of visualizing DGE data.

}

\keyword{algebra}




