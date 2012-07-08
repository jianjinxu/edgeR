\name{gof}
\alias{gof}

\title{Goodness of Fit Tests for Multiple GLM Fits}

\description{Conducts deviance goodness of fit tests for each fit in a \code{DGEGLM} object}

\usage{
gof(glmfit, pcutoff=0.1, adjust="holm", plot=FALSE, main="qq-plot of genewise goodness of fit", ...)
}
\arguments{ 

\item{glmfit}{\code{DGEGLM} object containing results from fitting NB GLMs to genes in a DGE dataset. Output from \code{glmFit}.}

\item{pcutoff}{scalar giving the cut-off value for the Holm-adjusted p-value. Genes with Holm-adjusted p-values lower than this cutoff value are flagged as `dispersion outlier' genes.}

\item{adjust}{method used to adjust goodness of fit p-values for multiple testing.}
\item{plot}{logical, if \code{TRUE} a qq-plot is produced.}
\item{main}{character, title for the plot.}
\item{\dots}{other arguments are passed to \code{qqnorm}.}
}

\details{
If \code{plot=TRUE}, produces a plot similar to Figure 2 of McCarthy et al (2012).
}

\value{
This function returns a list with the following components:
\item{gof.statistics}{numeric vector of deviance statistics, which are the statistics used for the goodness of fit test}
\item{gof.pvalues}{numeric vector of p-values providing evidence of poor fit; computed from the chi-square distribution on the residual degrees of freedom from the GLM fits.}
\item{outlier}{logical vector indicating whether or not each gene is a `dispersion outlier' (i.e., the model fit is poor for that gene indicating that the dispersion estimate is not good for that gene).} 
\item{df}{scalar, the residual degrees of freedom from the GLM fit for which the goodness of fit statistics have been computed. Also the degrees of freedom for the goodness of fit statistics for the LR (chi-quare) test for significance.}
}

\author{Davis McCarthy and Gordon Smyth}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\examples{
nlibs <- 3
ntags <- 100
dispersion.true <- 0.1

# Make first transcript respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ntags-1)))
mu.true <- 2^(beta.true \%*\% t(design))

# Generate count data
y <- rnbinom(ntags*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ntags,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("Gene",1:ntags,sep="")
d <- DGEList(y)

# Normalize
d <- calcNormFactors(d)

# Fit the NB GLMs
fit <- glmFit(d, design, dispersion=dispersion.true)
# Check how good the fit is for each gene
gof(fit)
}

\seealso{
\code{\link{qqnorm}}.

\code{\link{glmFit}} for more information on fitting NB GLMs to DGE data.
}

\keyword{algebra}
