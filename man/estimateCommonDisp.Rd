\name{estimateCommonDisp}
\alias{estimateCommonDisp}

\title{Estimates the Negative Binomial Common Dispersion by Maximizing the Negative Binomial Conditional Common Likelihood}

\description{Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags for the unadjusted counts provided. }


\usage{
estimateCommonDisp(object, tol=1e-06, rowsum.filter=5)
}

\arguments{ 

\item{object}{\code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts), and \code{samples} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{tol}{numeric scalar providing the tolerance to be passed to \code{optimize}; default value is \code{1e-06}}

\item{rowsum.filter}{numeric scalar giving a value for the filtering out of low abundance tags in the estimation of the common dispersion. Only tags with total sum of counts above this value are used in the estimation of the common dispersion. Low abundance tags can adversely affect the estimation of the common dispersion, so this argument allows the user to select an appropriate filter threshold for the tag abundance.}
}

\value{ \code{estimateCommonDisp} produces an object of class \code{DGEList} with the following components.
	\item{common.dispersion}{estimate of the common dispersion; the value for \code{phi}, the dispersion parameter in the NB model, that maximizes the negative binomial common likelihood on the \code{phi} scale}
	\item{counts}{table of unadjusted counts}
	\item{group}{vector indicating the group to which each library belongs}
	\item{lib.size}{vector containing the unadjusted size of each library}
	\item{pseudo.alt}{table of adjusted counts; quantile-to-quantile method (see \code{q2qnbinom}) used to adjust the raw counts so that library sizes are equal; adjustment here done under the alternative hypothesis that there is a true difference between groups}
	\item{conc}{list containing the estimates of the concentration of each tag in the underlying sample; \code{conc$p.common} gives estimates under the null hypothesis of no difference between groups; \code{conc$p.group} gives the estimate of the concentration for each tag within each group; concentration is a measure of abundance and thus expression level for the tags}
	\item{common.lib.size}{the common library size to which the count libraries have been adjusted}
}

\details{
The method of conditional maximum likelihood assumes that library sizes are equal, which is not true in general, so pseudocounts (counts adjusted so that the library sizes are equal) need to be calculated. The function \code{equalizeLibSizes} is called to adjust the counts using a quantile-to-quantile method, but this requires a fixed value for the common dispersion parameter. To obtain a good estimate for the common dispersion, pseudocounts are calculated under the Poisson model (dispersion is zero) and these pseudocounts are used to give an estimate of the common dispersion. This estimate of the common dispersion is then used to recalculate the pseudocounts, which are used to provide a final estimate of the common dispersion.
}


\references{
Robinson MD and Smyth GK (2008). Small-sample estimation of negative
binomial dispersion, with applications to SAGE data. \emph{Biostatistics},
9, 321-332
}

\author{Mark Robinson, Davis McCarthy}
\examples{
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d<-DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
cmdisp<-estimateCommonDisp(d)
}

\seealso{
\code{\link{estimateTagwiseDisp}} can be used to estimate a value for the dispersion parameter for each tag/transcript. The estimates are stabilized by squeezing the estimates towards the common value calculated by \code{estimateCommonDisp}.
}

\keyword{algebra}
