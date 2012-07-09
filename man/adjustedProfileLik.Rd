\name{adjustedProfileLik}
\alias{adjustedProfileLik}

\title{Compute Cox-Reid Adjusted Profile Likelihood for Negative Binomial GLMs}

\description{
Compute the Cox-Reid Adjusted Profile-likelihood for many negative binomial (NB) GLMs.
}

\usage{
adjustedProfileLik(dispersion, y, design, offset, adjust=TRUE)
}

\arguments{ 
\item{dispersion}{numeric scalar or vector giving the dispersion(s) towards which the tagwise dispersion parameters are shrunk.}

\item{y}{numeric matrix of counts}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{adjust}{logical, if \code{TRUE} then Cox-Reid adjustment is made to the log-likelihood, if \code{FALSE} then the log-likelihood is returned without adjustment. Default is \code{TRUE}.}

}

\value{
vector of adjusted profile log-likelihood values, one for each row of \code{y}.
}

\details{
In the \code{edgeR} context, \code{adjustedProfileLik} is a low-level function necessary for estimating dispersion parameters for NB GLMs. 
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\author{Yunshun Chen, Gordon Smyth}
\examples{
y <- matrix(rnbinom(1000, mu=10, size=2), ncol=4)
design <- matrix(1, 4, 1)
dispersion <- 0.5
apl <- adjustedProfileLik(dispersion, y, design, offset=0)
apl
}

\seealso{
\code{\link{dispCoxReidInterpolateTagwise}}, \code{\link{estimateGLMTagwiseDisp}}, \code{\link{maximizeInterpolant}}
}

\keyword{algebra}
