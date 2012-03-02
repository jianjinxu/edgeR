\title{Plots biological coefficient of variation against abundance (counts per million)}
\name{plotBCV}
\alias{plotBCV}
\description{
Calculate genewise biological coefficient of variation (BCV) and gene abundance (using counts per million) from a \code{DGEList} object, then plot BCV against abundance.
}
\usage{
plotBCV(object, xlab="Abundance (log2 counts per million)", ylab="Biological coefficient of variation (BCV)", ... )
}
\arguments{
  \item{object}{a \code{DGEList} object.}
  \item{xlab}{title for the x-axis}
  \item{ylab}{title for the y-axis}
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
This function provides a way to examine the relationship between biological coefficient of variation (BCV) and gene (tag) abundance visually. Under the Negative Binomial model, the square root of dispersion(s), which is provided from the \code{DGEList} object, can be considered as the BCV. The counts per million is calculated using the abundance from the \code{DGEList} object, and it is on log2 scale. 
}

\value{
A plot is created on the current graphics device.
}

\author{Davis McCarthy and Yunshun Chen}

\examples{
y <- DGEList(matrix(rnbinom(6000, size = 1/2, mu = 10),1000,6))
y <- estimateCommonDisp(y)
y <- estimateTrendedDisp(y)
y <- estimateTagwiseDisp(y)
plotBCV(y)
}

\keyword{plot}
