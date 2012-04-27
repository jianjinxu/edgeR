\title{Plot Biological Coefficient of Variation}
\name{plotBCV}
\alias{plotBCV}
\description{
Plot genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
}
\usage{
plotBCV(object, xlab="logCPM", ylab="Biological coefficient of variation", pch=16, cex=0.2, ... )
}
\arguments{
  \item{object}{a \code{DGEList} object.}
  \item{xlab}{label for the x-axis.}
  \item{ylab}{label for the y-axis.}
  \item{pch}{the plotting symbol. See \code{\link{points}} for more details.}
  \item{cex}{plot symbol expansion factor. See \code{\link{points}} for more details.} 
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
The BCV is the square root of the negative binomial dispersion.
This function displays the common, trended and tagwise BCV estimates.
}

\value{
A plot is created on the current graphics device.
}

\author{Davis McCarthy, Yunshun Chen, Gordon Smyth}

\examples{
BCV.true <- 0.1
y <- DGEList(matrix(rnbinom(6000, size = 1/BCV.true^2, mu = 10),1000,6))
y <- estimateCommonDisp(y)
y <- estimateTrendedDisp(y)
y <- estimateTagwiseDisp(y)
plotBCV(y)
}

\keyword{plot}
