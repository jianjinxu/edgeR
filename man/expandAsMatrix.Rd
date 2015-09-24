\title{expandAsMatrix}
\name{expandAsMatrix}
\alias{expandAsMatrix}
\description{
Expand scalar or vector to a matrix.
}
\usage{
expandAsMatrix(x, dim=NULL, byrow=TRUE)
}
\arguments{
  \item{x}{scalar, vector or matrix. If a vector, length must match one of the output dimensions.}
  \item{dim}{required dimension for the output matrix.}
  \item{byrow}{logical. Should the matrix be filled by columns or by rows (the default) if the length of the input vector is equal to both dimensions?}
}

\details{
This function expands a row or column vector to be a matrix.
It is used internally in edgeR to convert offsets to a matrix.
}

\value{
Numeric matrix of dimension \code{dim}.
}
\author{Gordon Smyth}
\examples{
expandAsMatrix(1:3,c(4,3))
expandAsMatrix(1:4,c(4,3))
}
\keyword{hplot}
