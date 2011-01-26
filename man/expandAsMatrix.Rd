\title{expandAsMatrix}
\name{expandAsMatrix}
\alias{expandAsMatrix}
\description{
Expand scalar or vector to a matrix.
}
\usage{
expandAsMatrix(x, dim)
}
\arguments{
  \item{x}{scalar, vector or matrix. If a vector, length must match one of the output dimensions.}
  \item{dim}{required dimension for the output matrix.}
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
\seealso{
\code{\link{mglmLS}}.
}
\keyword{hplot}
