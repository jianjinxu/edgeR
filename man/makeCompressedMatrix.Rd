\title{makeCompressedMatrix}
\name{makeCompressedMatrix}
\alias{makeCompressedMatrix}
\alias{compressedMatrix}
\alias{[.compressedMatrix}
\alias{as.matrix.compressedMatrix}

\description{
Construct a compressedMatrix object from a scalar, vector or matrix.
}

\usage{
makeCompressedMatrix(x, byrow=TRUE)

\method{[}{compressedMatrix}(x, i, j, ...)
\method{as.matrix}{compressedMatrix}(x, ...)
}

\arguments{
  \item{x}{For \code{makeCompressedMatrix}, a scalar, vector, matrix or compressedMatrix object.
        For the S3 methods, a compressedMatrix object.}
  \item{byrow}{logical. If \code{x} is a vector, should it be repeated across rows (default) or across columns?}
  \item{i, j}{subset indices to apply to \code{x}.}
  \item{...}{additional arguments, ignored.}
}

\details{
This function creates a compressedMatrix object from \code{x}. 
The compressedMatrix class inherits from a matrix and holds two logical scalar attributes \code{repeat.row} and \code{repeat.col}.
Each attribute specifies whether the values are to be repeated across rows and/or across columns.
This avoids the need to store redundant values in a full-sized matrix of dimensions \code{dim}, as would be done with \code{\link{expandAsMatrix}}.

To illustrate, consider that rows usually correspond to genes while columns usually correspond to libraries.
If we have a vector of library sizes, this will hold one unique value per library that is the same for all genes.
Thus, we should use \code{byrow=TRUE}, which will construct a compressedMatrix object storing one row containing this vector.
Here, \code{repeat.row=TRUE} and \code{repeat.col=FALSE}, indicating that the row is to be repeated for all genes.

On the other hand, we may have a vector of gene-specific values that is the same for all libraries (e.g., dispersions).
In this case, we should use \code{byrow=FALSE} to construct the compressedMatrix object.
This will store one column with \code{repeat.row=FALSE} and \code{repeat.col=TRUE}, indicating that the column should be repeated across libraries.

In cases where \code{x} is a scalar, \code{byrow} is ignored and both \code{repeat.row} and \code{repeat.col} will be \code{TRUE} by default.
If \code{x} is a matrix, both attributes will be \code{FALSE}.
If \code{x} is a compressedMatrix, it will be returned without modification.

Subsetting of a compressedMatrix object depends on the values of \code{repeat.row} and \code{repeat.col}.
If the rows are repeated, any subsetting by row will be ignored.
Similarly, if the columns are repeated, any subsetting by column will be ignored.
This reflects the fact that the repeated dimension has no fixed size, so subsetting on it is meaningless.
If neither are repeated, subsetting behaves as it would for a normal matrix.

Calling \code{as.matrix} will return the raw matrix without attributes or classes.
If either the columns or rows are repeated, the corresponding dimension in the returned matrix will be of length 1.
Otherwise, it will be of arbitrary length depending on the size/length of \code{x} used originally to construct \code{y}.
A compressedMatrix object can also be used as input to \code{\link{expandAsMatrix}}, which will expand it to the specified dimensions.

The compressedMatrix is used throughout edgeR to save space in storing offsets and (to a lesser extent) weights.
This is because, for routine analyses, offsets are the same for all genes so it makes little sense to expand it to the full dimensions of the count matrix.
Most functions will accept a compressedMatrix as input to \code{offset} or \code{weights} arguments.
}

\value{
A object of class compressedMatrix, containing \code{x} and the additional attributes \code{repeat.row} and \code{repeat.col}.
}

\author{Aaron Lun}

\examples{
# Repeated rows:
library.sizes <- runif(4, 1e6, 2e6)
lib.mat <- makeCompressedMatrix(library.sizes, byrow=TRUE)
lib.mat

lib.mat[,1:2] # subset by column works as expected
lib.mat[1:10,] # subset by row has no effect (see Details)
as.matrix(lib.mat)
expandAsMatrix(lib.mat, dim=c(10, 4))

# Repeated columns:
gene.disp <- runif(10, 0.01, 0.1)
disp.mat <- makeCompressedMatrix(gene.disp, byrow=FALSE)
disp.mat

disp.mat[,1:2] # subset by column has no effect
disp.mat[1:5,] # subset by row works as expected
as.matrix(disp.mat)
expandAsMatrix(disp.mat, dim=c(10, 4), byrow=FALSE)

# Scalar:
weights <- makeCompressedMatrix(1)
weights[1:10,] # subsetting has no effect
weights[,1:10]
as.matrix(weights)
expandAsMatrix(weights, dim=c(10, 4))

# Matrix:
offsets <- makeCompressedMatrix(matrix(runif(40), 10, 4))
offsets[1:5,]
offsets[,1:2]
as.matrix(offsets)
expandAsMatrix(offsets, dim=c(10, 4))
}

\seealso{
\code{\link{expandAsMatrix}}
}

\keyword{hplot}
