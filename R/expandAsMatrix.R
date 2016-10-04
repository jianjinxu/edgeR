expandAsMatrix <- function(x,dim=NULL,byrow=TRUE)
#	Convert scalar, row or column vector, or matrix, to be a matrix
#	Gordon Smyth, Yunshun Chen, Aaron Lun
#	26 Jan 2011.  Last modified 03 Oct 2016.
{
#	Check dim argument
	if(is.null(dim)) return(as.matrix(x))
	if(length(dim)<2) stop("dim must be numeric vector of length 2")
	dim <- round(dim[1:2])
	if(any(dim<0)) stop("negative dimensions not allowed")

#	x is a vector
	dx <- dim(x)
	if(is.null(dx)) {
		lx <- length(x)
		if(lx==1) return(matrix(x,dim[1],dim[2]))
		if(lx==dim[1] & lx==dim[2]) return(matrix(x,dim[1],dim[2],byrow=byrow))
		if(lx==dim[2]) return(matrix(x,dim[1],dim[2],byrow=TRUE))
		if(lx==dim[1]) return(matrix(x,dim[1],dim[2],byrow=FALSE))
		stop("x of unexpected length")
	}

#	x is a matrix or data.frame
	if(length(dx)<2) stop("x has less than 2 dimensions")
	if(length(dx)>2) stop("x has more than 2 dimensions")
	if(all(dx==dim)) return(as.matrix(x))

#	x is a compressedMatrix
	if(is(x, "compressedMatrix")) {
		if (attributes(x)$repeat.row && dim[2]==ncol(x)) {
			if (!byrow) warning("'byrow=FALSE' is not compatible with compressedMatrix settings")
			return(Recall(as.vector(x),dim=dim,byrow=TRUE))
		} else if (attributes(x)$repeat.col && dim[1]==nrow(x)) {
			if (byrow) warning("'byrow=TRUE' is not compatible with compressedMatrix settings")
			return(Recall(as.vector(x),dim=dim,byrow=FALSE))
		} else if (all(dx==1)) {
			return(Recall(as.vector(x),dim=dim))
		}
	}
	stop("x is matrix of wrong size")
}

