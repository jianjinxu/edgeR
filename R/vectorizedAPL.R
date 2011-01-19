
adjustedProfileLik <- function(dispersion, y, design, offset)
## tagwise Cox-Reid adjusted profile likelihoods for the dispersion
## dispersion can be scalar or tagwise vector
## y is matrix: rows are genes/tags/transcripts, columns are samples/libraries
## offset is matrix of the same dimensions as y
## Yunshun Chen, Gordon Smyth
## Created June 2010. Last modified 19 Jan 2011.
{
	if(any(dim(y)!=dim(offset))) stop("offset must be a matrix of same dimensions as y, the matrix of counts.")
	ntags <- nrow(y)
	nlibs <- ncol(y)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

#	Fit tagwise linear models
	ls <- mglmLS(y, design, dispersion, offset = offset)

#	Compute log-likelihood
	mu <- ls$fitted
	if(dispersion[1] == 0){
		loglik <- rowSums(dpois(y,lambda=mu,log = TRUE))
	} else {
		loglik <- rowSums(dnbinom(y,size=1/dispersion,mu=mu,log = TRUE))
	}

#	Cox-Reid adjustment
	A <- .vectorizedXWX(design, mu, dispersion)
	D <- .vectorizedLDL(A)
	cr <- 0.5*rowSums(log(abs(D)))

	loglik - cr
}


.vectorizedXWX <- function(design, fitted, dispersion)
#	XWX in packed form, tagwise in rows of matrix
{
	ntags <- nrow(fitted)
	nlibs <- ncol(fitted)
	ncoef <- ncol(design)
	W <- fitted/(1+dispersion*fitted)
	A <- matrix(0,ntags,ncoef*(ncoef+1)/2)
	colstart <- 0
	for (i in 1:ncoef) {
		WX <- t(design[,i] * t(W))
		colstart <- colstart
		A[,colstart+(1:i)] <- WX %*% design[,1:i]
		colstart <- colstart+i
	}
	A
}

.vectorizedLDL <- function(A)
## LDL decomposition of XWX in packed form
## Each row of A represents a XWX matrix
## Only the elements of D are returned in corresponding rows
{
	p <- floor(sqrt(2*ncol(A)))
	cum <- c(0,0,cumsum(1:p))
	if(cum[p+2] != ncol(A))
		stop("Dimension doesn't match!")
	index.l <- (1:cum[p+2])[-cum]
	d <- matrix(0, nrow(A), p)
	l <- matrix(0, nrow(A), cum[p+1])
	for(j in 1:(p-1)){		
		d[,j] <- A[,cum[j+2]] - rowSums(as.matrix(l[,(cum[j]+1):(cum[j]+j-1)]^2 * d[,1:(j-1)]))
		for(i in (j+1):p){
			if(j == 1){
				l[,cum[i]+j] <- A[,cum[i+1]+j]/d[,j]
			} else {
				l[,cum[i]+j] <- (A[,cum[i+1]+j] - rowSums(as.matrix(l[,(cum[i]+1):(cum[i]+j-1)]*l[,(cum[j]+1):(cum[j]+j-1)]*d[,1:(j-1)])))/d[,j]
			}
		}
	}
	d[,p] <- A[,cum[p+2]] - rowSums(as.matrix(l[,(cum[p]+1):(cum[p]+p-1)]^2 * d[,1:(p-1)]))
	return(d)
}
