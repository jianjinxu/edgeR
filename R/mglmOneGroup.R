mglmOneGroup <- function(y,dispersion=0,offset=0,maxit=50,trace=FALSE,tol=1e-6)
#	Fit null (single-group) negative-binomial glm with log-link to DGE data
#  Gordon Smyth
#	18 Aug 2010. Last modified 25 July 2012.
{
#	Check input values for y
	y <- as.matrix(y)
	if(any(y<0)) stop("y must be non-negative")
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Treat all zero rows as special case
	beta <- rep(-Inf,ntags)
	names(beta) <- rownames(y)

#	Check input values for dispersion
	if(any(dispersion<0)) stop("dispersion must be non-negative")

#	Poisson special case
	N <- exp(offset)
	if(all(dispersion==0)) {
		if(is.null(dim(N)))
			m <- mean(N)
		else
			m <- .rowMeans(N,ntags,nlibs)
		return(log(.rowMeans(y/m,ntags,nlibs)))
	}
	dispersion <- rep(dispersion,length=ntags)

#	Check input values for offset
	offset <- expandAsMatrix(offset,dim(y))
	N <- expandAsMatrix(N,dim(y))

#	Exact solution for gamma limit
	beta <- log(.rowMeans(y/N,ntags,nlibs))

#	Single library as special case
	if(nlibs==1) return(beta)

#	Fisher scoring iteration	
	iter <- 0
	i <- is.finite(beta)
	while(any(i)) {
		iter <- iter+1
		if(iter>maxit) {
			warning("max iterations exceeded")
			return(beta)
		}
		if(trace) cat("Iter=",iter,"Still converging=",sum(i),"\n")
		mu <- exp(beta[i]+offset[i,,drop=FALSE])
		var.div.mu <- 1+dispersion[i]*mu
		m <- nrow(mu)
		dl <- .rowSums((y[i,,drop=FALSE]-mu)/var.div.mu,m,nlibs)
		info <- .rowSums(mu/var.div.mu,m,nlibs)
		step <- dl/info
		beta[i] <- beta[i]+step
		i[i] <- abs(step)>tol
	}
	beta
}
