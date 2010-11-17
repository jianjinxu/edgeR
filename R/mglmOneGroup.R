mglmOneGroup <- function(y,offset=0,dispersion=0,maxit=50,verbose=FALSE)
#	Fit null (single-group) negative-binomial glm with log-link to DGE data
#	18 Aug 2010. Last modified 19 Aug 2010.
{
	if(any(y<0)) stop("y must be non-negative")
	y <- as.matrix(y)
	if(any(dispersion<0)) stop("dispersion must be non-negative")
	ntags <- nrow(y)
	nlib <- ncol(y)
	offset <- matrix(offset,ntags,nlib,byrow=TRUE)
	dispersion <- rep(dispersion,length=ntags)

	beta <- rep(-Inf,nrow(y))
	names(beta) <- rownames(y)
	i <- rowSums(y)>0
	if(!any(i)) return(beta)

#	Starting values
	betaold <- log(rowMeans(y[i,,drop=FALSE]/exp(offset[i,,drop=FALSE])))
	mu <- exp(betaold+offset[i,,drop=FALSE])

#	Fisher scoring iteration	
	iter <- 0
	while(any(i)) {
		iter <- iter+1
		if(iter>maxit) {
			warning("max iterations exceeded")
			return(beta)
		}
		dl <- rowSums((y[i,,drop=FALSE]-mu)/(1+dispersion[i]*mu))
		info <- rowSums(mu/(1+dispersion[i]*mu))
		step <- dl/info
		beta[i] <- betaold+step
		i[i] <- abs(step)>1e-6
		if(verbose) cat("Iter=",iter,"Still converging=",sum(i),"\n")
		mu <- exp(beta[i]+offset[i,,drop=FALSE])
		betaold <- beta[i]
	}
	beta
}
