exactTestByDeviance <- function(y1,y2,dispersion=0,big.count=900)
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.

#	Rejection region is defined by large deviance statistics,
#  making this a conditional likelihood ratio test.

# Davis McCarthy, Gordon Smyth.
# Created 8 August 2011.  Last modified 2 October 2011.
{
	y1 <- as.matrix(y1)
	y2 <- as.matrix(y2)
	ntags <- nrow(y1)
	if(ntags!=nrow(y2)) stop("Number of rows of y1 not equal to number of rows of y2")
	if(any(is.na(y1)) || any(is.na(y2))) stop("NAs not allowed")
	n1 <- ncol(y1)
	n2 <- ncol(y2)

	if(n1==n2) return(exactTestDoubleTail(y1=y1,y2=y2,dispersion=dispersion,big.count=big.count))

	sum1 <- round(rowSums(y1))
	sum2 <- round(rowSums(y2))
	N <- sum1+sum2
	mu <- N/(n1+n2)
	if(all(dispersion==0)) return(binomTest(sum1,sum2,p=n1/(n1+n2)))
	if(any(dispersion==0)) stop("dispersion must be either all zero or all positive")
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)
	r <- 1/dispersion
	all.zeros <- N==0

	pvals <- rep(1,ntags)
	if(ntags==0) return(pvals)

#	Eliminate all zero rows
	if(any(all.zeros)) {
		pvals[!all.zeros] <- Recall(y1=y1[!all.zeros,,drop=FALSE],y2=y2[!all.zeros,,drop=FALSE],dispersion=dispersion[!all.zeros],big.count=big.count)
		return(pvals)
	}
	for (i in 1:ntags) {
		ind <- 0:N[i]
		p.top <- dnbinom(ind,size=n1*r[i],mu=n1*mu[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mu[i])
		p.obs <- dnbinom(sum1[i],size=n1*r[i],mu=n1*mu[i]) * dnbinom(sum2[i],size=n2*r[i],mu=n2*mu[i])
		## We need a robust way to choose the appropriate prob masses to sum over
		## We use the LR stat to choose---choose all values of support with LR stat >= obs LR stat
		## Then sum prob masses over these chosen values of support
		ind.mat <- cbind(ind, N[i]-ind)
		rmat.i <- matrix(cbind(n1*r[i], n2*r[i]), nrow=length(ind), ncol=2, byrow=TRUE)
		mumat.i <- matrix(cbind(n1*mu[i], n2*mu[i]), nrow=length(ind), ncol=2, byrow=TRUE)
		dev.obs <- .nbdev(cbind(sum1[i],sum2[i]), rmat.i[1,,drop=FALSE], mumat.i[1,,drop=FALSE])
		dev.all <- .nbdev(ind.mat, rmat.i, mumat.i)
		keep <- dev.all>=dev.obs
		p.bot <- dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mu[i])
		pvals[i] <- sum(p.top[keep]/p.bot)
	}
	pmin(pvals,1)
}

.nbdev <- function(y, r, mu)
	## Compute NB deviances for observations, with potentially different size (r = 1/dispersion) and mean (mu)
	## Arguments can be vectors or matrices
	## Created by Davis McCarthy, 5 August 2011.
	## Last modified 8 August 2011.
{
	if(!identical(dim(y),dim(r)) | !identical(dim(y),dim(mu)))
		stop("Matrices y, r and m are not all the same size.\n")
	y[y==0] <- 1e-08
	dev <- 2*rowSums( y*log( y / mu ) - (y + r)*log( (1+y/r) / (1+mu/r) ) )
	dev
}
