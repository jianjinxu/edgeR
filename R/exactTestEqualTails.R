exactTestEqualTails <- function(y1,y2,dispersion)
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.
#	Left and right rejection regions have equal probability.

#	DOESN'T YET HANDLE DISPERSION=0 CASE

#	Gordon Smyth
#	28 Sep 2019.  Last modified 28 Sep 2011.
{
#	Convert matrices to vectors
	ntags <- NROW(y1)
	n1 <- NCOL(y1)
	n2 <- NCOL(y2)
	if(n1>1) y1 <- round(rowSums(y1))
	if(n2>1) y2 <- round(rowSums(y2))
	if(length(dispersion==1)) dispersion <- rep(dispersion,ntags)

#	Null fitted values
	y <- y1+y2
	mu <- y/(n1+n2)
	mu1 <- n1*mu
	mu2 <- n2*mu

	pvals <- rep(1,ntags)
	names(pvals) <- names(y1)

	p.bot <- dnbinom(y,size=(n1+n2)/dispersion,mu=y)
	size1 <- n1/dispersion
	size2 <- n2/dispersion
	left <- y1<mu1
	if(any(left)) {
		for (g in which(left)) {
			x <- 0:y1[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(y[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[left] <- pvals[left]/p.bot[left]
	}
	right <- y1>mu1
	if(any(right)) {
		for (g in which(right)) {
			x <- y1[g]:y[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(y[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[right] <- pvals[right]/p.bot[right]
	}
	pvals
}
