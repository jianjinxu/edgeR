#  SCORE.R

zscoreNBinom <- function(q, size, mu) 
#  Z-score equivalents for negative binomial deviates
#  Gordon Smyth, Aaron Lun
#  created 10 December 2011
#  last modified 9 March 2015
{
	z <- round(q)
	n <- length(q)
	size <- rep(size,length.out=n)
	mu <- rep(mu,length.out=n)
	d <- dnbinom(q,size=size,mu=mu, log=TRUE)
	up <- (q >= mu)

	if(any(up)) {
		p <- pnbinom(q[up],size=size[up],mu=mu[up],lower.tail=FALSE, log.p=TRUE)
		refined.p <- p + suppressWarnings(log1p(exp(d[up] - p)/2))
		z[up] <- suppressWarnings(qnorm(refined.p, lower.tail=FALSE, log.p=TRUE))
	}
	if(any(!up)) {
		p <- pnbinom(q[!up],size=size[!up],mu=mu[!up],lower.tail=TRUE, log.p=TRUE)
		refined.p <- p + suppressWarnings(log1p(-exp(d[!up] - p)/2))
		z[!up] <- suppressWarnings(qnorm(refined.p,lower.tail=TRUE, log.p=TRUE))
	}

	not.fin <- !is.finite(z)
	if (any(not.fin)) { 
		z[not.fin] <- (2L*up[not.fin] -1L) * sqrt(nbinomUnitDeviance(q[not.fin], 
			mean=mu[not.fin], dispersion=1/size[not.fin])) 
	}
	z
}
