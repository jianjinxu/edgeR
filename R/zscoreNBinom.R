#  SCORE.R

zscoreNBinom <- function(q, size, mu) 
#  Z-score equivalents for negative binomial deviates
#  Gordon Smyth
#  10 December 2011
{
	z <- round(q)
	n <- length(q)
	size <- rep(size,length.out=n)
	mu <- rep(mu,length.out=n)
	d <- dnbinom(q,size=size,mu=mu)
	up <- (q >= mu)
	if(any(up)) {
		p <- pnbinom(q[up],size=size[up],mu=mu[up],lower.tail=FALSE)
		z[up] <- qnorm(p+d[up]/2,lower.tail=FALSE)
	}
	if(any(!up)) {
		p <- pnbinom(q[!up],size=size[!up],mu=mu[!up],lower.tail=TRUE)
		z[!up] <- qnorm(p-d[!up]/2,lower.tail=TRUE)
	}
	z
}
