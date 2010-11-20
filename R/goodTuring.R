goodTuring <- function(x,plot=FALSE)
#	Simple Good-turing algorithm for frequency estimation
#	as described by Gale and Sampson (1995)

#	Has been tested against Sampson's C code from
#	http://www.grsampson.net/RGoodTur.html
#	and gives identical results.

#	Gordon Smyth
#	9 Nov 2010.  Last modified 21 Nov 2010.
{
#	Raw frequencies
	n <- tabulate(x+1L)
	n0 <- n[1]
	n <- n[-1]
	max.x <- length(n)
	r <- seq.int(from=1L,to=max.x)

#	Fit a linear trend to log-frequencies
	n.pos <- n[n>0]
	r.pos <- r[n>0]
	l <- length(n.pos)
	q <- diff(c(0L,r.pos,2L*r.pos[l]-r.pos[l-1]),lag=2)/2
	z <- n.pos/q
	design <- cbind(1,log(r.pos))
	fit <- lm.fit(x=design,y=log(z))
	if(plot) {
		plot(log(r.pos),log(z),xlab="log count",ylab="log frequency")
		abline(fit)
	}

#	Smoothed counts and ratios
	r.long <- c(r,max.x+1L)
	n.long.smooth <- exp(fit$coef[1]+fit$coef[2]*log(r.long))
	ratio.smooth <- n.long.smooth[r+1L]/n.long.smooth[r]

#	Empirical ratios
	n.long <- c(n,n[max.x])
	ratio <- n.long[r+1L]/n

#	Posterior expectations
	r.post.smooth <- (r+1L)*ratio.smooth
	r.post <- (r+1L)*ratio

#	Combine empirical and smoothed expectations
	se <- (r+1)*sqrt(ratio/n*(1+ratio))
	z.stat <- (r.post-r.post.smooth)/se
	first.r.equivalent <- which(abs(z.stat)<1.96)[1]
	if(first.r.equivalent>1) r.post.smooth[1:(first.r.equivalent-1)] <- r.post[1:(first.r.equivalent-1)]
	r.pos.post <- r.post.smooth[n>0]

#	Estimated frequencies
	N <- sum(n.pos*r.pos)
	P0 <- n[1]/N
	N.post <- sum(n.pos*r.pos.post)
	list(count=r.pos,estimated.p=(1-P0)*r.pos.post/N.post,P0=P0,n0=n0)
}
