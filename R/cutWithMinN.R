cutWithMinN <- function(x, intervals=2, min.n=1)
#	Cut numeric x into intervals, as equally spaced as possible subject
#	to including a minimum number of values in each interval
#	Gordon Smyth
#	7 May 2011.  Last modified 9 May 2011.
{
#	Check input
	x <- as.numeric(x)
	isna <- is.na(x)
	if(any(isna)) {
		group <- rep.int(NA,length(x))
		out <- Recall(x=x[!isna],intervals=intervals,min.n=min.n)
		group[!isna] <- out$group
		out$group <- group
		return(out)
	}
	intervals <- as.integer(intervals)
	min.n <- as.integer(min.n)
	if(length(x) < intervals*min.n) stop("too few observations: length(x) < intervals*min.n")

#	if(length(unique(x)) < intervals*min.n) stop("too few unique values for x")

	breaks.eqx <- seq(from=min(x),to=max(x),length.out=intervals+1L)
	breaks.eqx[1] <- breaks.eqx[1]-1
	breaks.eqx[intervals+1L] <- breaks.eqx[intervals+1L]+1
	breaks.eqn <- quantile(x,p=seq(from=0,to=1,length.out=intervals+1L))
	breaks.eqn[1] <- breaks.eqn[1]-1
	breaks.eqn[intervals+1L] <- breaks.eqn[intervals+1L]+1

#	First try
	z <- cut(x,breaks=intervals,labels=FALSE)
	n <- tabulate(z)
	if(all(n>=min.n)) return(list(group=z,breaks=breaks.eqx))

	for (i in 1:10) {
		breaks <- (i*breaks.eqn+(10-i)*breaks.eqx)/10
		z <- cut(x,breaks=breaks,labels=FALSE)
		n <- tabulate(z)
		if(all(n>=min.n)) return(list(group=z,breaks=breaks))
	}

	stop("function has failed, perhaps because of too many tied values?")
}