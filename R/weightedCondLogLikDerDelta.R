weightedCondLogLikDerDelta <- function(y, delta, tag, prior.n=10, ntags=nrow(y[[1]]), der=0)
# Calculates weighted conditional log-likelihood for a tag - necessary to estimate tagwise dispersions
{
	l0<-rep(0,ntags)
	onev<-rep(1,ntags)
	for(i in seq_len(length(y))) {
		l0<-condLogLikDerDelta(y[[i]],delta,der=der)+l0
	}
	m0<-sum(l0)
	l0a<-l0 + (prior.n/ntags)*m0
	l0a[tag]
}