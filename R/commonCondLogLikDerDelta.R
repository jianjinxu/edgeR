commonCondLogLikDerDelta <- function(y, delta, der=0, doSum=FALSE) 
# Calculates the common conditional log-likelihood (i.e. summed over all tags) - necessary so that optimize can be applied in estimateCommonDisp
# Davis McCarthy, July 2009
{
	l0<-0
	for(i in 1:length(y)) {
		l0<-condLogLikDerDelta(y[[i]],delta,der=der,doSum=doSum)+l0
	}
	m0<-matrix(colSums(l0),nrow=1,ncol=1)
	m0
}