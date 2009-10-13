approx.expected.info<-function(object,d,pseudo,robust=FALSE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to calculate the approximate expected information
{
	group<-object$samples$group
	levs.group<-levels(group)
	obs.inf<-rep(0,nrow(object$counts))
	for(i in 1:length(levs.group)) {
		if (sum( group==levs.group[i] ) > 1) {
			obs.inf<-obs.inf+condLogLikDerDelta(pseudo[,group==levs.group[i]],d,der=2,doSum=FALSE)*(-1)
		}
	}
	z<-rowSums(pseudo)
	if (robust) {
		require(MASS)
		inf.lm<-rlm(obs.inf~-1+z)
	} else {
		inf.lm<-lm(obs.inf~-1+z)
	}
	exp.inf.approx<-fitted(inf.lm)
	exp.inf.approx
}

