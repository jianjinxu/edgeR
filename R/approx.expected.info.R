# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to calculate the approximate expected information
approx.expected.info<-function(object,d,qA,robust=FALSE) {
	group<-object$group
	levs.group<-levels(group)
	obs.inf<-rep(0,nrow(object$data))
	for(i in 1:length(levs.group)) {
		if (sum( group==levs.group[i] ) > 1) {
			obs.inf<-obs.inf+condLogLikDerDelta(qA$pseudo[,group==levs.group[i]],d,der=2,doSum=FALSE)*(-1)
		}
	}
	t<-rowSums(qA$pseudo)
	if (robust) {
		require(MASS)
		inf.lm<-rlm(obs.inf~-1+t)
	} else {
		inf.lm<-lm(obs.inf~-1+t)
	}
	exp.inf.approx<-fitted(inf.lm)
	exp.inf.approx
}

