equalizeLibSizes <- function(object, disp=0, N=exp(mean(log(object$samples$lib.size*object$samples$norm.factors))))
#	Normalize counts so that the library sizes can be treated as equal.
#	Uses a quantile-to-quantile transformation so that new count counts are equivalent deviates on the equalized scale.
#	Created by Davis McCarthy, July 2009. Last modified 19 June 2012.
{
	nrows<-nrow(object$counts)
	ncols<-ncol(object$counts)
	lib.size <- object$samples$lib.size * object$samples$norm.factors
	group<-as.factor(object$samples$group)
	levs.group<-levels(group)
	if(length(disp)==1 && disp==0) {
		maxr <- 1e06
		conc<-estimatePs(object,maxr)
	} else {
		zerodisp <- disp==0
		disp[zerodisp] <- 1e-06
		conc<-estimatePs(object,1/disp)
	}
	input.mean<-matrix(0,nrow=nrows,ncol=ncol(object$counts))
	output.mean<-matrix(0,nrow=nrows,ncol=ncol(object$counts))
	if(length(disp)==1) {
		disp <- matrix(disp,nrow=nrows,ncol=ncol(object$counts))
	}
	for(i in 1:length(levs.group)) {
		input.mean[,group==levs.group[i]]<-outer(conc$conc.group[,i],lib.size[group==levs.group[i]])
		output.mean[,group==levs.group[i]]<-outer(conc$conc.group[,i],rep(N,sum(group==levs.group[i])))
	}
	pseudo <- q2qnbinom(object$counts, input.mean=input.mean, output.mean=output.mean, dispersion=disp)
	pseudo[pseudo<0]<-0
	list(pseudo=pseudo,conc=conc,N=N)
}
