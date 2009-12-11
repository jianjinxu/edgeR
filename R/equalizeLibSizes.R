equalizeLibSizes <- function(object, disp=0, N=prod(object$samples$lib.size)^(1/ncol(object$counts)), null.hypothesis=FALSE)
# Davis McCarthy, July 2009
# A function that simply adjusts the counts for library size for a fixed value of the dispersion parameter
{
	nrows<-nrow(object$counts)
	ncols<-ncol(object$counts)
	lib.size<-object$samples$lib.size
	group<-as.factor(object$samples$group)
	levs.group<-levels(group)
	y<-splitIntoGroups(object)
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
	if (null.hypothesis) {
			input.mean<-outer(conc$conc.common,lib.size)
			output.mean<-outer(conc$conc.common,rep(N,ncols))
	} else {
		for(i in 1:length(levs.group)) {
			input.mean[,group==levs.group[i]]<-outer(conc$conc.group[,i],lib.size[group==levs.group[i]])
			output.mean[,group==levs.group[i]]<-outer(conc$conc.group[,i],rep(N,sum(group==levs.group[i])))
		}
	}
	pseudo <- q2qnbinom(object$counts, input.mean=input.mean, output.mean=output.mean, dispersion=disp)
	pseudo[pseudo<0]<-0 
	return(list(pseudo=pseudo,conc=conc,N=N))
}

