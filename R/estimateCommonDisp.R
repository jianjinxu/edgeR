estimateCommonDisp <- function(object,tol=1e-06,rowsum.filter=5,verbose=FALSE)
# Do two iterations of calculating pseudodata and estimating common dispersion, first one uses Poisson
# Davis McCarthy, Mark Robinson, Gordon Smyth.
# Created 2009. Last modified 26 March 2012.
{
	if(!is(object,"DGEList")) stop("Currently supports DGEList objects")
	group <- object$samples$group <- as.factor(object$samples$group)

	if( all(tabulate(group)<=1) ) {
		warning("There is no replication, setting dispersion to NA.")
		object$common.dispersion <- NA
		return(object)
	}

	tags.used <- rowSums(object$counts) > rowsum.filter
	pseudo.obj <- object[tags.used,]
	disp <- 0
	for(i in 1:2) {
		q2q.out <- equalizeLibSizes(object,disp=disp)
		pseudo.obj$counts <- q2q.out$pseudo[tags.used,,drop=FALSE]
		y <- splitIntoGroups(pseudo.obj)
		delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, der=0, doSum=FALSE)
		delta <- delta$maximum
		disp <- delta/(1-delta)
	}
	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	object$common.dispersion <- disp
	object$pseudo.alt <- q2q.out$pseudo
	object$conc <- q2q.out$conc
	object$common.lib.size <- q2q.out$N
	object
}

