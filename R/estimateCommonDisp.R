estimateCommonDisp <- function(object,tol=1e-06,rowsum.filter=5,verbose=FALSE)
#	Estimate common dispersion using exact conditional likelihood
#	Davis McCarthy, Mark Robinson, Gordon Smyth.
#	Created 2009. Last modified 20 Nov 2013.
{
	if(!is(object,"DGEList")) stop("Currently supports DGEList objects")
	object <- validDGEList(object)
	group <- object$samples$group <- as.factor(object$samples$group)

	if( all(tabulate(group)<=1) ) {
		warning("There is no replication, setting dispersion to NA.")
		object$common.dispersion <- NA
		return(object)
	}

	tags.used <- rowSums(object$counts) > rowsum.filter
	pseudo.obj <- object[tags.used,]

#	Start from small dispersion
	disp <- 0.01
	for(i in 1:2) {
		out <- equalizeLibSizes(object,dispersion=disp)
		pseudo.obj$counts <- out$pseudo.counts[tags.used,,drop=FALSE]
		y <- splitIntoGroups(pseudo.obj)
		delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, der=0)
		delta <- delta$maximum
		disp <- delta/(1-delta)
	}
	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	object$common.dispersion <- disp
	object$pseudo.counts <- out$pseudo.counts
	object$pseudo.lib.size <- out$common.lib.size

#	Average logCPM
#	Note different behaviour to estimateGLMCommonDisp:
#	Here the estimated common.dispersion is used to compute AveLogCPM whereas
#	estimateGLMCommonDisp calculates AveLogCPM prior to estimating the common
#	dispersion using a pre-set dispersion of 0.05
	object$AveLogCPM <- aveLogCPM(object)

	object
}

