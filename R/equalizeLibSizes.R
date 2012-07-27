equalizeLibSizes <- function(object, dispersion=0, common.lib.size=NULL)
#	Normalize counts so that the library sizes can be treated as equal.
#	Uses a quantile-to-quantile transformation so that new count counts are equivalent deviates on the equalized scale.
#	Davis McCarthy, Gordon Smyth.
#	Created July 2009. Last modified 25 July 2012.
{
	d <- dim(object)
	ntags <- d[1]
	nlibs <- d[2]
	lib.size <- object$samples$lib.size * object$samples$norm.factors
	if(is.null(common.lib.size)) common.lib.size <- exp(mean(log(lib.size)))
	levs.group <- unique(object$samples$group)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

	input.mean <- output.mean <- matrix(0,ntags,nlibs)
	for(i in 1:length(levs.group)) {
		j <- object$samples$group==levs.group[i]
		beta <- mglmOneGroup(object$counts[,j,drop=FALSE],dispersion=dispersion,offset=log(lib.size[j]))
		lambda <- exp(beta)
		input.mean[,j] <- matrix(lambda,ncol=1) %*% matrix(lib.size[j],nrow=1)
		output.mean[,j] <- matrix(lambda, ncol=1) %*% matrix(common.lib.size, nrow=1, ncol=sum(j));
	}
	pseudo <- q2qnbinom(object$counts, input.mean=input.mean, output.mean=output.mean, dispersion=dispersion)
	pseudo[pseudo<0] <- 0
	list(pseudo.counts=pseudo, common.lib.size=common.lib.size)
}
