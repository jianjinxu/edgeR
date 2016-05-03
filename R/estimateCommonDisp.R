# Estimate common dispersion using exact conditional likelihood

estimateCommonDisp <- function(y, ...)
UseMethod("estimateCommonDisp")

estimateCommonDisp.DGEList <- function(y, tol=1e-06, rowsum.filter=5, verbose=FALSE, ...)
# Yunshun Chen. Created 18 March 2016.
{
	y <- validDGEList(y)
	group <- y$samples$group
	lib.size <- y$samples$lib.size * y$samples$norm.factors
	
	if( all(tabulate(group)<=1) ) {
		warning("There is no replication, setting dispersion to NA.")
		y$common.dispersion <- NA
		return(y)
	}

	out <- estimateCommonDisp(y$counts, group=group, lib.size=lib.size, tol=tol, rowsum.filter=rowsum.filter, verbose=verbose, ...)	
	y$common.dispersion <- out
	y <- equalizeLibSizes(y, dispersion=out)
	y$AveLogCPM <- aveLogCPM(y, dispersion=out)
	y
}


estimateCommonDisp.default <- function(y, group=NULL, lib.size=NULL, tol=1e-06, rowsum.filter=5, verbose=FALSE, ...)
#	Davis McCarthy, Mark Robinson, Gordon Smyth.
#	Created 2009. Last modified 18 March 2016.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check group
	if(is.null(group)) group <- rep(1, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")

	sel <- rowSums(y) > rowsum.filter
	
#	Start from small dispersion
	disp <- 0.01
	for(i in 1:2) {
		out <- equalizeLibSizes(y, group=group, dispersion=disp, lib.size=lib.size)
		y.pseudo <- out$pseudo.counts[sel, , drop=FALSE]
		y.split <- splitIntoGroups(y.pseudo, group=group)
		delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y.split, der=0)
		delta <- delta$maximum
		disp <- delta/(1-delta)
	}
	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	
	common.dispersion <- disp
	common.dispersion
}

