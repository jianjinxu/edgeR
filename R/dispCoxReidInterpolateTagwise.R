dispCoxReidInterpolateTagwise <- function(y, design, offset=NULL, dispersion, abundance=NULL, npts=11, min.row.sum=5, prior.n=10, span=0.3)
#	Estimate tagwise NB dispersions
#	using weighted Cox-Reid Adjusted Profile-likelihood
#	and cubic spline interpolation over a tagwise grid.
#	Yunshun Chen and Gordon Smyth
#	Created August 2010. Last modified 12 Feb 2011.
{
#	Check input arguments
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)
	design <- as.matrix(design)
	if(!is.fullrank(design)) stop("design matrix must be full column rank")
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	ldisp <- length(dispersion)
	if(ldisp==1) {
		dispersion <- rep(dispersion,ntags)
	} else {
		if(ldisp != ntags) stop("length of dispersion doesn't match nrow(y)")
	}
	if(is.null(abundance)) abundance <- mglmOneGroup(y,offset=offset)

#	Apply rowsum.filter and use input dispersion for small count tags
	i <- (rowSums(y) >= min.row.sum)
	if(any(!i)) {
		if(any(i)) dispersion[i] <- Recall(y=y[i,],design=design,offset=offset[i,],dispersion=dispersion[i],abundance=abundance[i],npts=npts,rowsum.filter=0,prior.n=prior.n,span=span)
		return(dispersion)
	}

#	Posterior profile likelihood
	spline.pts <- seq(from=-4,to=4,length=npts)
	abundance.rank <- rank(abundance)
	apl <- apl.smooth <- matrix(0, nrow=ntags, ncol=npts)
	for(i in 1:npts){
		spline.disp <- dispersion * 2^spline.pts[i]
		apl[,i] <- adjustedProfileLik(spline.disp, y=y, design=design, offset=offset)
		apl.smooth[,i] <- loessFit(apl[,i],abundance.rank,span=span,iterations=1)$fitted
	}
	apl.smooth <- (apl+prior.n*apl.smooth)/(1+prior.n)

#	Tagwise maximization
	d <- dispersion
	for(j in 1:ntags) d[j] <- maximizeInterpolant(spline.pts, apl.smooth[j,])

	dispersion * 2^d
}

