dispCoxReidInterpolateTagwise <- function(y, design, offset=NULL, dispersion, trend=TRUE, abundance=NULL, min.row.sum=5, prior.n=getPriorN(y, design), span=2/3, grid.npts=11, grid.range=c(-6,6))
#	Estimate tagwise NB dispersions
#	using weighted Cox-Reid Adjusted Profile-likelihood
#	and cubic spline interpolation over a tagwise grid.
#	Yunshun Chen and Gordon Smyth
#	Created August 2010. Last modified 3 Oct 2011.
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

#	Apply min.row.sum and use input dispersion for small count tags
	i <- (rowSums(y) >= min.row.sum)
	if(any(!i)) {
		if(any(i)) dispersion[i] <- Recall(y=y[i,],design=design,offset=offset[i,],dispersion=dispersion[i],abundance=abundance[i],grid.npts=grid.npts,min.row.sum=0,prior.n=prior.n,span=span,trend=trend)
		return(dispersion)
	}

#	Posterior profile likelihood
	spline.pts <- seq(from=grid.range[1],to=grid.range[2],length=grid.npts)
	apl <- matrix(0, nrow=ntags, ncol=grid.npts)
	for(i in 1:grid.npts){
		spline.disp <- dispersion * 2^spline.pts[i]
		apl[,i] <- adjustedProfileLik(spline.disp, y=y, design=design, offset=offset)
	}
	if(trend==FALSE) {
		apl.smooth <- matrix(colMeans(apl),ntags,grid.npts,byrow=TRUE)
	} else {
		o <- order(abundance)
		oo <- order(o)
		width <- min(1000,ntags)
		apl.smooth <- movingAverageByCol(apl[o,],width=1000)[oo,]
	}
	apl.smooth <- (apl+prior.n*apl.smooth)/(1+prior.n)

#	Tagwise maximization
	d <- dispersion
	for(j in 1:ntags) d[j] <- maximizeInterpolant(spline.pts, apl.smooth[j,])

	dispersion * 2^d
}

#dispCoxReidInterpolateTagwise0 <- function(y, design, offset=NULL, dispersion, method="trend", abundance=NULL, min.row.sum=5, prior.n=getPriorN(y, design), span=2/3, grid.npts=11, grid.range=c(-6,6))
##	Estimate tagwise NB dispersions
##	using weighted Cox-Reid Adjusted Profile-likelihood
##	and cubic spline interpolation over a tagwise grid.
##	Yunshun Chen and Gordon Smyth
##	Created August 2010. Last modified 09 Aug 2011.
#{
##	Check input arguments
#	y <- as.matrix(y)
#	ntags <- nrow(y)
#	nlibs <- ncol(y)
#	design <- as.matrix(design)
#	if(!is.fullrank(design)) stop("design matrix must be full column rank")
#	if(is.null(offset)) offset <- 0
#	offset <- expandAsMatrix(offset,dim(y))
#	ldisp <- length(dispersion)
#	if(ldisp==1) {
#		dispersion <- rep(dispersion,ntags)
#	} else {
#		if(ldisp != ntags) stop("length of dispersion doesn't match nrow(y)")
#	}
#	if(is.null(abundance)) abundance <- mglmOneGroup(y,offset=offset)
#
##	Apply min.row.sum and use input dispersion for small count tags
#	i <- (rowSums(y) >= min.row.sum)
#	if(any(!i)) {
#		if(any(i)) dispersion[i] <- Recall(y=y[i,],design=design,offset=offset[i,],dispersion=dispersion[i],abundance=abundance[i],grid.npts=grid.npts,min.row.sum=0,prior.n=prior.n,span=span,method=method)
#		return(dispersion)
#	}
#
##	Posterior profile likelihood
#	spline.pts <- seq(from=grid.range[1],to=grid.range[2],length=grid.npts)
#	abundance.rank <- rank(abundance)
#	apl <- apl.smooth <- matrix(0, nrow=ntags, ncol=grid.npts)
#
#	for(i in 1:grid.npts){
#		spline.disp <- dispersion * 2^spline.pts[i]
#		apl[,i] <- adjustedProfileLik(spline.disp, y=y, design=design, offset=offset)
#		if(method=="common"){
#			apl.smooth[,i] <- mean(apl[,i])
#		} else if(method=="trend"){
#			apl.smooth[,i] <- loessFit(apl[,i],abundance.rank,span=span,iterations=1)$fitted
#		} else stop("Wrong method")
#	}
#	apl.smooth <- (apl+prior.n*apl.smooth)/(1+prior.n)
#
##	Tagwise maximization
#	d <- dispersion
#	for(j in 1:ntags) d[j] <- maximizeInterpolant(spline.pts, apl.smooth[j,])
#
#	dispersion * 2^d
#}
