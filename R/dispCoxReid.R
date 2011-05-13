dispCoxReid <- function(y, design, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000)
#	Cox-Reid APL estimator of common dispersion
#	Gordon Smyth, Davis McCarthy
#	26 Jan 2011.  Last modified 24 Mar 2011.
{
	y <- as.matrix(y)
	design <- as.matrix(design)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	small.row.sum <- rowSums(y)<min.row.sum
	if(any(small.row.sum)) {
		y <- y[!small.row.sum,,drop=FALSE]
		offset <- offset[!small.row.sum,,drop=FALSE]
	}
	if(nrow(y)<1) stop("no data rows with required number of counts")
	if(!is.null(subset) && subset<=nrow(y)/2) {
		A <- mglmOneGroup(y,offset=offset)
		i <- systematicSubset(subset,A)
		y <- y[i,,drop=FALSE]
		offset <- offset[i,,drop=FALSE]
	}

	fun <- function(par,y,design,offset) {
		sum(adjustedProfileLik(par^4,y,design,offset))
	}

	out <- optimize(f=fun,interval=interval^0.25,y=y,design=design,offset=offset,maximum=TRUE,tol=tol)
	out$maximum^4
}
