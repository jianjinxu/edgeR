dispCoxReid <- function(y, design=NULL, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000, AveLogCPM=NULL)
#	Cox-Reid APL estimator of common dispersion
#	Gordon Smyth, Davis McCarthy
#	26 Jan 2011.  Last modified 4 Feb 2013.
{
#	Check y
	y <- as.matrix(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}

#	Check offseet
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	if(min(interval)<0) stop("please give a non-negative interval for the dispersion")

#	Apply min row count
	small.row.sum <- rowSums(y)<min.row.sum
	if(any(small.row.sum)) {
		y <- y[!small.row.sum,,drop=FALSE]
		offset <- offset[!small.row.sum,,drop=FALSE]
		if(!is.null(AveLogCPM)) AveLogCPM <- AveLogCPM[!small.row.sum]
	}
	if(nrow(y)<1) stop("no data rows with required number of counts")

#	Subsetting
	if(!is.null(subset) && subset<=nrow(y)/2) {
		if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,offset=offset)
		i <- systematicSubset(subset,AveLogCPM)
		y <- y[i,,drop=FALSE]
		offset <- offset[i,,drop=FALSE]
	}

#	Function for optimizing
	fun <- function(par,y,design,offset) {
		sum(adjustedProfileLik(par^4,y,design,offset))
	}

	out <- optimize(f=fun,interval=interval^0.25,y=y,design=design,offset=offset,maximum=TRUE,tol=tol)
	out$maximum^4
}
