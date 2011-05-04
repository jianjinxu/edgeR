dispDeviance <- function(y, design, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=1000, robust=FALSE, trace=FALSE)
#	Deviance estimator of common dispersion
#	Gordon Smyth, Davis McCarthy
#	26 Jan 2011.  Last modified 4 May 2011.
{
	y <- as.matrix(y)
	design <- as.matrix(design)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	if(min(interval)<0) stop("please give a non-negative interval for the dispersion")
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

	Df <- ncol(y)-ncol(design)
	if(robust) {
		p <- pchisq(Df, df=Df)
		bias <- function(x) quantile(x,p)-Df
	} else {
		bias <- function(x) mean(x)-Df
	}

	fun <- function(par,y,design,offset) {
		fit <- glmFit(y,design,dispersion=par^4,offset=offset)
		if(trace) cat(par^4,bias(fit$deviance),"\n")
		bias(fit$deviance)
	}

	if(fun(interval[1],y,design,offset)<=0) {
		return(interval[1])
	}

	if(fun(interval[2],y,design,offset)>0) {
		warning("dispersion estimate above interval upper limit")
		return(interval[2])
	}

	if(trace) cat("Dispersion, mean(deviance)-df\n")
	out <- uniroot(f=fun,interval=interval^0.25,y=y,design=design,offset=offset,tol=tol)
	out$root^4
}
