dispCoxReid <- function(y, design, offset=NULL, interval=c(0,4), tol=0.001, ...) 
UseMethod("dispCoxReid")

dispCoxReid.DGEList <- function(y, design, offset=NULL, interval=c(0,4), tol=0.001)
{
	if(is.null(offset)) offset <- getOffsets(y)
	dispCoxReid(y=y$counts, design=design, offset=offset, interval=interval)
}

dispCoxReid.default <- function(y, design, offset=0, interval=c(0,4), tol=0.001)
#	Estimate common dispersion by optimize
#	Gordon Smyth
#	26 Jan 2011.  Last modified 19 Jan 2011.
{
	y <- as.matrix(y)
	offset <- expandAsMatrix(offset,dim(y))
	
	fun <- function(par,y,design,offset) {
		tryCatch(sum(adjustedProfileLik(par,y,design,offset),na.rm=TRUE),error=function(e) -1e10)
	}

	out <- optimize(f=fun,interval=interval,y=y,design=design,offset=offset,maximum=TRUE,tol=tol)
	out$maximum
}
