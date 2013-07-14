dispBinTrend <- function(y, design=NULL, offset=NULL, df=5, span=0.3, min.n=400, method.bin="CoxReid", method.trend="spline", AveLogCPM=NULL, ...)
#	Estimate common dispersion in bins based on AveLogCPM,
#	then fit a curve through the dispersions
#	Davis McCarthy, Gordon Smyth
#	Created 10 Feb 2011.  Last modified 17 April 2013.
{
#	Check y
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	pos <- rowSums(y)>0
	if(!any(pos)) return(AveLogCPM=AveLogCPM, dispersion=rep(0,ntags))
	npostags <- sum(pos)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,nlibs,1)
	} else {
		design <- as.matrix(design)
	}

#	Check offset
	if(is.null(offset)) offset <- expandAsMatrix(log(colSums(y)),dim(y))

#	Check methods
	method.bin <- match.arg(method.bin, c("CoxReid", "Pearson", "deviance"))
	method.trend <- match.arg(method.trend, c("spline", "loess"))

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y)

#	Define bins of genes based on min.n in each bin
#	All zero rows are marked as group==0
	group <- as.numeric(pos)

	if(npostags < 100)
		nbins <- 1
	else {
		nbins <- floor(npostags^0.4)
		nbins <- min(nbins,1000)
		min.n <- min(min.n,floor(npostags/nbins))
	}
	if(min.n < 50) {
		nbins <- floor(npostags/50)
		min.n <- 50
	}

#	nbins <- floor(npostags/min.n)
#	nbins <- min(max(nbins,1),1000)

	if(nbins>1) group[pos] <- cutWithMinN(AveLogCPM[pos],intervals=nbins,min.n=min.n)$group

#	Estimate dispersion in each bin
	bin.d <- bin.A <- rep(0,nbins)
	for(i in 1:nbins) {
		bin <- group==i
		bin.d[i] <- estimateGLMCommonDisp(y[bin,], design, method=method.bin, offset[bin,], min.row.sum=0, ...)
		bin.A[i] <- mean(AveLogCPM[bin])
	}

#	If just one bin, trended dispersion is constant
	if(nbins==1) {
		dispersion <- rep.int(bin.d,ntags)
		return(list(AveLogCPM=AveLogCPM, dispersion=dispersion, bin.AveLogCPM=bin.A, bin.dispersion=bin.d))
	}

#	If few bins, use linear interpolation
	if(nbins<7) {
		f <- approxfun(bin.A,sqrt(bin.d),rule=2)
		dispersion <- f(AveLogCPM)^2
		return(list(AveLogCPM=AveLogCPM, dispersion=dispersion, bin.AveLogCPM=bin.A, bin.dispersion=bin.d))
	}

#	Spline smoother through binned dispersions
	if( method.trend=="spline" ) {
		require("splines")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(bin.A,p=p1)
		r <- range(bin.A)
		knots2 <- r[1]+p1*(r[2]-r[1])
		knots <- 0.3*knots1+0.7*knots2
		basisbins <- ns(bin.A,df=df,knots=knots,intercept=TRUE)
		beta <- coefficients(lm.fit(basisbins, sqrt(bin.d)))
		basisall <- predict(basisbins,newx=AveLogCPM)
		dispersion <- drop(basisall %*% beta)^2
	}

#	Loess smoother though binned dispersions
	if( method.trend=="loess" ) {
		fit <- loessFit(sqrt(bin.d), bin.A, span=span, iterations=1)
		f <- approxfun(bin.A, fit$fitted, rule=2)
		dispersion <- f(AveLogCPM)^2
	}

	list(AveLogCPM=AveLogCPM, dispersion=dispersion, bin.AveLogCPM=bin.A, bin.dispersion=bin.d)
}

