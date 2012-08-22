### BIN.DISPERSION.R

binCMLDispersion <- function(y, nbins=50)
#	Estimate common dispersion in bins based on abundance
#	Created Oct 2010. Last modified 2 Aug 2012.
{
	if(!is(y,"DGEList")) stop("y must be DGEList object.")
	ntags <- nrow(y)
	if(nbins>ntags) stop("nbins greater than number of rows of data")

	logCPM <- y$logCPM
	if(is.null(logCPM)) {
		logCPM <- mglmOneGroup(y$counts,offset=getOffset(y),dispersion=0.1)
		logCPM <- log2(exp(logCPM+log(1e6))+0.5)
	}
	bins <- cutWithMinN(logCPM,intervals=nbins,min.n=floor(ntags/nbins))

	disp.bins <- logCPM.bins <- rep(NA,nbins)
	for(i in 1:nbins) {
		tagsinbin <- bins$group==i
		disp.bins[i] <- estimateCommonDisp(y[tagsinbin,])$common.dispersion
		logCPM.bins[i] <- mean(logCPM[tagsinbin])
	}
	list(dispersion.bins=disp.bins, logCPM.bins=logCPM.bins)
}


binGLMDispersion <- function(y, design, offset=NULL, min.n=100, method="CoxReid", abundance=NULL, ... )
#	Estimate common dispersion in bins based on abundance.
#	Davis McCarthy, Gordon Smyth.
#	Created 7 Feb 2011. Last modified 3 Aug 2012.
{
	if(is(y,"DGEList")) {
		if(is.null(offset)) offset <- getOffset(y)
		y <- y$counts
	} else {
		y <- as.matrix(y)
		lib.size <- colSums(y)
		if(is.null(offset)) offset <- log(lib.size)
	}
	if(is.null(abundance)) {
		abundance <- mglmOneGroup(y,offset=offset,dispersion=0.01)
		abundance <- log2(exp(abundance+log(1e6))+0.5)
	}
	offset <- expandAsMatrix(offset,dim(y))
	method <- match.arg(method, c("CoxReid", "Pearson", "deviance"))

#	Remove all zero rows
	all.zero <- rowSums(y)==0
	if(any(all.zero)) {
		y <- y[!all.zero,,drop=FALSE]
		offset <- offset[!all.zero,,drop=FALSE]
		abundance <- abundance[!all.zero]
	}

	ntags <- nrow(y)
	if(ntags==0) return(list(dispersion=numeric(0),abundance=numeric(0)))
 
#	Define bins of genes based on min.n in each bin
	nbins <- floor(ntags/min.n)
	nbins <- min(max(nbins,1),1000)
	if(nbins==1)
		group <- rep(1,ntags)
	else
		group <- cutWithMinN(abundance, intervals=nbins, min.n=min.n)$group

#	Estimate dispersion in each bin
	dispersion <- ave.abundance <- rep(0,nbins)
	for(i in 1:nbins) {
		bin <- group==i
		dispersion[i] <- estimateGLMCommonDisp(y[bin,], design, method=method, offset[bin,], min.row.sum=0, ...)
		ave.abundance[i] <- mean(abundance[bin])
	}
	
	list(dispersion=dispersion, abundance=ave.abundance)
}
