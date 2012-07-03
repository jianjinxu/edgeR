# Created March 2011. Last modified 3 July 2012.

estimateGLMTagwiseDisp <- function(y, design=NULL, offset=NULL, dispersion=NULL, trend=TRUE, prior.df=20, span=NULL, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design=NULL, offset=NULL, dispersion=NULL, trend=!is.null(y$trended.dispersion), prior.df=20, span=NULL, ...)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(is.null(dispersion)) if(trend) {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) stop("No trended.dispersion found in data object. Run estimateGLMTrendedDisp first.")
		if(is.null(y$abundance)) y$abundance <- mglmOneGroup(y,offset=offset)
	} else {
		dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No common.dispersion found in data object. Run estimateGLMCommonDisp first.")
	}
	out <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=offset, dispersion=dispersion, trend=trend, abundance=y$abundance, prior.df=prior.df, ...)
	y$prior.df <- prior.df
	y$span <- out$span
	y$tagwise.dispersion <- out$tagwise.dispersion
	y
}

estimateGLMTagwiseDisp.default <- function(y, design=NULL, offset=NULL, dispersion, trend=TRUE, prior.df=20, span=NULL, ...)
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(ntags==0) return(numeric(0))
	nlibs <- ncol(y)
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: setting dispersion to NA")
		return(rep(NA,ntags))
	}
#	prior.n <- getPriorN(y=y,design=design,prior.df=prior.df)
	if(is.null(span)) if(ntags>10) span <- (10/ntags)^0.23 else span <- 1
	tagwise.dispersion <- dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, prior.df=prior.df, span=span, ...)
	list(tagwise.dispersion=tagwise.dispersion,span=span)
}
