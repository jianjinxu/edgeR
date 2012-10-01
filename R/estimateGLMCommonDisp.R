#  Last modified 1 Oct 2012

estimateGLMCommonDisp <- function(y, design=NULL, offset=NULL, method="CoxReid", verbose=FALSE, ...) 
UseMethod("estimateGLMCommonDisp")

estimateGLMCommonDisp.DGEList <- function(y, design=NULL, offset=NULL, method="CoxReid", verbose=FALSE, ...)
{
    if( is.null(offset) )
        offset <- getOffset(y)
	y$abundance <- mglmOneGroup(y$counts,offset=getOffset(y),dispersion=0.05)
	y$logCPM <- log1p(exp(y$abundance+log(1e6)))/log(2)
	d <- estimateGLMCommonDisp(y=y$counts, design=design, offset=offset, method=method, verbose=verbose, ...)
	y$common.dispersion <- d
	y
}

estimateGLMCommonDisp.default <- function(y, design=NULL, offset=NULL, method="CoxReid", verbose=FALSE, ...)
{
	y <- as.matrix(y)
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: setting dispersion to NA")
		return(NA)
	}
	method <- match.arg(method, c("CoxReid","Pearson","deviance"))
	disp <- switch(method,
		CoxReid=dispCoxReid(y, design=design, offset=offset, ...),
		Pearson=dispPearson(y, design=design, offset=offset, ...),
		deviance=dispDeviance(y, design=design, offset=offset, ...)
	)
	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	disp
}

