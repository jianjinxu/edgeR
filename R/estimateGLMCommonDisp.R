#  Last modified 13 March 2013

estimateGLMCommonDisp <- function(y, design=NULL, offset=NULL, method="CoxReid", subset=10000, AveLogCPM=NULL, verbose=FALSE, ...) 
UseMethod("estimateGLMCommonDisp")

estimateGLMCommonDisp.DGEList <- function(y, design=NULL, offset=NULL, method="CoxReid", subset=10000, AveLogCPM=NULL, verbose=FALSE, ...)
{
#	If provided as arguments, offset and AveLogCPM over-rule the values stored in y
	if(!is.null(AveLogCPM)) y$AveLogCPM <- AveLogCPM
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	if(!is.null(offset)) y$offset <- expandAsMatrix(offset,dim(y))

	disp <- estimateGLMCommonDisp(y=y$counts, design=design, offset=getOffset(y), method=method, subset=subset, AveLogCPM=y$AveLogCPM, verbose=verbose, ...)
	y$common.dispersion <- disp
	y$design <- design
	y
}

estimateGLMCommonDisp.default <- function(y, design=NULL, offset=NULL, method="CoxReid", subset=10000, AveLogCPM=NULL, verbose=FALSE, ...)
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
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: setting dispersion to NA")
		return(NA)
	}

#	Check method
	method <- match.arg(method, c("CoxReid","Pearson","Pearson2","deviance"))

#	Check offset
	if(is.null(offset)) offset <- log(colSums(y))

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y)

#	Call lower-level function
	disp <- switch(method,
		CoxReid=dispCoxReid(y, design=design, offset=offset, subset=subset, AveLogCPM=AveLogCPM, ...),
		Pearson=dispPearson(y, design=design, offset=offset, subset=subset, AveLogCPM=AveLogCPM, ...),
		deviance=dispDeviance(y, design=design, offset=offset, subset=subset, AveLogCPM=AveLogCPM, ...)
	)
	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	disp
}

