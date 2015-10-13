plotMD.DGEList <- function(object, column=1, xlab="Average log CPM (this sample and others)", ylab="log-ratio (this sample vs others)", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, prior.count=3, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015. Last modified 24 June 2015.
{
	nlib <- ncol(object)
	if(nlib < 2L) stop("Need at least two columns")

#	Convert column to integer if not already
	j <- 1:nlib
	names(j) <- colnames(object)
	column <- j[column[1]]

	logCPM <- cpm(object, log=TRUE, prior.count=prior.count)
	AveOfOthers <- rowMeans(logCPM[,-column,drop=FALSE],na.rm=TRUE)
	Diff <- logCPM[,column]-AveOfOthers
	Mean <- (logCPM[,column]+AveOfOthers)/2

	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		Diff[ is.na(w) | (w <= 0) ] <- NA
	}

	plotWithHighlights(x=Mean,y=Diff,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.DGEGLM <- function(object, column=ncol(object), coef=NULL, xlab="Average log CPM", ylab="log-fold-change", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015. Last modified 24 June 2015.
{
	if(!is.null(coef)) column <- coef
	if(is.null(object$AveLogCPM)) stop("AveLogCPM component is absent.")
	logFC <- as.matrix(object$coefficients)[,column]
	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		logFC[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=object$AveLogCPM,y=logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.DGEExact <- plotMD.DGELRT <- function(object, xlab="Average log CPM", ylab="log-fold-change", main=object$comparison, status=object$genes$Status, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Last modified 24 June 2015.
{
	plotWithHighlights(x=object$table$logCPM,y=object$table$logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}
