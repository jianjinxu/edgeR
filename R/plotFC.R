plotFC<-function(object,xlab="log2Conc",ylab="log2FC",ylim=NULL,pch=19,cex=0.2,...) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009, revised August 2009
# A function/method to plot differential expression against overall expression (equivalent to plots used in microarray analyis)
{
	if(is.null(ylim))
		ylim<-c(-1,1)*max(abs(object$table$logFC))
	plot(object$table$logConc, object$table$logFC,pch=pch, cex=cex, xlab=xlab,ylab=ylab,ylim=ylim,panel.first=grid(),...)
}

