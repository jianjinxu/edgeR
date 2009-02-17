# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function/method to plot differential expression against overall expression (equivalent to plots used in microarray analyis)
plotMA<-function(object,pair=c(1,2),xlab="A",ylab="M",ylim=NULL,pch=19,...) {
	g1<-pair[1]
	g2<-pair[2]
	M<-log2(object$ps$p.group[,g2]/object$ps$p.group[,g1])
	A<-(log2(object$ps$p.group[,g1])+log2(object$ps$p.group[,g2]))/2
	if(is.null(ylim))
		ylim<-c(-1,1)*max(abs(M))
	plot( A, M, pch=pch, xlab=xlab,ylab=ylab,ylim=ylim,...)
}



# A function/method to plot differential expression against overall expression (equivalent to plots used in microarray analyis)
#setMethod("plotMA","deDGEList",
#function(object,xlab="A",ylab="M",ylim=NULL,pch=19,...) {
#	M<-log2(object$ps$p2/object$ps$p1)
#	A<-(log2(object$ps$p1)+log2(object$ps$p2))/2
#	if(is.null(ylim))
#		ylim<-c(-1,1)*max(abs(M))
#	plot( A, M, pch=pch, xlab=xlab,ylab=ylab,ylim=ylim,...)
#})

