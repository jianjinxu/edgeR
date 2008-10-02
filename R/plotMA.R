


setMethod("plotMA","deDGEList",
function(object,xlab="A",ylab="M",ylim=NULL,pch=19,...) {
  M<-log2(object$ps$p2/object$ps$p1)
  A<-(log2(object$ps$p1)+log2(object$ps$p2))/2
  if(is.null(ylim))
    ylim<-c(-1,1)*max(abs(M))
  plot( A, M, pch=pch, xlab=xlab,ylab=ylab,ylim=ylim,...)
})

