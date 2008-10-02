


setMethod("plotMA","deDGEList",
function(object,xlab="A",ylab="M",pch=19,...) {
  plot( (log2(object$ps$p1)+log2(object$ps$p2))/2, log2(object$ps$p2/object$ps$p1), pch=pch, xlab=xlab,ylab=ylab,...)
})

