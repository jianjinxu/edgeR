require(methods)

setClass("deDGEList",
#  Linear model fit
representation("list")
)

setClass("DGEList",
#  Linear model fit
representation("list")
)

setClass("EBList",
#  Linear model fit
representation("list")
)

setMethod("show","DGEList",
  function(object) {
    cat(class(object),": ",nrow(object$data)," rows, ",ncol(object$data)," libraries\n",sep="")
  })

setMethod("show","EBList",
  function(object) {
    cat(class(object),": alpha=",object$alpha,"\n",sep="")
  })

setMethod("show","deDGEList",
  function(object) {
    cat(class(object),": ",ncol(ms$pseudo)," samples, adjusted to library size of ",object$M,"\n",sep="")
  })

