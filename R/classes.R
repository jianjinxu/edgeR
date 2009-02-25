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
    cat(class(object),": ",ncol(object$pseudo)," samples, ",length(levels(object$group))," groups, adjusted to library size of ",object$M,"\n",sep="")
  })

#setGeneric("plotMA", function(object, pair=c(1,2), xlab = "A", ylab = "M", ylim=NULL, pch = 19, ...) standardGeneric("plotMA"))


DGEList <- function(data=matrix(0), lib.size=integer(0), group=factor(), ...)
{
  if (ncol(data) != length(lib.size))
    stop("Length of 'lib.size' must equal number of columns in 'data'")
  if (ncol(data) != length(group))
    stop("Length of 'group' must equal number of columns in 'data'")
  if (!is.factor(group))
    group<-as.factor(group)
  new("DGEList",list(data=as.matrix(data), lib.size=lib.size, group=group, ...))
}


getData <- function(object)
{
  object$data
}

