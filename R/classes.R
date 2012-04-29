require(methods)

setClass("DGEExact",
#  Linear model fit
representation("list")
)

setClass("DGEList",
#  Linear model fit
representation("list")
)

setClass("DGEGLM",
representation("list")
)

setClass("DGELRT",
representation("list")
)

setIs("DGEList","LargeDataObject")
setIs("DGEExact","LargeDataObject")
setIs("DGEGLM","LargeDataObject")
setIs("DGELRT","LargeDataObject")

dim.DGEList <- function(x) if (is.null(x$counts)) c(0, 0) else dim(as.matrix(x$counts))
dim.DGEExact <- dim.TopTags <- dim.DGEGLM <- dim.DGELRT <- function(x) if (is.null(x$table)) c(0, 0) else dim(as.matrix(x$table))

length.DGEList <- length.DGEExact <- length.TopTags <- length.DGEGLM <- length.DGELRT <- function(x) prod(dim(x))

dimnames.DGEList <- function(x) dimnames(x$counts)
assign("dimnames<-.DGEList",function(x,value)
{
	dimnames(x$counts) <- value
	if(!is.null(x$samples)) row.names(x$samples) <- value[[2]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

dimnames.DGEExact <- function(x) dimnames(x$coefficients)
assign("dimnames<-.DGEExact",function(x,value)
{
	dimnames(x$coefficients) <- value
	if(!is.null(x$samples)) row.names(x$samples) <- value[[2]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

dimnames.DGEGLM <- function(x) dimnames(x$coefficients)
assign("dimnames<-.DGEGLM",function(x,value)
{
	dimnames(x$coefficients) <- value
	if(!is.null(x$samples)) row.names(x$samples) <- value[[2]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

dimnames.TopTags <- function(x) dimnames(x$table)
assign("dimnames<-.TopTags",function(x,value)
{
	dimnames(x$table) <- value
	x
})

as.matrix.DGEList <- function(x,...) as.matrix(x$counts)

