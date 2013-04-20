# S3 as.matrix method

as.matrix.DGEList <- function(x,...) as.matrix(x$counts)

# S3 as.data.frame method

as.data.frame.TopTags <- function(x,row.names=NULL,optional=FALSE,...)
{
	if(!is.null(row.names)) row.names(x$table) <- row.names
	x$table
}

# S3 dim methods
# These enable nrow() and ncol() as well

dim.DGEList <- function(x) if(is.null(x$counts)) c(0,0) else dim(as.matrix(x$counts))
dim.DGEGLM <- function(x) if(is.null(x$coefficients)) c(0,0) else dim(as.matrix(x$coefficients))
dim.DGEExact <- dim.TopTags <- dim.DGELRT <- function(x) if(is.null(x$table)) c(0,0) else dim(as.matrix(x$table))

# S3 length methods

length.DGEList <- length.DGEExact <- length.TopTags <- length.DGEGLM <- length.DGELRT <- function(x) prod(dim(x))

# S3 dimnames methods
# These enable rownames() and colnames() as well

dimnames.DGEList <- function(x) dimnames(x$counts)
dimnames.DGEGLM <- function(x) dimnames(x$coefficients)
dimnames.DGEExact <- dimnames.DGELRT <- dimnames.TopTags <- function(x) dimnames(x$table)

# S3 dimnames<- methods
# These enable rownames()<- and colnames()<- as well

assign("dimnames<-.DGEList",function(x,value)
{
	dimnames(x$counts) <- value
	if(!is.null(x$samples)) row.names(x$samples) <- value[[2]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

assign("dimnames<-.DGEGLM",function(x,value)
{
	dimnames(x$coefficients) <- value
	if(!is.null(x$samples)) row.names(x$samples) <- value[[2]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

