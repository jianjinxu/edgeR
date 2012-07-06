loessByCol <- function(y, x=NULL, span=0.5)
# Calls a C++ function to do the dirty work of fitting a degree-0,
# non-robustified LOWESS curve through each column of a matrix
# of data points.
# C++ version by Aaron Lun, 26 June 2012.
# Replaces:
# Rcode version by Davis McCarthy, May 2010.
# loess version by Yunshun Chen, 08 May 2012.
{
    y <- as.matrix(y)
    ntags <- nrow(y)
    if(is.null(x)) x <- 1:ntags

    # Sorting by the x-value.
    x.order<-order(x)
    y<-y[x.order,,drop=FALSE]
    x<-x[x.order]

    nspan <- floor(span*ntags)
    nspan <- min(ntags,nspan)
    if(nspan<=1) {
    	return(list(fitted.values=y,leverages=rep(1,ntags)))
    }

    # Passing to the compiled code.
    fitted<-.Call("lowess_by_col", x, y, ncol(y), nspan, PACKAGE="edgeR")
   
    # Unsorting them to recover the original order.    
    fitted[[1]][x.order,]<-fitted[[1]]
    fitted[[2]][x.order]<-fitted[[2]]

    # Beautifying.
    names(fitted)<-c("fitted.values", "leverages")

    return(fitted)
}
