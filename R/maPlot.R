maPlot <- function(x,y, normalize=FALSE, smearWidth=1, col=NULL, allCol="black", lowCol="orange", deCol="red",de.tags=NULL, ...) {
    if(normalize) {
      x <- x/sum(x)
      y <- y/sum(y)
    }
    A <- (log2(x)+log2(y))/2
    M <- log2(y) - log2(x)
    if( range(x)[2] > 1 ) {
      # assume x and y are counts
      w <- x==min(x) | y==min(y)
      A[w] <- min(A[!w]) - runif(sum(w),min=0,max=smearWidth)
      M[w] <- log2(y[w]+min(y[!w])) - log2(x[w]+min(x[!w])) 
    } else {
      # assume x and y are proportions
      w <- x==min(x) | y==min(y)
      A[w] <- min(A[!w]) - runif(sum(w),min=0,max=smearWidth)
      M[w] <- log2(y[w]+min(y[!w])) - log2(x[w]+min(x[!w])) 
    }
    if( is.null(col) ) {
      col <- rep(allCol, length(A))
      col[w] <- lowCol
    }
    plot(A,M,col=col,...)
    points(A[de.tags],M[de.tags],col=deCol,...)
    invisible(list(A=A,M=M,w=w))
}