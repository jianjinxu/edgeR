maPlot <- function(x,y, logAbundance=NULL, logFC=NULL, normalize=FALSE, smearWidth=1, col=NULL, allCol="black", lowCol="orange", deCol="red",de.tags=NULL, smooth.scatter=FALSE, lowess=FALSE, ...) {
    ## Low-level function for creating an MA-plot for DGE data.
    ## Created by Mark Robinson. Last modified by Davis McCarthy, 19 November 2010.
    if( !is.null(logAbundance) & !is.null(logFC) ) {
        A <- logAbundance
        M <- logFC
        w <- rep(FALSE, length(A))
        w <- A < -25
        if( any(w) ) {
            shift <- max(abs(M[w])) - max(abs(M[!w]))
            A[w] <- min(A[!w]) - runif(sum(w),min=0,max=smearWidth)
            M[w] <- sign(M[w]) * (abs(M[w]) - shift)
        }
    } else {
        if(normalize) {
            x <- x/sum(x)
            y <- y/sum(y)
        }
        A <- (log2(x)+log2(y))/2
        M <- log2(y) - log2(x)
        w <- x==min(x) | y==min(y)
        if( any(w) ) {
            A[w] <- min(A[!w]) - runif(sum(w),min=0,max=smearWidth)
            M[w] <- log2(y[w]+min(y[!w])) - log2(x[w]+min(x[!w]))
        }
    }
    if( is.null(col) ) {
      col <- rep(allCol, length(A))
      if( any(w) )
          col[w] <- lowCol
    }
    if(smooth.scatter) {
        smoothScatter(A, M, col=col, ...)
        grid()
        if( any(w) )
            points(A[w], M[w], col=lowCol, ...)
    }
    else
        plot(A,M,col=col,...)
    points(A[de.tags],M[de.tags],col=deCol,...)
    if(lowess) {
        keep <- A > min(A[!w]) + 1
        low <- lowess(A[keep],M[keep], f=1/4)
        lines(low,col="red",lwd=4)
    }
    invisible(list(A=A,M=M,w=w))
}

