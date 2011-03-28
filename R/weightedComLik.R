## Functions to calculate the weights for the common likelihood when calculating tagwise dispersions and allowing for dependence in the dispersion estimate on tag abundance

weightedComLik <- function(object,l0,prop.used=0.25) {
    ## Function for calculating weights to do LOESS-like weighted local calculations of the common dispersion
    ## Written by Davis McCarthy, May 2010. Last modified 3 June 2010.
    ## We order the tags based on their average abundance across all groups
    ## l0 is a matrix of ntags rows and number of columns given by grid.length
    if(is.null(object$conc))
        stop("estimateCommonDisp() must be run before using this function.\n")
    o <- order(object$conc$conc.common)
    x <- object$conc$conc.common
    xord <- x[o]
    ntags <- nrow(object$counts)
    ntags.either.side <- ceiling(prop.used*ntags/2)
    weighted.common.lik <- matrix(0,nrow=ntags,ncol=ncol(l0))
    for(i in 1:ntags) {
        wts.vec <- rep(0, length=ntags)
        if(i < ntags.either.side+2)
            tags.used <- 1:(i+ntags.either.side)
        else {
            if(i > ntags-ntags.either.side)
                tags.used <- (i-ntags.either.side):ntags
            else
                tags.used <- (i-ntags.either.side):(i+ntags.either.side)
        }
        d <- max(abs(xord[i]-xord[tags.used]))
        tag.weight <- .tricube( (xord[tags.used]-xord[i])/d )
        #weights[o[i],o[tags.used]] <- tag.weight/sum(tag.weight)
        wts.vec[o[tags.used]] <- tag.weight/sum(tag.weight)
        #weighted.common.lik[o[i],] <-  weights[o[i],] %*% l0 ## or colSums(weights[o[i],] * l0)
        weighted.common.lik[o[i],] <-  wts.vec %*% l0 # colSums(wts.vec * l0) #
    }
    weighted.common.lik
}

.tricube <- function(x) {
    (1 - abs(x)^3)^3
}

weightedComLikMA <- function(object, l0, prop.used=0.05) {
    ## Function for calculating weights to do LOESS-like weighted local calculations of the common dispersion
    ## Written by Davis McCarthy, May 2010. Last modified 3 June 2010.
    ## We order the tags based on their average abundance across all groups
    ## l0 is a matrix of ntags rows and number of columns given by grid.length
    if(is.null(object$conc))
        stop("estimateCommonDisp() must be run before using this function.\n")
    o <- order(object$conc$conc.common)
    ntags <- nrow(object$counts)
    weighted.common.lik <- matrix(0,nrow=ntags,ncol=ncol(l0))
    width <- floor(prop.used*ntags)
    weighted.common.lik <- movingAverageByCol(l0[o,], width=width)
    weighted.common.lik[o,] <- weighted.common.lik
    weighted.common.lik
}




.weightedCommonLik <- function(object,l0,prop.used=0.4) {
## Function for calculating weights to do LOESS-like weighted local calculations of the common dispersion
## Written by Davis McCarthy, May 2010. Last modified 3 June 2010.
    ## We order the tags based on their average abundance across all groups
    ## l0 is a matrix of ntags rows and number of columns given by grid.length
    o <- order(object$conc$conc.common)
    o2 <- order(o) ## Ordering to regain original order from ordered object
    a <- object$conc$conc.common[o] # Ordered abundance
    l0.ord <- l0[o,]
    ntags <- nrow(object$counts)
    weighted.common.lik <- apply(l0.ord,2,.lo.weights,a=a,prop.used=prop.used)
    weighted.common.lik[o2,]
}

.lo.weights <- function(y,a,prop.used) {
## Function to give lowess smooth for weights for the common likelihood
## Written by Davis McCarthy, May 2010. Last modified 3 June 2010.
    lowess(a, y, iter=0, f=prop.used)$y
}


