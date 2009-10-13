deDGE <- function(object, alpha = 500, doPoisson = FALSE, verbose = TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to carry out the DE analysis of DGE data
{
    if (!is(object, "DGEList")) 
        stop("Currently supports DGEList objects")
    object$counts <- as.matrix(object$counts)
    if (doPoisson) {
        if (verbose) 
            cat("Quantile adjusting as Poisson.\n")
        qA <- quantileAdjust(object, r.init = 1000, n.iter = 1)
        qA$r <- rep(1000, length(qA$r))
    }
    else {
        if (verbose) 
            cat("Calculating shrinkage overdispersion parameters.\n")
        qA <- quantileAdjust(object, alpha = alpha, verbose = verbose)
    }
    rownames(qA$pseudo) <- rownames(object$counts)
    colnames(qA$pseudo) <- paste("pseudo", colnames(object$counts), sep = ".")
    new("deDGEList", (list(ps = qA$ps, r = qA$r, pseudo = qA$pseudo, group = object$samples$group, M = qA$N)))
}