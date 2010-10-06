binMeanVar <- function(x, conc=NULL, group, nbins=100) {
    ## Function to bin DGE data based on abundance and calculate the mean and pooled variance for each tag, as well as the average mean and variance for each bin. Allows us to investigate the mean-variance relationship in the data.
    ## Expect x to be a matrix of counts or pseudocounts---pseudocounts preferable as this adjusts for library size.
    means <- rowMeans(x)
    vars <- apply(x,1,pooledVar,group=group)
    means.quantiles <- quantile(means, probs=seq(0,1,length=nbins+1))
    if(!is.null(conc))
       conc.quantiles <- quantile(conc, probs=seq(0,1,length=nbins+1))
    if(any(duplicated(means.quantiles))) {
        cat("Duplicated quantiles for the means, so using quantiles of conc instead\n")
        if(is.null(conc))
            stop("conc is NULL, so cannot bin data.\n")
        else
            breaks <- conc.quantiles
    }
    else
        breaks <- means.quantiles
    f <- cut(means,breaks=breaks)
    bins <- split(1:nrow(x),f)
    var.bins <- means.bins <- list()
    for(i in 1:nbins){
        means.bins[[i]] <- means[bins[[i]]]
        var.bins[[i]] <- vars[bins[[i]]]
    }
    ave.means <- sapply(means.bins, mean)
    ave.vars <- sapply(var.bins, mean)
    list(avemeans=ave.means,avevars=ave.vars,bin.means=means.bins, bin.vars=var.bins, means=means, vars=vars)
}

pooledVar <- function(y,group) {
## Function to calculate the pooled variance for a vector of data with factor indicating group supplied
    numerator <- 0
    denominator <- 0
    for(i in 1:nlevels(group) ) {
        choosei <- as.logical(match(group,levels(group)[i]))
        choosei[is.na(choosei)] <- FALSE
        ni <- sum(choosei)
        numerator <- numerator + ni*var(y[choosei])
        denominator <- denominator + ni - 1
    }
    numerator/denominator
}

plotMeanVar <- function(object, meanvar=NULL, show.raw.vars=FALSE, show.tagwise.vars=FALSE, dispersion.method="coxreid", scalar=NULL, NBline=FALSE, nbins=100, ...) {
    ## Creates a mean-variance plot (with binned values) for a given DGEList object
    ## Uses the pooledVar and binMeanVar functions and operates on pseudo-counts to account for differences in library sizes
    if(!is(object,"DGEList"))
        stop("This function requires a DGEList object/n")
    if(!is.null(meanvar))
        if(is.null(meanvar$means) | is.null(meanvar$vars) | is.null(meanvar$avemeans) | is.null(meanvar$avevars)) {
            cat("Cannot extract all required elements of meanvar object, so will recompute it.\n")
            meanvar <- NULL
        }
    if(!is.null(scalar)) {
        if( length(scalar) != nrow(object$counts) )
            stop("The supplied argument scalar must have length equal to the number of columns of the count matrix in the DGEList object.")
        scalar <- scalar
    }
    else
        scalar <- object$samples$lib.size*object$samples$norm.factors/exp(mean(log(object$samples$lib.size*object$samples$norm.factors)))
    scalingmatrix <- outer(rep(1,nrow(object$counts)), scalar)
    x <- object$counts/scalingmatrix
    dispersion.method <- match.arg(dispersion.method,c("coxreid","qcml"))
    if(NBline | show.tagwise.vars) {
        if( dispersion.method=="coxreid" ) {
            common.dispersion <- object$CR.common.dispersion
            tagwise.dispersion <- object$CR.tagwise.dispersion
        }
        else {
            common.dispersion <- object$common.dispersion
            tagwise.dispersion <- object$tagwise.dispersion
        }
        if(is.null(common.dispersion)) {
            if(dispersion.method=="coxreid")
                stop("Could not extract Cox-Reid common dispersion. Try running CRDisp on the DGEList object before plotMeanVar.\n")
            else
                stop("Could not extract qCML common dispersion. Try running estimateCommonDisp on the DGEList object before plotMeanVar.\n")
        }
        lmu <- seq(1e-5,10,length.out=1000)
        nb.var <- 10^lmu + (10^lmu)^2*common.dispersion
    }
    if(is.null(meanvar) & dispersion.method=="qcml")
        meanvar <- binMeanVar(x, object$conc$conc.common, object$samples$group, nbins=nbins)
    if(is.null(meanvar) & dispersion.method=="coxreid")
        meanvar <- binMeanVar(x, conc=NULL, object$samples$group, nbins=nbins)
    if(show.tagwise.vars) {
        if(dispersion.method=="coxreid" & is.null(object$CR.tagwise.dispersion))
            stop("Cannot extract Cox-Reid tagwise dispersions. Try running CRDisp on your object first.")
        if(dispersion.method=="qcml" & is.null(object$tagwise.dispersion))
            stop("Cannot extract tagwise dispersions. Try running estimateTagwiseDisp() on your object first.")
        tagvars <- meanvar$means + meanvar$means^2*tagwise.dispersion
    }
    if(show.raw.vars) {
        plot(meanvar$means, meanvar$vars, log="xy", col="gray60", cex=0.6, xlab="Mean gene expression level (log10 scale)", ylab="Pooled gene-level variance (log10 scale)",plot.first=grid(), ...)
        if(show.tagwise.vars)
            points(meanvar$means, tagvars, col="lightskyblue", cex=0.6)
        points(meanvar$avemeans, meanvar$avevars, pch="x", col="red", cex=1.5)
    }
    else {
        if(show.tagwise.vars) {
            plot(meanvar$means, tagvars, col="lightskyblue", log="xy", cex=0.6, xlab="Mean gene expression level (log10 scale)", ylab="Estimated tagwise variance (log10 scale)",plot.first=grid(), ...)
            points(meanvar$avemeans, meanvar$avevars, pch="x", col="red", cex=1.5)
        }
        else {
            if( any(!is.finite(meanvar$avevars)) )
                maxy <- max(meanvar$vars)
            else
                maxy <- max(meanvar$avevars)
            plot(meanvar$avemeans, meanvar$avevars, pch="x", col="red", cex=1.5, ylim=c(0.1,maxy), log="xy", xlab="Mean gene expression level (log10 scale)", ylab="Pooled gene-level variance (log10 scale)",plot.first=grid(), ...)
        }
    }
    abline(0,1,lwd=2)
    if(NBline)
        lines(10^lmu,nb.var,col="dodgerblue3",lwd=3)
    return(invisible(meanvar))
}


