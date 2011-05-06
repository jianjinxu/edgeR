### BIN.DISPERSION.R

binCMLDispersion <- function(y, nbins=50) {
    ## Function to bin DGE data based on abundance and calculate the common dispersion for each bin. Allows us to investigate any dependence of the dispersion on the abundance.
    ## Written by Davis McCarthy.
    ## Oct 2010. Last modified 10 Feb 2011.
    if( !is(y, "DGEList") )
        stop("Function only operates on DGEList objects.\n")
    if( is.null( y$conc$conc.common ) ) {
        y <- estimateCommonDisp(y)
        cat("Running estimateCommonDisp() on object to get concentration/abundance.\n")
    }
    ntags <- nrow(y)
    abundance <- log(y$conc$conc.common)

    disp.bins <- abundance.bins <- rep(NA,nbins)
    o <- order(abundance)
    ntagsinbin <- floor(ntags / nbins)

    for(i in 1:nbins) {
        if( i==nbins )
            bin <- o[ (1 + (i-1)*ntagsinbin):ntags]
        else
            bin <- o[ (1 + (i-1)*ntagsinbin):( i*ntagsinbin)]

        disp.bins[i] <- estimateCommonDisp(y[bin,])$common.dispersion
        abundance.bins[i] <- mean(abundance[bin])
    }
    list(dispersion=disp.bins, abundance=abundance.bins)
}


binGLMDispersion <- function( y, design, nbins=50, offset=NULL, method="CoxReid", ... )
    ## Bin DGE data based on abundance and compute the Cox-Reid estimate of the common dispersion in each bin
    ## Written by Davis McCarthy.
    ## 7 Feb 2011. Last modified 3 May 2011.
{
    if( is(y, "DGEList") ) {
        if(is.null(offset))
           offset <- getOffsets(y)
        y <- y$counts
    }
    else {
        y <- as.matrix(y)
        lib.size <- colSums(y)
        if( is.null(offset) )
            offset <- log(lib.size)
    }
    nlibs <- ncol(y)
    ntags <- nrow(y)
    offset <- expandAsMatrix(offset,dim(y))

    method <- match.arg(method, c("CoxReid", "Pearson", "deviance"))
    all.zero <- rowSums(y)==0
    y[all.zero,1] <- 1
    abundance <- mglmOneGroup(y,offset=offset)
    o <- order(abundance)
    ntagsinbin <- floor(ntags / nbins)
    dispersion <- ave.abundance <- rep(NA,nbins)
    
    for(i in 1:nbins) {
        if( i==nbins )
            bin <- o[ (1 + (i-1)*ntagsinbin):ntags]
        else
            bin <- o[ (1 + (i-1)*ntagsinbin):( i*ntagsinbin)]
        if( any(!all.zero[bin]) ) {
            dispersion[i] <- estimateGLMCommonDisp(y[bin,], design, method=method, offset[bin,], min.row.sum=0, ...)
            ave.abundance[i] <- mean(abundance[bin])
        }
    }
    if( any(is.na(dispersion)) ) {
        warning("Some bins contained only tags/genes with zero counts in all libraries. These bins were ignored.\n")
        keep <- !is.na(dispersion)
        dispersion <- dispersion[keep]
        ave.abundance <- ave.abundance[keep]
    }
    
    new("list", list(dispersion=dispersion, abundance=ave.abundance))
}





