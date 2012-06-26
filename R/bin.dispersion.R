### BIN.DISPERSION.R

binCMLDispersion <- function(y, nbins=50) {
    ## Function to bin DGE data based on abundance and calculate the common dispersion for each bin. Allows us to investigate any dependence of the dispersion on the abundance.
    ## Written by Davis McCarthy.
    ## Oct 2010. Last modified 10 Feb 2011.
    if( !is(y, "DGEList") )
        stop("Function only operates on DGEList objects.\n")
    if( is.null( y$common.dispersion ) ) {
        y <- estimateCommonDisp(y)
        cat("Running estimateCommonDisp() on object to get concentration/abundance.\n")
    }
    ntags <- nrow(y)

    disp.bins <- logCPM.bins <- rep(NA,nbins)
    o <- order(y$logCPM)
    ntagsinbin <- floor(ntags / nbins)

    for(i in 1:nbins) {
        if( i==nbins )
            bin <- o[ (1 + (i-1)*ntagsinbin):ntags]
        else
            bin <- o[ (1 + (i-1)*ntagsinbin):( i*ntagsinbin)]

        disp.bins[i] <- estimateCommonDisp(y[bin,])$common.dispersion
        logCPM.bins[i] <- mean(y$logCPM[bin])
    }
    list(dispersion.bins=disp.bins, logCPM.bins=logCPM.bins)
}


binGLMDispersion <- function( y, design, min.n=500, offset=NULL, method="CoxReid", ... )
    ## Bin DGE data based on abundance and compute the Cox-Reid estimate of the common dispersion in each bin
    ## Written by Davis McCarthy.
    ## 7 Feb 2011. Last modified 13 May 2011.
{
    if( is(y, "DGEList") ) {
        if(is.null(offset))
           offset <- getOffset(y)
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
    y.nonzero <- y[!all.zero,]
    abundance <- mglmOneGroup(y,offset=offset)
#    o <- order(abundance)
#    ntagsinbin <- floor(ntags / nbins)
 
    ## Define bins of genes; based on min.n in each bin
    nbins <- floor(length(unique(abundance))/min.n)  ## This allows for cutting to be properly determined (no ties on abundance)
    if( nbins < 2 ) {
        nbins <- 2
        warning("With ",nrow(y)," genes and setting the parameter minimum number (min.n) of genes per bin to ",min.n,",  there should technically be fewer than 2 bins. To make estimation of trended dispersions possible we set the number of bins to be 2.\n")
    }
    if( nbins < 8) {
        min.n.new <- floor(length(unique(abundance))/nbins) 
        warning("With ",nrow(y)," genes and setting the parameter minimum number (min.n) of genes per bin to ",min.n,",  there are only ", nbins, " bins. Using ",nbins," bins here means that the minimum number of genes in each of the ", nbins, " bins is in fact ", min.n.new,". This number of bins and minimum number of genes per bin may not be sufficient for reliable estimation of a trend on the dispersions.\n")
        min.n <- min.n.new
    }
    bins <- cutWithMinN(abundance[!all.zero], intervals=nbins, min.n=min.n)
    dispersion <- ave.abundance <- rep(NA,nbins)
   
    for(i in 1:nbins) {
        ##if( i==nbins )
        ##    bin <- o[ (1 + (i-1)*ntagsinbin):ntags]
        ##else
        ##    bin <- o[ (1 + (i-1)*ntagsinbin):( i*ntagsinbin)]
        bin <- bins$group==i
        dispersion[i] <- estimateGLMCommonDisp(y[bin,], design, method=method, offset[bin,], min.row.sum=0, ...)
        ave.abundance[i] <- mean(abundance[bin])
    }
    if( any(is.na(all.zero)) ) {
        warning("Some tags/genes have zero counts in all libraries. These genes were ignored when computing binned dispersions.\n")
#        keep <- !is.na(dispersion)
#        dispersion <- dispersion[keep]
#        ave.abundance <- ave.abundance[keep]
    }
    
    new("list", list(dispersion=dispersion, abundance=ave.abundance))
}





