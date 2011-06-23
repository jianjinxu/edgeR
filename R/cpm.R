cpm <- function(x, normalized.lib.sizes=FALSE) {
    ## Calculate counts per million from a DGEList object or a matrix
    ## Written by Davis McCarthy and Gordon Smyth.
    ## Created 20 June 2011. Last modified 20 June 2011
    if(is(x, "DGEList")) {
        if(normalized.lib.sizes)
            lib.size <- x$samples$lib.size*x$samples$norm.factors
        else
            lib.size <- x$samples$lib.size
        x <- x$counts
    }
    else {
        if(normalized.lib.sizes)
            warning("Matrix of counts supplied, so normalized library sizes are not known. Library sizes are column sums of the count matrix.\n")
        lib.size <- colSums(x)
    }
    cpm <- 1e06*t(t(x)/lib.size)
    cpm
}


