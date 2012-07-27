#############################################################
# maximizeInterpolant: written by Aaron Lun
#
# This function does the same thing as maximizeInterpolant
# except that compiled C++ code lies at its core. It takes an
# ordered set of spline points and a likelihood matrix where
# each row corresponds to a tag and each column corresponds
# to a spline point. It then calculates the position at 
# which the maximum interpolated likelihood occurs for each
# row. It is faster than the original and also calculates
# exact solutions (avoiding errors in the Newton-Raphson
# algorithm where narrow peaks are 'jumped' over).

maximizeInterpolant <- function( x, y ) {
    if (is.vector(y)) {
        y<-rbind(y);
        warning("Converting vector of likelihoods to matrix format.");
    }
    if (length(x)!=ncol(y)) { 
        stop("Number of columns must equal number of spline points.");
    } else if (is.unsorted(x) || anyDuplicated(x)) {
        stop("Spline points must be unique and sorted.");
    }
    out<-.Call("maximize_interpolant", x, y, PACKAGE="edgeR");
    return(out);
}
