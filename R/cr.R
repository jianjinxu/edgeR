CRDisp <- function(y, design=NULL, offset=0, npts=10, min.disp=1e-2, max.disp=2, nrandom=200, rowsum.filter = 5, prior.n=10, lib.size=NULL, verbose=TRUE)
    ## Function to estimate the common dispersion and tagwise disperisons using CoxReid Adjusted Profile-likelihood
    ## The function uses cubic spline interpolation in finding the MLEs.
    ## Written by Yunshun Chen, August 2010. Last modified by Yunshun Chen, 10 Sept 2010

{
    if( is(y,"DGEList") ) {
        if(is.null(lib.size))
            lib.size <- y$samples$lib.size*y$samples$norm.factors
        y.mat <- y$counts
    } else {
        y.mat <- as.matrix(y)
        if(is.null(lib.size)) {
            lib.size <- colSums(y.mat)
            if(verbose)
                cat("No lib.size supplied, so lib.size is taken as the column sums of the matrix of counts.")
        }
    }
    if(is.null(design)) {
        if( is(y, "DGEList") ) {
           design <- model.matrix(~y$samples$group)
           if(verbose)
               cat("Design matrix is being formed from the DGEList object, using y$samples$group.\n")
       }
       else
           stop("No design matrix supplied as an argument with matrix of counts.")
    }
    ngenes <- nrow(y.mat)
    narrays <- ncol(y.mat)
    if( length(offset)==length(y.mat) ) {
        offset.mat <- as.matrix(offset, nrow=ngenes, ncol=narrays)
    } else {
        offset.mat <- matrix(0, nrow=ngenes, ncol=narrays)
        if(length(offset)==narrays | length(offset)==1)
            offset.mat <- matrix(offset, nrow=ngenes, ncol=narrays, byrow=TRUE)
        else 
            stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
    }
        lib.size.mat <- outer(rep(1,nrow(y.mat)), log(lib.size))
        tags.used <- rowSums(y.mat) > rowsum.filter
        y.filt <- y.mat[tags.used,]
        offset.mat.filt <- offset.mat[tags.used,]
        ntags <- nrow(y.filt)

        tags.random <- sample(1:nrow(y.filt), nrandom)
        offset.mat.random <- offset.mat.filt[tags.random,]
        y.random <- y.filt[tags.random,]
        lib.size.mat.random <- outer(rep(1,nrow(y.random)), log(lib.size))
        lib.size.mat.filt <- outer(rep(1,nrow(y.filt)), log(lib.size))  

        a1 <- min.disp^(0.25)
        b1 <- max.disp^(0.25)
        spline.pts <- a1 + (0:(20-1))*(b1-a1)/(20-1)
        spline.disp <- (spline.pts)^4
        apl.com.random <- c()
        apl.tgw.random <- matrix(0, nrow=20, ncol=nrandom)
        ## mustart <- NULL
        for(i in 1:20){
                ## y.apl <- adjustedProfileLik(spline.disp[i], y.random, design=design, offset=offset.mat.random+lib.size.mat.random, mustart=mustart)
                y.apl <- adjustedProfileLik(spline.disp[i], y.random, design=design, offset=offset.mat.random+lib.size.mat.random)
                ## apl.tgw.random[i,] <- y.apl[[1]]
                apl.tgw.random[i,] <- y.apl
                ## mustart <- y.apl[[2]]
        }
        apl.com.random <- rowSums(apl.tgw.random)/nrandom
        bound <- max(apl.com.random) - 5/prior.n
        select <- apl.com.random > bound
        min.disp.new <- min(spline.disp[select])
        max.disp.new <- max(spline.disp[select])

        a <- min.disp.new^(0.25)
        b <- max.disp.new^(0.25)
        spline.pts <- a + (0:(npts-1))*(b-a)/(npts-1)
        spline.disp <- (spline.pts)^4
        apl.com <- c()
        apl.tgw <- matrix(0, nrow=npts, ncol=ntags)
        ## mustart <- NULL
        for(i in 1:npts){
                ## y.apl <- adjustedProfileLik(spline.disp[i], y.filt, design=design, offset=offset.mat.filt+lib.size.mat.filt, mustart=mustart)
                y.apl <- adjustedProfileLik(spline.disp[i], y.filt, design=design, offset=offset.mat.filt+lib.size.mat.filt)
                ## apl.tgw[i,] <- y.apl[[1]]
                apl.tgw[i,] <- y.apl
                ## mustart <- y.apl[[2]]
        }
        apl.com <- rowSums(apl.tgw)/ntags
        cr.com <- (maximize.by.interpolation(spline.pts, apl.com))^4
        cr.tgw.all <- rep(0,ngenes)
        cr.tgw.filt <- rep(0,ntags)
        for(i in 1:ntags) cr.tgw.filt[i] <- (maximize.by.interpolation(spline.pts, apl.tgw[,i]+ prior.n*apl.com))^4
        cr.tgw.all[tags.used] <- cr.tgw.filt
        new("DGEList",list(samples=y$samples, counts=y$counts, genes=y$genes, CR.common.dispersion=cr.com,
                CR.tagwise.dispersion=cr.tgw.all))
}






############################################
####### *. Adjusted Profile likelihood #####
############################################

adjustedProfileLik <- function(phi, y, design, offset) {
    ## Function to calculate the adjusted profile-likelihood given dispersion, design matrix and response.
    ## Created by Andy Chen, June 2010. Last modified by Yunshun Chen, 3 Sept 2010.
    ## y is simply a table of counts: rows are genes/tags/transcripts, columns are samples/libraries
    ## offset needs to be a matrix of offsets of the same dimensions as y
    require(MASS)
    if(!identical(dim(y), dim(offset)))
        stop("offset must be a matrix with the same dimensions as y, the table of counts.\n")
    tgw.apl <- c()
    ## mu <- matrix(0, nrow=nrow(y), ncol=ncol(y))
    theta <- 1/phi
    f <- negative.binomial(theta)
    for(i in 1:nrow(y)){
        ## fit <- glm.fit(design, y[i,], offset = offset[i,], mustart=mustart[i,], family = f)
          fit <- glm.fit(design, y[i,], offset = offset[i,], family = f)
        ## mu[i,] <- fitted(fit)
          loglik <- -fit$aic/2+fit$qr$rank
        cr <- sum(log(abs(diag(fit$qr$qr)[1:fit$qr$rank])))
          tgw.apl[i] <- loglik - cr
    }
    ## return(list(tgw.apl, mu))
    return(tgw.apl)
}




maximize.by.interpolation <- function(x,z,maxit=10,eps=1e-7,plot=FALSE)
#       Maximize a function given a table of values
#       by spline interpolation
#       Gordon Smyth
#       26 August 2010. Modified 1 Sept 2010.
{
        n <- length(z)
        imax <- which.max(z)
        r <- range(x)
        x0 <- x[imax]

#       If maximum occurs at end point, return that value
        if(x0==r[1] || x0==r[2]) return(x0)

        f <- splinefun(x,z)
        if(plot) {
                xx <- seq(from=r[1],to=r[2],length=100)
                zz <- f(xx)
                plot(xx,zz,type="l")
                points(x,z)
        }
        x <- x0
        for (iter in 1:maxit) {
                step <- f(x,deriv=1)/f(x,deriv=2)
                x <- x-step
                if(x<r[1] || x>r[2]) {
                        warning("Divergence")
                        return(x0)
                }
                if(abs(step) < eps) return(x)
        }
        warning("max iterations exceeded")
        x
}
