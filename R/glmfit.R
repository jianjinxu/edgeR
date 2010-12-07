#  GENERALIZED LINEAR MODELS

glmFit <- function(y, design, dispersion, offset=NULL, weights=NULL, lib.size=NULL)
	##	Fit negative binomial generalized linear model for each transcript
	##  to a series of digital expression libraries
	##	Davis McCarthy and Gordon Smyth
	##	Created 17 August 2010. Last modified 24 Nov 2010.
	##  User-level function. Takes a matrix of counts or a DGElist object (y).
{
	if(is(y,"DGEList")) {
		y.mat <- y$counts
		samples <- y$samples
		genes <- y$genes
        if( !is.null(offset) & !is.null(lib.size))
            warning("If non-null offset and non-null lib.size are provided, then the offset argument will be used and the supplied lib.size ignored.\n")
		if(is.null(lib.size)) lib.size <- y$samples$lib.size
		lib.size <- lib.size*y$samples$norm.factors
	} else {
		y.mat <- as.matrix(y)
		samples <- genes <- NULL
		if(is.null(lib.size)) lib.size <- colSums(y.mat)
	}
	if(any(is.na(y.mat))) stop("Not currently supporting NAs in y")
	nlibs <- ncol(y.mat)
	design <- as.matrix(design)
	if(!is.null(weights)) {
		warning("weights not currently supported")
		weights <- asMatrixWeights(weights,dim(y.mat))
		weights[weights <= 0] <- NA
		y.mat[!is.finite(weights)] <- NA
	}
	ngenes <- nrow(y.mat)
										# Define objects in which to store various results from the glm fits
										#  Check length of dispersion vector and that dispersion is not equal to zero
	if( !( length(dispersion)==1 | length(dispersion)==nrow(y.mat) ) )
		stop("Length of 'dispersion' is not compatible with size of 'y'. Dispersion must have length either equal to one or equal to the number of rows of y.\n")
	if( length(dispersion==1) )
		dispersion <- rep(dispersion, length=nrow(y.mat))
	if( !is.null(offset) & length(offset)!=nlibs & length(offset)!=1 & length(offset)!=length(y.mat) )
		stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
	else {
        if( !is.null(offset) )
            if( length(offset)==length(y.mat) )
                offset <- matrix(offset, nrow=ngenes, ncol=nlibs)
            else
                offset <- matrix(offset, nrow=ngenes, ncol=nlibs, byrow=TRUE)
        else
            offset <- matrix(log(lib.size), nrow=ngenes, ncol=nlibs, byrow=TRUE)
    }
                                        # Fit a glm to each gene sequentially using the line-search algorithm implemented in mglmLS()
	fit <- mglmLS(y.mat, design=design, dispersion=dispersion, start=NULL, offset=offset, tol=1e-5, maxit=50, trace=FALSE)
	coefficients <- fit$coefficients
	fitted.values <- fit$fitted
	colnames(coefficients) <- colnames(design)
	rownames(coefficients) <- rownames(y.mat)
	dimnames(fitted.values) <- dimnames(y.mat)
										# Compute deviances
	deviances <- deviances.function(dispersion)
	dev <- deviances(y.mat,fitted.values,dispersion)
										# Compute residual degrees of freedom
	df.residual <- rep(nlibs-ncol(design),ngenes)
										# What do we do about weights with this new function?
##	abundance <- log2( rowMeans(t( t(y.mat)/lib.size )) )
	abundance <- mglmOneGroup(y.mat, offset=offset, dispersion=dispersion, maxit=50) # Is this the best way to get the abundance?
	new("DGEGLM",list(coefficients=coefficients, df.residual=df.residual, deviance=dev, design=design, 
					 offset=offset, samples=samples, genes=genes, dispersion=dispersion, 
					 lib.size=lib.size, weights=weights, fitted.values=fitted.values, abundance=abundance))

}


glmLRT <- function(y,glmfit,coef=ncol(glmfit$design),contrast=NULL)
#	Tagwise likelihood ratio tests for DGEGLM
#	Gordon Smyth and Davis McCarthy.
#	Created 1 July 2010. Last modified 24 November 2010.
{
	if(is(y,"DGEList"))
		y.mat <- y$counts 
	else
		y.mat <- as.matrix(y)
	if(!is(glmfit,"DGEGLM"))
		stop("The glmfit argument must be a DGEGLM object for the full model. Run glmFit with the design matrix of the full model before LR testing.\n")

#	Full design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

#	contrast takes precedence over coef
#	Evaluate contrast
#	Reform design matrix so that contrast of interest is last column
	if(is.null(contrast)) {
		logFC <- glmfit$coefficients[,coef]/log(2)
		if(is.character(coef))
			coef.name <- coef
		else
			coef.name <- coef.names[coef]
	} else {
		logFC <- (glmfit$coefficients %*% contrast)/log(2)
		i <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		qr <- qr(contrast)
		Q <- qr.Q(qr,complete=TRUE)
		sign1 <- sign(qr$qr[1,1])
		Q <- cbind(Q[,-1],Q[,1])
		design <- design %*% Q
		coef <- nbeta
	}

#	Null design matrix
	design0 <- design[,-coef,drop=FALSE]

#	Null fit
	fit.null <- glmFit(y,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,lib.size=NULL)

	LR <- fit.null$deviance - glmfit$deviance
	LRT.pvalue <- pchisq(LR, df=( fit.null$df.residual - glmfit$df.residual ), lower.tail = FALSE, log.p = FALSE)
	tab <- data.frame(
		logConc=glmfit$abundance,
		logFC=logFC,
		LR.statistic=LR,
		p.value=LRT.pvalue
	)
	rownames(tab) <- rownames(y.mat)
	if(is(y,"DGEList")) {
		y$counts <- NULL
		y$pseudo.alt <- NULL
		y$table <- tab 
		y$coefficients.full <- glmfit$coefficients
		y$coefficients.null <- fit.null$coefficients
		y$design.full <- glmfit$design
		y$design.null <- design0
		y$dispersion.used <- glmfit$dispersion
	} else {
		y <- list(table=tab, coefficients.full=glmfit$coefficients, coefficients.null=fit.null$coefficients, design.full=glmfit$design, dispersion.used=glmfit$dispersion)
	}
	y$comparison <- coef.name
	new("DGELRT",unclass(y))
}

