#  FIT GENERALIZED LINEAR MODELS

glmFit <- function(y, design, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL, prior.count.total=0.5, start=NULL, method="auto", ...)
UseMethod("glmFit")

glmFit.DGEList <- function(y, design=NULL, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL, prior.count.total=0.5, start=NULL, method="auto", ...)
{
	if( is.null(dispersion) ) {
		if( !is.null(y$tagwise.dispersion) )
			dispersion <- y$tagwise.dispersion
		else {
			if( !is.null(y$trended.dispersion) )
				dispersion <- y$trended.dispersion
			else {
				if( !is.null(y$common.dispersion) )
					dispersion <- y$common.dispersion
				else
					stop("No dispersion values found in DGEList object. Run dispersion estimation functions such as estimateGLMCommonDisp, estimateGLMTrendedDisp and estimateGLMTagwiseDisp before using glmFit.\n")
			}
		}
	}
	if(is.null(offset) && is.null(lib.size)) offset <- getOffset(y)
	fit <- glmFit(y=y$counts,design=design,dispersion=dispersion,offset=offset,weights=weights,lib.size=lib.size,prior.count.total=prior.count.total,start=start,method=method,...)
	fit$samples <- y$samples
	fit$genes <- y$genes
	new("DGEGLM",fit)
}

glmFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL, prior.count.total=0.5, start=NULL, method="auto", ...)
#	Fit negative binomial generalized linear model for each transcript
#	to a series of digital expression libraries
#	Davis McCarthy and Gordon Smyth

#	Created 17 August 2010. Last modified 1 May 2012.
{
#	Check input
	y <- as.matrix(y)
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(is.null(dispersion)) {
		stop("No dispersion values provided.")
	} else {
		if(!( length(dispersion)==1 | length(dispersion)==nrow(y) ))
			stop("Length of dispersion vector incompatible with count matrix. Dispersion argument must be either of length 1 (i.e. common dispersion) or length equal to the number of rows of y (i.e. individual dispersion value for each tag/gene).")
	}
	if(!is.null(offset) && !is.null(lib.size)) warning("offset and lib.size both supplied: offset takes precedence, lib.size ignored.")
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(is.null(offset)) offset <- log(lib.size)
	offset <- expandAsMatrix(offset,dim(y))
	iswt <- !is.null(weights)
	if(iswt) {
		weights <- expandAsMatrix(weights,dim(y))
		weights[weights <= 0] <- NA
		y[!is.finite(weights)] <- NA
	}
	method <- match.arg(method,c("auto","linesearch","levenberg","simple"))
#	End of input checking

	ngenes <- nrow(y)
	nlibs <- ncol(y)
	isna <- any(is.na(y))

#	Choose fitting algorithm
	if(method=="auto") {
		if(isna || iswt) {
			method <- "simple"
		} else {
			group <- designAsFactor(design)
			if(nlevels(group)==ncol(design)) {
				method <- "oneway"
			} else {
				method <- "linesearch"
			}
		}
	}
	if(method!="simple") {
		if(iswt) stop("weights only supported by simple fitting method")
		if(isna) stop("NAs only supported by simple fitting method")
	}

#	Fit a glm to each gene
	fit <- switch(method,
		linesearch=mglmLS(y,design=design,dispersion=dispersion,start=start,offset=offset,...),
		oneway=mglmOneWay(y,design=design,dispersion=dispersion,offset=offset),
		levenberg=mglmLevenberg(y,design=design,dispersion=dispersion,offset=offset),
		simple=mglmSimple(y,design=design,dispersion=dispersion,offset=offset,weights=weights)
	)

#	Prepare output
	fit$counts <- y
	if(prior.count.total>0)
		fit$coefficients <- predFC(y,design,offset=offset,dispersion=dispersion,prior.count.total=prior.count.total)
	else
		fit$coefficients <- as.matrix(fit$coefficients)
	colnames(fit$coefficients) <- colnames(design)
	rownames(fit$coefficients) <- rownames(y)
	fit$fitted.values <- as.matrix(fit$fitted.values)
	dimnames(fit$fitted.values) <- dimnames(y)
	if(is.null(fit$deviance)) {
		deviances <- deviances.function(dispersion)
		fit$deviance <- deviances(y,fit$fitted.values,dispersion)
	}
	if(is.null(fit$df.residual)) fit$df.residual <- rep(nlibs-ncol(design),ngenes)
	if(is.null(fit$abundance)) fit$abundance <- mglmOneGroup(y, offset=offset, dispersion=dispersion)
	if(is.null(fit$design)) fit$design <- design
	if(is.null(fit$offset)) fit$offset <- offset
	if(is.null(fit$dispersion)) fit$dispersion <- dispersion
	fit$method <- method
	new("DGEGLM",fit)
}


glmLRT <- function(glmfit,coef=ncol(glmfit$design),contrast=NULL)
#	Tagwise likelihood ratio tests for DGEGLM
#	Gordon Smyth and Davis McCarthy.
#	Created 1 July 2010. Last modified 3 July 2012.
{
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef or contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
	}

#	Full design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

#	contrast takes precedence over coef
#	Evaluate contrast
#	Reform design matrix so that contrast of interest is last column
	if(is.null(contrast)) {
        if(length(coef) > 1)
            coef <- unique(coef)
		if(is.character(coef)) {
            check.coef <- coef %in% colnames(design)
            if( any(!check.coef) )
                stop("One or more named coef arguments do not match a column of the design matrix.\n")
			coef.name <- coef
            coef <- match(coef, colnames(design))
        }
        else
            coef.name <- coef.names[coef]
        logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
        if(length(coef)==1) logFC <- as.vector(logFC)
	} else {
		logFC <- (glmfit$coefficients %*% contrast)/log(2)
		i <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		qr <- qr(contrast)
		Q <- qr.Q(qr,complete=TRUE)
		Q <- cbind(Q[,-1],Q[,1]*qr$qr[1,1])
		design <- design %*% Q
		coef <- nbeta
	}

#	Null design matrix
	design0 <- design[,-coef,drop=FALSE]

#	Null fit
	fit.null <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion)

	LR <- fit.null$deviance - glmfit$deviance
	df.test <- fit.null$df.residual - glmfit$df.residual
	LRT.pvalue <- pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)
	rn <- rownames(glmfit)
	if(is.null(rn))
		rn <- 1:nrow(glmfit)
	else
		rn <- make.unique(rn)
	tab <- data.frame(
		logFC=logFC,
		logCPM=(glmfit$abundance+log(1e6))/log(2),
		LR=LR,
		PValue=LRT.pvalue,
		row.names=rn
	)
	glmfit$counts <- NULL
	glmfit$table <- tab 
	glmfit$comparison <- coef.name
	glmfit$df.test <- df.test
	new("DGELRT",unclass(glmfit))
}

