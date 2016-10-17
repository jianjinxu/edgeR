#  FIT GENERALIZED LINEAR MODELS

glmFit <- function(y, ...)
UseMethod("glmFit")

glmFit.DGEList <- function(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, ...)
#	Created 11 May 2011.  Last modified 14 Aug 2016.
{
	if(is.null(design)) design <- model.matrix(~y$samples$group)
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("No dispersion values found in DGEList object.")

	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)

	fit <- glmFit(y=y$counts,design=design,dispersion=dispersion,offset=getOffset(y),lib.size=NULL,weights=y$weights,prior.count=prior.count,start=start,...)
	fit$samples <- y$samples
	fit$genes <- y$genes
	fit$prior.df <- y$prior.df
	fit$AveLogCPM <- y$AveLogCPM
	new("DGEGLM",fit)
}

glmFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, weights=NULL, prior.count=0.125, start=NULL, ...)
#	Fit negative binomial generalized linear model for each transcript
#	to a series of digital expression libraries
#	Davis McCarthy and Gordon Smyth
#	Created 17 August 2010. Last modified 03 Oct 2016.
{
#	Check y
	y <- as.matrix(y)
	ntag <- nrow(y)
	nlib <- ncol(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
		ne <- nonEstimable(design)
		if(!is.null(ne)) stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n", paste(ne, collapse = " ")))
	}

#	Check dispersion
	if(is.null(dispersion)) stop("No dispersion values provided.")
	dispersion.mat <- makeCompressedMatrix(dispersion, byrow=FALSE)

#	Check offset and lib.size
	offset <- .compressOffsets(y=y, lib.size=lib.size, offset=offset)

#	weights are checked in lower-level functions

#	Fit the tagwise glms
#	If the design is equivalent to a oneway layout, use a shortcut algorithm
	group <- designAsFactor(design)
	if(nlevels(group)==ncol(design)) {
		fit <- mglmOneWay(y,design=design,dispersion=dispersion.mat,offset=offset,weights=weights,coef.start=start)
		fit$deviance <- nbinomDeviance(y=y,mean=fit$fitted.values,dispersion=dispersion.mat,weights=weights)
		fit$method <- "oneway"
	} else {
		fit <- mglmLevenberg(y,design=design,dispersion=dispersion.mat,offset=offset,weights=weights,coef.start=start,maxit=250)
		fit$method <- "levenberg"
	}

#	Prepare output
	fit$counts <- y
	if(prior.count>0) {
		fit$unshrunk.coefficients <- fit$coefficients
		colnames(fit$unshrunk.coefficients) <- colnames(design)
		rownames(fit$unshrunk.coefficients) <- rownames(y)
		fit$coefficients <- predFC(y,design,offset=offset,dispersion=dispersion.mat,prior.count=prior.count,weights=weights,...)*log(2)
	}
	colnames(fit$coefficients) <- colnames(design)
	rownames(fit$coefficients) <- rownames(y)
	dimnames(fit$fitted.values) <- dimnames(y)
#	FIXME: we are not allowing missing values, so df.residual must be same for all tags
	fit$df.residual <- rep(nlib-ncol(design),ntag)
	fit$design <- design
	fit$offset <- offset
	fit$dispersion <- dispersion
	fit$weights <- weights
	fit$prior.count <- prior.count
	new("DGEGLM",fit)
}


glmLRT <- function(glmfit,coef=ncol(glmfit$design),contrast=NULL)
#	Tagwise likelihood ratio tests for DGEGLM
#	Gordon Smyth, Davis McCarthy and Yunshun Chen.
#	Created 1 July 2010.  Last modified 11 June 2015.
{
#	Check glmfit
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
	}
	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
	nlibs <- ncol(glmfit)
	
#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

#	Evaluate logFC for coef to be tested
#	Note that contrast takes precedence over coef: if contrast is given
#	then reform design matrix so that contrast of interest is last column.
	if(is.null(contrast)) {
		if(length(coef) > 1) coef <- unique(coef)
		if(is.character(coef)) {
			check.coef <- coef %in% colnames(design)
			if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
			coef.name <- coef
			coef <- match(coef, colnames(design))
		}
		else
			coef.name <- coef.names[coef]
		logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
	} else {
		contrast <- as.matrix(contrast)
		qrc <- qr(contrast)
		ncontrasts <- qrc$rank
		if(ncontrasts==0) stop("contrasts are all zero")
		coef <- 1:ncontrasts
		if(ncontrasts < ncol(contrast)) contrast <- contrast[,qrc$pivot[coef]]
		logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
		if(ncontrasts>1) {
			coef.name <- paste("LR test of",ncontrasts,"contrasts")
		} else {
			contrast <- drop(contrast)
			i <- contrast!=0
			coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		}
		Dvec <- rep.int(1,nlibs)
		Dvec[coef] <- diag(qrc$qr)[coef]
		Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
		design <- design %*% Q
	}
	if(length(coef)==1) logFC <- as.vector(logFC)

#	Null design matrix
	design0 <- design[,-coef,drop=FALSE]

#	Null fit
	fit.null <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,prior.count=0)

#	Likelihood ratio statistic
	LR <- fit.null$deviance - glmfit$deviance
	df.test <- fit.null$df.residual - glmfit$df.residual
	LRT.pvalue <-  pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)

	rn <- rownames(glmfit)
	if(is.null(rn))
		rn <- 1:nrow(glmfit)
	else
		rn <- make.unique(rn)
	tab <- data.frame(
		logFC=logFC,
		logCPM=glmfit$AveLogCPM,
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

