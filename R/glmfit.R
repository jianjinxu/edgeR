#  FIT GENERALIZED LINEAR MODELS

glmFit <- function(y, design, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL, prior.count=0.125, start=NULL, method="auto", ...)
UseMethod("glmFit")

glmFit.DGEList <- function(y, design=NULL, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL, prior.count=0.125, start=NULL, method="auto", ...)
#	Created 11 May 2011.  Last modified 11 March 2013.
{
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("No dispersion values found in DGEList object.")
	if(is.null(offset) && is.null(lib.size)) offset <- getOffset(y)
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	fit <- glmFit(y=y$counts,design=design,dispersion=dispersion,offset=offset,weights=weights,lib.size=lib.size,prior.count=prior.count,start=start,method=method,...)
	fit$samples <- y$samples
	fit$genes <- y$genes
	fit$prior.df <- y$prior.df
	fit$AveLogCPM <- y$AveLogCPM
	new("DGEGLM",fit)
}

glmFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, weights=NULL, lib.size=NULL, prior.count=0.125, start=NULL, method="auto", ...)
#	Fit negative binomial generalized linear model for each transcript
#	to a series of digital expression libraries
#	Davis McCarthy and Gordon Smyth
#	Created 17 August 2010. Last modified 13 Nov 2012.
{
#	Check input
	y <- as.matrix(y)
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
		ne <- nonEstimable(design)
		if(!is.null(ne)) stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n", paste(ne, collapse = " ")))
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
				method <- "levenberg"
			}
		}
	}
	if(method!="simple") {
		if(iswt) stop("weights only supported by simple fitting method")
		if(isna) stop("NAs only supported by simple fitting method")
	}

#	Fit a glm to each gene
	fit <- switch(method,
		linesearch=mglmLS(y,design=design,dispersion=dispersion,coef.start=start,offset=offset,...),
		oneway=mglmOneWay(y,design=design,dispersion=dispersion,offset=offset),
		levenberg=mglmLevenberg(y,design=design,dispersion=dispersion,offset=offset,coef.start=start,maxit=250,...),
		simple=mglmSimple(y,design=design,dispersion=dispersion,offset=offset,weights=weights)
	)

#	Prepare output
	fit$counts <- y
	if(prior.count>0)
		fit$coefficients <- predFC(y,design,offset=offset,dispersion=dispersion,prior.count=prior.count)*log(2)
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
#	if(is.null(fit$abundance)) fit$abundance <- mglmOneGroup(y, offset=offset, dispersion=dispersion)
	if(is.null(fit$design)) fit$design <- design
	if(is.null(fit$offset)) fit$offset <- offset
	if(is.null(fit$dispersion)) fit$dispersion <- dispersion
	fit$method <- method
	new("DGEGLM",fit)
}


glmLRT <- function(glmfit,coef=ncol(glmfit$design),contrast=NULL,test="chisq")
#	Tagwise likelihood ratio tests for DGEGLM
#	Gordon Smyth, Davis McCarthy and Yunshun Chen.
#	Created 1 July 2010.  Last modified 14 Dec 2012.
{
#	Check glmfit
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
	}
	nlibs <- ncol(glmfit)

#	Check test
	test <- match.arg(test,c("F","f","chisq"))
	if(test=="f") test <- "F"
	
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

#	Chisquare or F-test	
	LRT.pvalue <- switch(test,
		"F" = {
			phi <- quantile(glmfit$dispersion,p=0.5)
			mu <- quantile(glmfit$fitted.values,p=0.5)
			gamma.prop <- (phi*mu/(1 + phi*mu))^2
			prior.df <- glmfit$prior.df
			if(is.null(prior.df)) prior.df <- 20
			glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
			pf(LR/df.test, df1=df.test, df2=glmfit$df.total, lower.tail = FALSE, log.p = FALSE)
		},
		"chisq" = pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)
	)

	rn <- rownames(glmfit)
	if(is.null(rn))
		rn <- 1:nrow(glmfit)
	else
		rn <- make.unique(rn)
	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
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

