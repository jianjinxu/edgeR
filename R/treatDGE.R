treatDGE <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, lfc=0)
#	Likelihood ratio test with threshold
#	Yunshun Chen and Gordon Smyth
#	05 May 2014.
{
	if(lfc < 0) stop("lfc has to be non-negative")
#	Check glmfit
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
	}
	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
	nlibs <- ncol(glmfit)
	ngenes <- nrow(glmfit)

#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)
	
	if(glmfit$prior.count!=0){
		coefficients.mle <- glmfit$unshrunk.coefficients
	} else {
		coefficients.mle <- glmfit$coefficients
	}

#	Evaluate logFC for coef to be tested
#	Note that contrast takes precedence over coef: if contrast is given
#	then reform design matrix so that contrast of interest is the first column
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
		logFC <- coefficients.mle[, coef, drop=FALSE]/log(2)
	} else {
		contrast <- as.matrix(contrast)
		reform <- contrastAsCoef(design, contrast=contrast, first=TRUE)
		coef <- 1
		logFC <- drop((coefficients.mle %*% contrast)/log(2))
		contrast <- drop(contrast)
		i <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		design <- reform$design
	}
	logFC <- as.vector(logFC)

#	Null design matrix
	design0 <- design[, -coef, drop=FALSE]

#	LRT of beta_0 = zero
	fit0.zero <- glmFit(glmfit$counts, design=design0, offset=glmfit$offset, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	X2.zero <- pmax(0, fit0.zero$deviance - glmfit$deviance)

#	LRT of beta_0 = tau
	offset.adj <- matrix(-lfc*log(2), ngenes, 1)
	up <- logFC >= 0
	offset.adj[up, ] <- lfc*log(2)
	offset.new <- glmfit$offset + offset.adj %*% t(design[, coef, drop=FALSE])
	fit0.tau <- glmFit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	fit1.tau <- glmFit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	X2.tau <- pmax(0, fit0.tau$deviance - fit1.tau$deviance)
	
	z.zero <- sqrt(X2.zero)
	z.tau <- sqrt(X2.tau)
	within <- abs(logFC) <= lfc
	sgn <- 2*within - 1
	p.value <- pnorm( z.tau*sgn ) + pnorm( -z.tau*sgn-2*z.zero )

	tab <- data.frame(
		logFC=logFC,
		logCPM=glmfit$AveLogCPM,
		PValue=p.value,
		row.names=rownames(glmfit)
	)
	glmfit$lfc <- lfc
	glmfit$counts <- NULL
	glmfit$table <- tab 
	glmfit$comparison <- coef.name
	new("DGELRT",unclass(glmfit))
}



	
	
	
	