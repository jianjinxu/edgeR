treatDGE <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, lfc=0)
{
	message("treatDGE() has been renamed to glmTreat().  Please use the latter instead.")
	glmTreat(glmfit=glmfit,coef=coef,contrast=contrast,lfc=lfc)
}

glmTreat <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, lfc=0)
#	Likelihood ratio test or quasi-likelihood F-test with a threshold
#	Yunshun Chen and Gordon Smyth
#	Created on 05 May 2014. Last modified on 29 April 2015
{
	if(lfc < 0) stop("lfc has to be non-negative")
	
#	Check if glmfit is from glmFit() or glmQLFit()
	isLRT <- is.null(glmfit$df.prior)

#	Switch to glmLRT() or glmQLFTest() if lfc is zero
	if(lfc==0) {
		fun <- ifelse(isLRT, "glmLRT", "glmQLFTest")
		cat( paste0("Zero log2-FC threshold detected. Switch to ", fun, "() instead."), "\n" )
		return( do.call(fun, args=list(glmfit, coef, contrast)) )
	}

#	Check if there is any log-FC shrinkage
	shrunk <- glmfit$prior.count!=0
	
#	Check glmfit
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit or glmQLFit).")
	}
	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
	nlibs <- ncol(glmfit)
	ngenes <- nrow(glmfit)

#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

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
		unshrunk.logFC <- glmfit$coefficients[, coef, drop=FALSE]/log(2)
		if(shrunk) {
			logFC <- as.vector(unshrunk.logFC)
			unshrunk.logFC <- glmfit$unshrunk.coefficients[, coef, drop=FALSE]/log(2)
		}
	} else {
		contrast <- as.matrix(contrast)
		reform <- contrastAsCoef(design, contrast=contrast, first=TRUE)
		coef <- 1
		unshrunk.logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
		if(shrunk) {
			logFC <- as.vector(unshrunk.logFC)
			unshrunk.logFC <- drop((glmfit$unshrunk.coefficients %*% contrast)/log(2))
		}
		contrast <- drop(contrast)
		i <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		design <- reform$design
	}
	unshrunk.logFC <- as.vector(unshrunk.logFC)
	up <- unshrunk.logFC >= 0

#	Null design matrix
	design0 <- design[, -coef, drop=FALSE]

	if(isLRT) {
#		Adjusted offset
		offset.adj <- matrix(-lfc*log(2), ngenes, 1)
		offset.adj[up, ] <- lfc*log(2)
		
#		Test statistics at beta_0 = tau
		offset.new <- glmfit$offset + offset.adj %*% t(design[, coef, drop=FALSE])
		fit0.tau <- glmFit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
		fit1.tau <- glmFit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
		X2.tau <- pmax(0, fit0.tau$deviance - fit1.tau$deviance)
		z.tau <- sqrt(X2.tau)
		
#		Test statistics at beta_0 = -tau
		offset.new <- glmfit$offset - offset.adj %*% t(design[, coef, drop=FALSE])
		fit0.tau <- glmFit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
		fit1.tau <- glmFit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
		X2.tau2 <- pmax(0, fit0.tau$deviance - fit1.tau$deviance)
		z.tau2 <- sqrt(X2.tau2)

		within <- abs(unshrunk.logFC) <= lfc
		sgn <- 2*within - 1

#		Integral of Normal CDF from a to b
		fun <- function(a,b) ifelse(a==b, pnorm(a), ( b*pnorm(b)+dnorm(b) - (a*pnorm(a)+dnorm(a)) )/(b-a) )
		p.value <- 2*fun(-z.tau2, z.tau*sgn)
	} else {
#		Guass quadrature
		if(!requireNamespace("statmod",quietly=TRUE)) stop("statmod required but is not available")
		nnodes <- 5
		gquad <- statmod::gauss.quad.prob(nnodes, dist="uniform", l=0, u=lfc*log(2))
	
		offset.adj <- X2.pos <- X2.neg <- matrix(-gquad$nodes, ngenes, nnodes, byrow=TRUE)
		offset.adj[up, ] <- -offset.adj[up, ] 

		for(k in 1:nnodes){
#			Test statistics at beta_0 = pos_nodes (tau)
			offset.new <- glmfit$offset + offset.adj[,k] %*% t(design[, coef, drop=FALSE])
			fit0.pos <- glmQLFit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
			fit1.pos <- glmQLFit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
			X2.pos[,k] <- pmax(0, fit0.pos$deviance - fit1.pos$deviance)
			
#			Test statistics at beta_0 = neg_nodes (-tau)
			offset.new <- glmfit$offset - offset.adj[,k] %*% t(design[, coef, drop=FALSE])
			fit0.neg <- glmQLFit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
			fit1.neg <- glmQLFit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
			X2.neg[,k] <- pmax(0, fit0.neg$deviance - fit1.neg$deviance)
		}
		z.pos <- sqrt(X2.pos)
		z.neg <- sqrt(X2.neg)
	
		within <- matrix(abs(unshrunk.logFC)*log(2), ngenes, nnodes) <= abs(offset.adj)
		sgn <- 2*within - 1

#		Calculate expected p-values by Gauss quadrature
		t.pos <- z.pos/sqrt(glmfit$var.post)
		t.neg <- z.neg/sqrt(glmfit$var.post)
		df.total <- glmfit$df.prior + glmfit$df.residual.zeros
		max.df.residual <- ncol(glmfit$counts)-ncol(glmfit$design)
		df.total <- pmin(df.total, nrow(glmfit)*max.df.residual)
		p.value <- ( pt(t.pos*sgn, df=df.total) + pt(t.neg, df=df.total, lower.tail=FALSE) ) %*% gquad$weights
	
#		Ensure is not more significant than z-test
		i <- glmfit$var.post < 1
		if(any(i)) {
			z.pvalue <- ( pnorm(z.pos[i,]*sgn[i,]) + pnorm(z.neg[i,], lower.tail=FALSE) ) %*% gquad$weights
			p.value[i] <- pmax(p.value[i], z.pvalue)
		}	
	}

#	Table output
	tab <- data.frame(
		logCPM=glmfit$AveLogCPM,
		PValue=p.value,
		row.names=rownames(glmfit)
	)
	if(shrunk) {
		tab <- cbind(logFC=logFC, unshrunk.logFC=unshrunk.logFC, tab)
	} else {
		tab <- cbind(logFC=unshrunk.logFC, tab)
	}

	glmfit$lfc <- lfc
	glmfit$counts <- NULL
	glmfit$table <- tab
	glmfit$comparison <- coef.name
	new("DGELRT",unclass(glmfit))
}
