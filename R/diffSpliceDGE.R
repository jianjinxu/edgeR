diffSpliceDGE <- function(fit.exon, coef=ncol(fit.exon$design), geneid, exonid=NULL, verbose=TRUE)
{
# Identify exons and genes with splice variants using negative binomial GLMs
# Yunshun Chen and Gordon Smyth
# Created 29 March 2014.  Last modified 25 Aug 2014. 

#	Check input (from diffSplice in limma)
	exon.genes <- fit.exon$genes
	nexons <- nrow(fit.exon)
	design <- fit.exon$design

	if(is.null(exon.genes)) exon.genes <- data.frame(ExonID=1:nrow(fit.exon))
	if(length(geneid)==1) {
		genecolname <- as.character(geneid)
		geneid <- exon.genes[[genecolname]]
	} else {
		exon.genes$GeneID <- geneid
		genecolname <- "GeneID"
	}
	if(!is.null(exonid))
		if(length(exonid)==1) {
			exoncolname <- as.character(exonid)
			exonid <- exon.genes[[exoncolname]]
		} else {
			exon.genes$ExonID <- exonid
			exoncolname <- "ExonID"
		}
	else
		exoncolname <- NULL

#	Sort by geneid (from diffSplice in limma)
	if(is.null(exonid))
		o <- order(geneid)
	else
		o <- order(geneid,exonid)
	geneid <- geneid[o]
	exon.genes <- exon.genes[o,,drop=FALSE]

	fit.exon <- fit.exon[o, ]

#	Gene level information
	gene.counts <- rowsum(fit.exon$counts, geneid, reorder=FALSE)
	gene.dge <- DGEList(counts=gene.counts, genes=unique(geneid))
	gene.dge <- estimateDisp(gene.dge, design, robust=FALSE)
	fit.gene <- glmFit(gene.dge, design)

# 	Count exons and get genewise variances
	gene.nexons <- rowsum(rep(1,nexons), geneid, reorder=FALSE)
	if(verbose) {
		cat("Total number of exons: ", nexons, "\n")
		cat("Total number of genes: ", length(gene.nexons), "\n")
		cat("Number of genes with 1 exon: ", sum(gene.nexons==1), "\n")
		cat("Mean number of exons in a gene: ", round(mean(gene.nexons),0), "\n")
		cat("Max number of exons in a gene: ", max(gene.nexons), "\n")
	}

#	Squeeze
	fit.gene.trend <- glmFit(gene.dge, design=design, dispersion=gene.dge$trended.dispersion)
	zerofit <- (fit.gene.trend$fitted.values < 1e-4) & (fit.gene.trend$counts < 1e-4)
	gene.df.residual <- .residDF(zerofit, design)
	s2 <- fit.gene.trend$deviance / gene.df.residual
	s2[gene.df.residual==0] <- 0
	s2 <- pmax(s2,0)
	s2.fit <- squeezeVar(s2, df=gene.df.residual, covariate=fit.gene.trend$AveLogCPM, robust=FALSE)

#	Remove genes with only 1 exon
	gene.keep <- gene.nexons > 1
	ngenes <- sum(gene.keep)
	if(ngenes==0) stop("No genes with more than one exon")

	exon.keep <- rep(gene.keep, gene.nexons)
	geneid <- geneid[exon.keep]
	exon.genes <- exon.genes[exon.keep, , drop=FALSE]
	fit.exon <- fit.exon[exon.keep, ]

	fit.gene <- fit.gene[gene.keep, ]
	gene.nexons <- gene.nexons[gene.keep]
	gene.df.test <- gene.nexons-1
	gene.df.residual <- gene.df.residual[gene.keep]
	
# 	Genewise betas
	g <- rep(1:ngenes, times=gene.nexons)
	gene.counts.exon <- fit.gene$counts[g, , drop=FALSE]
	gene.dispersion.exon <- fit.gene$dispersion[g]
	gene.fit.exon <- glmFit(gene.counts.exon, design=design, dispersion=gene.dispersion.exon, lib.size=gene.dge$samples$lib.size)
	gene.betabar <- gene.fit.exon$coefficients[, coef, drop=FALSE]
	offset.new <- fit.exon$offset + gene.betabar %*% t(design[, coef, drop=FALSE])
	coefficients <- fit.exon$coefficients[, coef, drop=FALSE] - gene.betabar

#	Testing
	design0 <- design[, -coef, drop=FALSE]
	fit.null <- glmFit(fit.exon$counts, design=design0, offset=offset.new, dispersion=fit.exon$dispersion)
	fit.alt <- glmFit(fit.exon$counts, design=design, offset=offset.new, dispersion=fit.exon$dispersion)

# 	Exon p-values
	exon.LR <- fit.null$deviance - fit.alt$deviance
	exon.df.test <- fit.null$df.residual - fit.alt$df.residual	
	exon.F <- exon.LR / exon.df.test / s2.fit$var.post[gene.keep][g]
	gene.df.total <- s2.fit$df.prior + gene.df.residual
	max.df.residual <- ncol(fit.exon$counts)-ncol(design)
	gene.df.total <- pmin(gene.df.total, ngenes*max.df.residual)
	exon.p.value <- pf(exon.F, df1=exon.df.test, df2=gene.df.total[g], lower.tail=FALSE, log.p=FALSE)

	#Ensure is not more significant than chisquare test
	i <- s2.fit$var.post[gene.keep][g] < 1
	if(any(i)) {
		chisq.pvalue <- pchisq(exon.LR[i], df=exon.df.test[i], lower.tail=FALSE, log.p=FALSE)
		exon.p.value[i] <- pmax(exon.p.value[i], chisq.pvalue)
	}

#	Gene p-values

#	Gene Simes' p-values
	o <- order(g, exon.p.value, decreasing=FALSE)
	p <- exon.p.value[o]
	q <- rep(1, sum(gene.nexons))
	r <- cumsum(q) - rep(cumsum(q)[cumsum(gene.nexons)]-gene.nexons, gene.nexons)
	pp <- p*rep(gene.nexons, gene.nexons)/r
	#pp <- p*rep(gene.nexons-1, gene.nexons)/pmax(r-1, 1)
	oo <- order(-g, pmin(pp,1), decreasing=TRUE)
	gene.Simes.p.value <- pp[oo][cumsum(gene.nexons)]

#	Gene F p-values
	gene.F <- rowsum(exon.F, geneid, reorder=FALSE) / (gene.df.test)
	gene.F.p.value <- pf(gene.F, df1=(gene.df.test), df2=gene.df.total, lower.tail=FALSE)

#	Output
	out <- new("DGELRT",list())
	out$comparison <- colnames(design)[coef]
	out$design <- design
	out$coefficients <- as.vector(coefficients)
	
#	Exon level output
	out$exon.df.test <- exon.df.test
	out$exon.df.prior <- s2.fit$df.prior[g]
	out$exon.df.residual <- gene.df.residual[g]
	out$exon.F <- exon.F
	out$exon.p.value <- exon.p.value
	out$genes <- exon.genes
	out$genecolname <- genecolname
	out$exoncolname <- exoncolname
	
#	Gene level output
	out$gene.df.test <- gene.df.test
	out$gene.df.prior <- s2.fit$df.prior
	out$gene.df.residual <- gene.df.residual
	out$gene.Simes.p.value <- gene.Simes.p.value
	out$gene.F <- gene.F
	out$gene.F.p.value <- gene.F.p.value

#	Which columns of exon.genes contain gene level annotation? (from diffSplice in limma)
	exon.lastexon <- cumsum(gene.nexons)
	exon.firstexon <- exon.lastexon-gene.nexons+1
	no <- logical(nrow(exon.genes))
	isdup <- vapply(exon.genes,duplicated,no)[-exon.firstexon,,drop=FALSE]
	isgenelevel <- apply(isdup,2,all)
	out$gene.genes <- exon.genes[exon.lastexon,isgenelevel, drop=FALSE]
	out$gene.genes$NExons <- gene.nexons
	
	out
}


topSpliceDGE <- function(lrt, level="gene", gene.test="Simes", number=10, FDR=1)
# Yunshun Chen and Gordon Smyth
# Created 29 March 2014.  Last modified 24 September 2014. 
{
	level <- match.arg(level,c("exon","gene"))
	gene.test <- match.arg(gene.test,c("Simes","F","f"))
	if(level=="exon") {
		number <- min(number, nrow(lrt$genes))
		P <- lrt$exon.p.value
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number, sum(BH<FDR))
		o <- order(P)[1:number]
		data.frame(lrt$genes[o,,drop=FALSE],logFC=lrt$coefficients[o],F=lrt$exon.F[o],P.Value=P[o],FDR=BH[o])
	} else {
		number <- min(number, nrow(lrt$gene.genes))
		if(gene.test == "Simes") P <- lrt$gene.Simes.p.value else P <- lrt$gene.F.p.value 
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number,sum(BH<FDR))
		o <- order(P)[1:number]
		if(gene.test=="Simes")
			data.frame(lrt$gene.genes[o,,drop=FALSE],P.Value=P[o],FDR=BH[o])
		else
			data.frame(lrt$gene.genes[o,,drop=FALSE],F=lrt$gene.F[o],P.Value=P[o],FDR=BH[o])
	}
}


plotSpliceDGE <- function(lrt, geneid=NULL, rank=1L, FDR = 0.05)
# Plot exons of most differentially spliced gene
# Yunshun Chen and Gordon Smyth
# Created 29 March 2014.  Last modified 24 September 2014.
{
	# Gene labelling including gene symbol
	genecolname <- lrt$genecolname
	genelab <- grep(paste0(genecolname,"|Symbol|symbol"), colnames(lrt$gene.genes), value = T)
	
	if(is.null(geneid)) {
		if(rank==1L)
			i <- which.min(lrt$gene.Simes.p.value)
		else
			i <- order(lrt$gene.Simes.p.value)[rank]
		geneid <- paste(lrt$gene.genes[i,genelab], collapse = ".")
	} else {
		i <- which(lrt$gene.genes[,lrt$genecolname]==geneid)
		geneid <- paste(lrt$gene.genes[i,genelab], collapse = ".")
		if(!length(i)) stop(paste("geneid",geneid,"not found"))
	}
	exon.lastexon <- cumsum(lrt$gene.genes$NExons[1:i])
	j <- (exon.lastexon[i]-lrt$gene.genes$NExons[i]+1):exon.lastexon[i]
	exoncolname <- lrt$exoncolname
	if(is.null(exoncolname)){
		plot(lrt$coefficients[j], xlab = "Exon", ylab = "logFC (this exon vs the average)", main = geneid, type = "b")
	}
	# Plot exons and mark exon ids on the axis
	if(!is.null(exoncolname)) {
		exon.id <- lrt$genes[j, exoncolname]
		xlab <- paste("Exon", exoncolname, sep = " ")
		
		plot(lrt$coefficients[j], xlab = "", ylab = "logFC (this exon vs the average)", main = geneid,type = "b", xaxt = "n")
		axis(1, at = 1:length(j), labels = exon.id, las = 2, cex.axis = 0.6)
		mtext(xlab, side = 1, padj = 5.2)

		# Mark the topSpliced exons
		top <- topSpliceDGE(lrt, number = Inf, level = "exon", FDR = FDR)
		m <- which(top[,genecolname] %in% lrt$gene.genes[i,genecolname])

		if(length(m) > 0){
			if(length(m) == 1) cex <- 1.5 else{
				abs.fdr <- abs(log10(top$FDR[m]))
				from <- range(abs.fdr)
				to <- c(1,2)
				cex <- (abs.fdr - from[1])/diff(from) * diff(to) + to[1]
			}	
			mark <- match(top[m, exoncolname], exon.id)
			points((1:length(j))[mark], lrt$coefficients[j[mark]], col = "red", pch = 16, cex = cex)
		}
	}
	abline(h=0,lty=2)
}











