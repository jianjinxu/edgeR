goana.DGELRT <- function(de, geneid = rownames(de), FDR = 0.05, species = "Hs", trend = FALSE, ...)
#  Gene ontology analysis of DE genes from linear model fit
#  Gordon Smyth and Yifang Hu
#  Created 25 August 2014.
{
# Check fit
	if(is.null(de$table$PValue)) stop("p value not found in fit object (from eBayes).")	
	ngenes <- nrow(de)

	# Check geneid
	# Can be either a vector of IDs or a column name
	geneid <- as.character(geneid)
	if(length(geneid) == ngenes) {
		universe <- geneid
	} else
		if(length(geneid) == 1L) {
			universe <- de$genes[[geneid]]
			if(is.null(universe)) stop(paste("Column",geneid,"not found in de$genes"))
		} else
			stop("geneid has incorrect length")

	# Check trend
	# Can be logical, or a numeric vector of covariate values, or the name of the column containing the covariate values
	if(is.logical(trend)) {
		if(trend) {
			covariate <- de$table$logCPM
			if(is.null(covariate)) stop("Amean not found in fit")
		}
	} else
		if(is.numeric(trend)) {
			if(length(trend) != ngenes) stop("If trend is numeric, then length must equal nrow(de)")
			covariate <- trend
			trend <- TRUE
		} else {
			if(is.character(trend)) {
				if(length(trend) != 1L) stop("If trend is character, then length must be 1")
				covariate <- de$genes[[trend]]
				if(is.null(covariate)) stop(paste("Column",trend,"not found in de$genes"))
				trend <- TRUE
			} else
				stop("trend is neither logical, numeric nor character")
		}

	# Check FDR
	if(!is.numeric(FDR) | length(FDR) != 1) stop("FDR must be numeric and of length 1.")
	if(FDR < 0 | FDR > 1) stop("FDR should be between 0 and 1.")

	# Get up and down DE genes
	fdr.coef <- p.adjust(de$table$PValue, method = "BH")
	EG.DE.UP <- universe[fdr.coef < FDR & de$table$logFC > 0]
	EG.DE.DN <- universe[fdr.coef < FDR & de$table$logFC < 0]
	de.gene <- list(Up=EG.DE.UP, Down=EG.DE.DN)

	# Fit monotonic cubic spline for DE genes vs. gene.weights
	if(trend) {
			PW <- isDE <- rep(0,ngenes)
			isDE[fdr.coef < FDR] <- 1
			o <- order(covariate)
			PW[o] <- tricubeMovingAverage(isDE[o],span=0.5,full.length=TRUE)
	}
	if(!trend) PW <- NULL

	NextMethod(de = de.gene, universe = universe, species = species, weights = PW, ...)
}
