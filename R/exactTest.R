exactTest <- function(object, pair=NULL, dispersion="auto", rejection.region="doubletail", big.count=900, prior.count.total=0.5)
#	Calculates exact p-values for the differential expression levels of tags in the two groups being compared.
#	Davis McCarthy, Gordon Smyth.
#	Created September 2009. Last modified 1 March 2012.
{
#	Check input
	if(!is(object,"DGEList")) stop("Currently only supports DGEList objects as the object argument.")
	if(is.null(rownames(object$counts))) rownames(object$counts) <- paste("tag",1:nrow(object$counts),sep=".")
	object$samples$group <- as.factor(object$samples$group)
	levs.group <- levels(object$samples$group)
	if(is.null(pair)) pair <- levs.group[1:2]
	if(length(pair)!=2) stop("Pair must be of length 2.")
	if(is.numeric(pair))
		pair <- levels(object$samples$group)[pair]
	else
		pair <- as.character(pair)	
	if(!all(pair %in% levs.group)) stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "), "\n")
#	cat("Comparison of groups: ",as.vector(pair[2]),"-",as.vector(pair[1]),"\n")

#	Normalized library sizes
	lib.size <- object$samples$lib.size * object$samples$norm.factors
	offset <- log(lib.size)
	lib.size.average <- exp(mean(offset))

#	Get dispersion vector
	if(is.null(dispersion)) dispersion <- "auto"
	if(is.character(dispersion)) {
		dispersion <- match.arg(dispersion,c("auto","common","trended","tagwise"))
		dispersion <- switch(dispersion,
			"common"=object$common.dispersion,
			"trended"=object$trended.dispersion,
			"tagwise"=object$tagwise.dispersion,
			"auto"=getDispersion(object)
		)
		if(is.null(dispersion)) stop("specified dispersion not found in object")
	}
	ldisp <- length(dispersion)
	ntags <- nrow(object$counts)
	if(ldisp!=1 && ldisp!=ntags) stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
	if(ldisp==1) dispersion <- rep(dispersion,ntags)

#	Average abundance
	j <- object$samples$group %in% pair
	abundance <- mglmOneGroup(object$counts[,j,drop=FALSE],dispersion=dispersion,offset=offset[j])
	logCPM <- (abundance+log(1e6))/log(2)

#	logFC
	prior.count <- lib.size[j]
	prior.count <- prior.count.total*prior.count/sum(prior.count)
	j1 <- object$samples$group==pair[1]
	n1 <- sum(j1)
	y1 <- object$counts[,j1,drop=FALSE]
	abundance1 <- mglmOneGroup(y1+matrix(prior.count[j1],ntags,n1,byrow=TRUE),offset=offset[j1])
	j2 <- object$samples$group==pair[2]
	n2 <- sum(j2)
	y2 <- object$counts[,j2,drop=FALSE]
	abundance2 <- mglmOneGroup(y2+matrix(prior.count[j2],ntags,n2,byrow=TRUE),offset=offset[j2])
	logFC <- (abundance2-abundance1)/log(2)

#	Equalize library sizes
	input.mean <- matrix(exp(abundance),ntags,n1)
	output.mean <- input.mean*lib.size.average
	input.mean <- t(t(input.mean)*lib.size[j1])
	y1 <- q2qnbinom(y1,input.mean=input.mean,output.mean=output.mean,dispersion=dispersion)
	input.mean <- matrix(exp(abundance),ntags,n2)
	output.mean <- input.mean*lib.size.average
	input.mean <- t(t(input.mean)*lib.size[j2])
	y2 <- q2qnbinom(y2,input.mean=input.mean,output.mean=output.mean,dispersion=dispersion)

	rejection.region <- match.arg(rejection.region,c("doubletail","deviance","smallp"))
	exact.pvals <- switch(rejection.region,
		doubletail=exactTestDoubleTail(y1,y2,dispersion=dispersion,big.count=big.count),
		deviance=exactTestByDeviance(y1,y2,dispersion=dispersion,big.count=big.count),
		smallp=exactTestBySmallP(y1,y2,dispersion=dispersion,big.count=big.count)
	)

	de.out <- data.frame(logFC=logFC, logCPM=logCPM, PValue=exact.pvals)
	rownames(de.out) <- rownames(object$counts)
	new("DGEExact",list(table=de.out, comparison=pair, genes=object$genes))
}
