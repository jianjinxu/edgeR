exactTest <- function(object, pair=NULL, dispersion="auto", rejection.region="doubletail", big.count=900)
#	Calculates exact p-values for the differential expression levels of tags in the two groups being compared.
#	Davis McCarthy, Gordon Smyth.
#	Created September 2009. Last modified 10 Jan 2012.
{
	if(!is(object,"DGEList")) stop("Currently only supports DGEList objects as the object argument.")

	object$samples$group <- as.factor(object$samples$group)
	levs.group <- levels(object$samples$group)
	if( is.null(rownames(object$counts)) )
		rownames(object$counts) <- paste("tag",1:nrow(object$counts),sep=".")
	if( is.null(pair) )
		pair <- levs.group[1:2]
	if( !all(pair %in% levs.group) )
		stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "), "\n")
	if(length(pair)!=2) stop("Pair must be of length 2.")
	if(is.numeric(pair)) pair <- levels(object$samples$group)[pair]
	else pair <- as.character(pair)	

	this.pair <- ( object$samples$group %in% pair )
	cat("Comparison of groups: ",as.vector(pair[2]),"-",as.vector(pair[1]),"\n")
	group.pair <- factor(as.vector(object$samples$group[this.pair]))
	levs.pair <- levels(group.pair)
	lib.size <- object$samples$lib.size * object$samples$norm.factors
	obj.pair <- DGEList(counts=object$counts[,this.pair], group=group.pair, lib.size=lib.size[this.pair])

	if(is.null(dispersion)) dispersion <- "auto"
	if(is.character(dispersion)) {
		dispersion <- match.arg(dispersion,c("auto","common","tagwise"))
		dispersion <- switch(dispersion,
			"common"=object$common.dispersion,
			"tagwise"=object$tagwise.dispersion,
			"auto"=if(is.null(object$tagwise.dispersion)) object$common.dispersion else object$tagwise.dispersion
		)
		if(is.null(dispersion)) stop("specified dispersion not found in object")
	}
	ldisp <- length(dispersion)
	ntags <- nrow(object$counts)
	if(ldisp!=1 && ldisp!=ntags) stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
	dispersion <- pmax(dispersion,1e-06)
	if(ldisp==1) dispersion <- rep(dispersion,ntags)

	q2q.pair <- equalizeLibSizes(obj.pair,disp=dispersion,null.hypothesis=TRUE)
	mus <- q2q.pair$N*q2q.pair$conc$conc.common
	y <- splitIntoGroupsPseudo(q2q.pair$pseudo,group.pair,pair)

	rejection.region <- match.arg(rejection.region,c("doubletail","deviance","smallp"))
	exact.pvals <- switch(rejection.region,
		doubletail=exactTestDoubleTail(y1=y$y1,y2=y$y2,dispersion=dispersion,big.count=big.count),
		deviance=exactTestByDeviance(y1=y$y1,y2=y$y2,dispersion=dispersion,big.count=big.count),
		smallp=exactTestBySmallP(y1=y$y1,y2=y$y2,dispersion=dispersion,big.count=big.count)
	)
	
	logConc<-(log2(q2q.pair$conc$conc.group[,pair[1]==levs.pair])+log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]))/2
	logFC<-log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]/q2q.pair$conc$conc.group[,pair[1]==levs.pair])
	logFC[obj.pair$all.zeros] <- 0

	de.out <- data.frame(logFC=logFC, logConc=logConc, PValue=exact.pvals)
	rownames(de.out) <- rownames(obj.pair$counts)
	new("DGEExact",list(table=de.out, comparison=pair, genes=object$genes))
}
