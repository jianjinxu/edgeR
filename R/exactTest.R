exactTest <- function(object, pair=NULL, dispersion=NULL, common.disp=TRUE, rejection.region="doubletail", big.count=1000)
#	Calculates exact p-values for the differential expression levels of tags in the two groups being compared.
#	Davis McCarthy, Gordon Smyth.
#	Created September 2009. Last modified 29 Sep 2011.
{
    if( !is(object,"DGEList") )
        stop("Currently only supports DGEList objects as the object argument.")
    if( is.null(dispersion) && is.null(object$common.dispersion) && is.null(object$tagwise.dispersion) )
        stop("Value(s) for the dispersion parameter must be specified. Try running estimateCommonDisp() and/or estimateTagwiseDisp() before exactTest().")
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
    if(!is.null(dispersion)) {
        if( length(dispersion)!=1 && length(dispersion)!=nrow(object$counts) )
            stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.\n")
        if( length(dispersion)==1 )
            dispersion <- rep(dispersion, length=nrow(object$counts))
        if(any(dispersion==0))
            dispersion[dispersion==0] <- 1e-06
    } else {
        if(common.disp)	
            dispersion <- rep(object$common.dispersion, length=nrow(object$counts))
        else {
            dispersion <- object$tagwise.dispersion
        }
    }

    q2q.pair <- equalizeLibSizes(obj.pair,disp=dispersion,null.hypothesis=TRUE)
    mus <- q2q.pair$N*q2q.pair$conc$conc.common
    y<-splitIntoGroupsPseudo(q2q.pair$pseudo,group.pair,pair)

	rejection.region <- match.arg(rejection.region,c("doubletail","deviance","smallp"))
	exact.pvals <- switch(rejection.region,
		doubletail=exactTestDoubleTail(y1=y$y1,y2=y$y2,dispersion=dispersion,big.count=big.count),
		deviance=exactTestDoubleTail(y1=y$y1,y2=y$y2,dispersion=dispersion,big.count=big.count),
		smallp=exactTestDoubleTail(y1=y$y1,y2=y$y2,dispersion=dispersion,big.count=big.count)
	)
	
    logConc<-(log2(q2q.pair$conc$conc.group[,pair[1]==levs.pair])+log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]))/2
    logFC<-log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]/q2q.pair$conc$conc.group[,pair[1]==levs.pair])
    logFC[obj.pair$all.zeros] <- 0

    de.out<-data.frame(logConc=logConc, logFC=logFC, p.value=exact.pvals)
    rownames(de.out) <- rownames(obj.pair$counts)
    new("DGEExact",list(table=de.out, comparison=pair, genes=object$genes))
}
