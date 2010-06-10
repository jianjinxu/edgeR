exactTest<-function(object, pair=NULL, dispersion=NULL, common.disp=TRUE)
# Written by by Davis McCarthy, September 2009, last modified 15 December 2009
# Calculates exact p-values for the differential expression levels of tags in the two groups being compared
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
	obj.pair <- DGEList(counts=object$counts[,this.pair], group=group.pair, lib.size=object$samples$lib.size[this.pair])
        tags.retained <- rownames(object$counts) %in% rownames(obj.pair$counts)
        if(!is.null(dispersion)) {
          if( length(dispersion)!=1 && length(dispersion)!=nrow(object$counts) )
            stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.\n")
          if( length(dispersion)==nrow(object$counts) ) {
            if( sum(tags.retained)!=nrow(object$counts) )
              dispersion <- dispersion[tags.retained]
          }
          if(any(dispersion==0))
            dispersion[dispersion==0] <- 1e-06
        } else {
            if(common.disp)	
              dispersion <- object$common.dispersion 
            else {
              tags.retained <- rownames(object$counts) %in% rownames(obj.pair$counts)
              dispersion <- object$tagwise.dispersion[tags.retained]
            }
          }
	q2q.pair <- equalizeLibSizes(obj.pair,disp=dispersion,null.hypothesis=TRUE)
	mus <- object$common.lib.size*q2q.pair$conc$conc.common
	y<-splitIntoGroupsPseudo(q2q.pair$pseudo,group.pair,pair)
	exact.pvals<- exactTest.matrix(y$y1,y$y2,mus,r=1/dispersion)
	logConc<-(log2(q2q.pair$conc$conc.group[,pair[1]==levs.pair])+log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]))/2
	logFC<-log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]/q2q.pair$conc$conc.group[,pair[1]==levs.pair])
	de.out<-data.frame(logConc=logConc, logFC=logFC, p.value=exact.pvals)
	rownames(de.out) <- rownames(obj.pair$counts)
	new("deDGEList",list(table=de.out, comparison=pair, genes=object$genes[tags.retained,]))
}
