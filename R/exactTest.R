exactTest<-function(object,pair=NULL, common.disp=TRUE)
# Written by by Davis McCarthy, September 2009
# Calculates exact p-values for the differential expression levels of tags in the two groups being compared
{
	if (!is(object,"DGEList")) stop("Currently only supports DGEList objects as the object argument.")
	levs.group<-levels(object$samples$group)
	if(is.null(pair)) {
		pair<-levels(as.factor(c(levs.group[1],levs.group[2])))
	} else if(!is.factor(pair)) {
		pair<-as.factor(pair)
		#pair<-levels(as.factor(pair))
	}
	if( sum(pair[1]==levs.group)==0 | sum(pair[2]==levs.group)==0 ) 
		stop("At least one element of given pair is not a group.\n Groups are: ",paste(levs.group, collapse=" "),"\n")
	cat("Comparison of groups: ",pair[2],"-",pair[1],"\n")
	
	this.pair<-( pair[1]==object$samples$group | pair[2]==object$samples$group )
	group.pair <- factor(as.vector(object$samples$group[this.pair]))
	levs.pair<-levels(group.pair)
	obj.pair <- DGEList(counts=object$counts[,this.pair], group=group.pair, lib.size=object$samples$lib.size[this.pair])
	if(common.disp)	
		dispersion <- object$common.dispersion 
	else 
		dispersion <- object$tagwise.dispersion
        if( is.null(dispersion) )
		stop("Need a non-null dispersion.  Perhaps you have not run estimateCommonDisp() or estimateTagwiseDisp()?")
	q2q.pair <- equalizeLibSizes(obj.pair,disp=dispersion,null.hypothesis=TRUE)
	mus <- object$common.lib.size*q2q.pair$conc$conc.common
	exact.pvals<-exactTestNB(q2q.pair$pseudo,obj.pair$samples$group,pair,mus,r=1/dispersion)
	
	logConc<-(log2(q2q.pair$conc$conc.group[,pair[1]==levs.pair])+log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]))/2
	logFC<-log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]/q2q.pair$conc$conc.group[,pair[1]==levs.pair])
	de.out<-data.frame(logConc=logConc, logFC=logFC, p.value=exact.pvals)
	rownames(de.out) <- rownames(obj.pair$counts)
	new("deDGEList",list(table=de.out, comparison=pair))
}
