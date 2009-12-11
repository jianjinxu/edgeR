exactTest<-function(object,pair=NULL, common.disp=TRUE)
# Written by by Davis McCarthy, September 2009
# Calculates exact p-values for the differential expression levels of tags in the two groups being compared
{
	if (!is(object,"DGEList")) stop("Currently only supports DGEList objects as the object argument.")
	object$samples$group <- as.factor(object$samples$group)
	levs.group<-levels(object$samples$group)
	if (is.null(pair))
		pair <- levs.group[1:2]
	if( !all(pair %in% levs.group) )
		stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "), "\n")
	stopifnot(length(pair)==2)
	
	this.pair<-( object$samples$group %in% pair )
	cat("Comparison of groups: ",as.vector(pair[2]),"-",as.vector(pair[1]),"\n")
	group.pair <- factor(as.vector(object$samples$group[this.pair]))
	levs.pair<-levels(group.pair)
	obj.pair <- DGEList(counts=object$counts[,this.pair], group=group.pair, lib.size=object$samples$lib.size[this.pair])
	if(common.disp)	
		dispersion <- object$common.dispersion 
	else 
		dispersion <- object$tagwise.dispersion
	q2q.pair <- equalizeLibSizes(obj.pair,disp=dispersion,null.hypothesis=TRUE)
	mus <- object$common.lib.size*q2q.pair$conc$conc.common
	exact.pvals<-exactTest.matrix(q2q.pair$pseudo,obj.pair$samples$group,pair,mus,r=1/dispersion)
	
	logConc<-(log2(q2q.pair$conc$conc.group[,pair[1]==levs.pair])+log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]))/2
	logFC<-log2(q2q.pair$conc$conc.group[,pair[2]==levs.pair]/q2q.pair$conc$conc.group[,pair[1]==levs.pair])
	de.out<-data.frame(logConc=logConc, logFC=logFC, p.value=exact.pvals)
	rownames(de.out) <- rownames(obj.pair$counts)
	new("deDGEList",list(table=de.out, comparison=pair, genes=object$genes))
}