
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function that displays a summary table of the n most differentially expressed tags, judged on p-value
topTags<-function(object,pair=NULL,n=10,adj.method= "BH",verbose=TRUE) {
	if (!is(object,"deDGEList"))
		stop("Currently only supports deDGEList objects.")
	levs.group<-levels(object$group)
	if(is.null(pair)) {
		pair<-levels(as.factor(c(levs.group[1],levs.group[2])))
		cat("No pair supplied, so comparing groups ",levs.group[1]," and ",levs.group[2],"\n")
	} else if(!is.factor(pair)) {
		pair<-levels(as.factor(pair))
	}
	if( sum(pair[1]==levs.group)==0 | sum(pair[2]==levs.group)==0 ) 
		stop("At least one element of given pair is not a group.\n Groups are: ",levs.group)
	g1<-pair[1]
	g2<-pair[2]
	fisher<-exactTestNB(object$pseudo,object$group,pair,object$M*object$ps$p.common,object$r,verbose)
	tab<-data.frame(A=(log2(object$ps$p.group[,g1==levs.group])+log2(object$ps$p.group[,g2==levs.group]))/2,M=log2(object$ps$p.group[,g2==levs.group]/object$ps$p.group[,g1==levs.group]),P.Value=fisher,adj.P.Val=p.adjust(fisher,adj.method))
	rownames(tab)<-rownames(object$pseudo)
	o<-order(tab$P.Value)
	tab[o[1:min(nrow(tab),n)],]
}