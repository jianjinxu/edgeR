# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function that displays a summary table of the n most differentially expressed tags, judged on p-value
topTags<-function(object,pair=c(1,2),n=10,adj.method= "BH",verbose=TRUE) {
	if (!is(object,"deDGEList"))
		stop("Currently only supports deDGEList objects.")
	g1<-pair[1]
	g2<-pair[2]
	fisher<-exactTestNB(object$pseudo,object$group,pair,object$M*object$ps$p.common,object$r,verbose)
	tab<-data.frame(A=(log2(object$ps$p.group[,g1])+log2(object$ps$p.group[,g2]))/2,M=log2(object$ps$p.group[,g2]/object$ps$p.group[,g1]),P.Value=fisher,adj.P.Val=p.adjust(fisher,adj.method))
	rownames(tab)<-rownames(object$pseudo)
	o<-order(tab$P.Value)
	tab[o[1:min(nrow(tab),n)],]
}