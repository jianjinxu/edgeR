topTags<-function(object,n=10,adjust.method="BH",sort.by="p.value") 
# Written by Mark Robinson, edited by Davis McCarthy, September 2009
# A function that displays a summary table of the n most differentially expressed tags, judged on p-value
{
	tab<-data.frame(logConc=object$table$logConc,logFC=object$table$logFC,p.value=object$table$p.value,adj.p.val=p.adjust(object$table$p.value,adjust.method))
	rownames(tab)<-rownames(object$table)
	if(sort.by=="logFC")  o<-order(abs(tab$logFC),decreasing=TRUE) else o<-order(tab$p.value)
	new("TopTags",list(table=tab[o[1:min(nrow(tab),n)],],adjust.method=as.character(adjust.method),comparison=as.character(object$comparison)))
}
