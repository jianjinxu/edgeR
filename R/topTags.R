topTags<-function(object,n=10,adjust.method="BH",sort.by="p.value") 
# Written by Mark Robinson, edited by Davis McCarthy, September 2009
# A function that displays a summary table of the n most differentially expressed tags, judged on p-value
{
	
	if(sort.by=="logFC")  o<-order(abs(object$table$logFC),decreasing=TRUE) 
	else o<-order(object$table$p.value)
	adj.p.val <- p.adjust(object$table$p.value,adjust.method)
	chosen <- o[1:min(nrow(object$table),n)]
	tab<-data.frame(logConc=object$table$logConc[chosen], logFC=object$table$logFC[chosen], p.value=object$table$p.value[chosen] ,adj.p.val=adj.p.val[chosen])
	if(!is.null(object$genes)){
		tab <- cbind(object$genes[chosen,], tab) # Assumes that object$genes is a data.frame
	}
	rownames(tab)<-rownames(object$table[chosen,])
	colnames(tab)<-c(colnames(object$genes), "logConc","logFC","PValue","FDR")
	new("TopTags",list(table=tab, adjust.method=as.character(adjust.method),comparison=as.character(object$comparison)))
}
