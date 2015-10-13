splitIntoGroups<-function(object) 
# Written by Davis McCarthy, February 2009, idea suggested by Mark Robinson
# A function to split the data from a DGEList object according to group
{
	nrows<-nrow(object$counts)
	y<-lapply(split(t(object$counts),object$samples$group), FUN=function(u) matrix(u,nrow=nrows,byrow=TRUE))
	for(i in 1:length(levels(object$samples$group))) {
		rownames(y[[i]])<-rownames(object$counts)
	}
	y
}