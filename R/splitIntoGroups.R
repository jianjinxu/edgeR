# Written by Davis McCarthy, February 2009, idea suggested by Mark Robinson
# A function to split the data from a DGEList object according to group
splitIntoGroups<-function(object) {
	nrows<-nrow(object$data)
	y<-lapply(split(t(object$data),object$group), FUN=function(u) matrix(u,nrow=nrows,byrow=TRUE))
	for(i in 1:length(unique(object$group))) {
		rownames(y[[i]])<-rownames(object$data)
	}
	y
}