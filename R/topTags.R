
topTags<-function(object,n=10,adj.method= "BH") {
  if (!is(object,"deDGEList"))
    stop("Currently only supports deDGEList objects.")
  tab<-data.frame(A=log2(object$ps$p1)+log2(object$ps$p2)/2,M=log2(object$ps$p2/object$ps$p1),P.Value=object$exact,adj.P.Val=p.adjust(object$exact,adj.method))
  rownames(tab)<-rownames(object$pseudo)
  o<-order(tab$P.Value)
  tab[o[1:n],]
}

