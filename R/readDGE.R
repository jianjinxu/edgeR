
readDGE<-function(files,...) {
  tags<-NULL
  for (fn in files) {
     x<-read.table(fn,...)
     tags<-union(tags,as.character(x[,1]))
  }
  counts<-as.data.frame(matrix(rep(0,length(tags)*length(files)),ncol=length(files)))
  for (i in seq(along=files)) {
     x<-read.table(files[i],...)
     aa<-match(as.character(x[,1]),tags)
     colnames(counts)[i]<-sub(".txt","",files)[i]
     counts[aa,i]<-x[,2]
  }
  rownames(counts)<-tags
  colS<-colSums(counts)
  return(list(data=as.matrix(counts),lib.size=colS))
}

