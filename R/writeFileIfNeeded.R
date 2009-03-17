`writeFileIfNeeded` <-
function(tabOrFile,label) {
  if( is.matrix(tabOrFile) | is.data.frame(tabOrFile) ) {
    fn<-paste(label,"data",sep=".")
    cat("Writing file:",fn,"\n")
    write.table(tabOrFile,fn,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
    return(list(file=fn,written=TRUE))
  } else if (file.exists(tabOrFile)) {
    return(list(file=tabOrFile,written=FALSE))
  } else {
    stop("Either 'tabOrFile' is not a matrix/data.frame or the file specified does not exist")
  }
}

