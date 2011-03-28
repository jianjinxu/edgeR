library(edgeR)

set.seed(0); u <- runif(100)

#y <- matrix(rnbinom(2*4,mu=4,size=1.5),2,4)
#lib.size <- rep(50000,4)
#group <- c(1,1,2,2)

#d<-list(data=y,lib.size=lib.size,group=group)


# generate raw counts from NB, create list object
y<-matrix(rnbinom(80,size=5,mu=10),nrow=20)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$counts)<-paste("tagno",1:nrow(d$counts),sep=".")

# estimate common dispersion and find differences in expression
d<-estimateCommonDisp(d)
de<-exactTest(d)

# example using exactTest.matrix directly
y<-matrix(rnbinom(20,mu=10,size=1.5),nrow=5)
group<-factor(c(1,1,2,2))
y<-splitIntoGroupsPseudo(y,group,pair=c(1,2))
mus<-rep(10,5)
f<-exactTest.matrix(y$y1,y$y2,mus,r=1.5,all.zeros=rep(FALSE,length=nrow(y$y1)))

# mglmOneWay
design <- model.matrix(~group)
mglmOneWay(d[1:10,],design,dispersion=0.2)
