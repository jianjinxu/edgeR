library(edgeR)

set.seed(0); u <- runif(100)

# generate raw counts from NB, create list object
y <- matrix(rnbinom(80,size=5,mu=10),nrow=20)
y <- rbind(0,c(0,0,2,2),y)
rownames(y) <- paste("Tag",1:nrow(y),sep=".")
d <- DGEList(counts=y,group=rep(1:2,each=2),lib.size=1001:1004)

# estimate common dispersion and find differences in expression
d <- estimateCommonDisp(d)
d$common.dispersion
de <- exactTest(d)
summary(de$table)
topTags(de)

# mglmOneWay
design <- model.matrix(~group,data=d$samples)
mglmOneWay(d[1:10,],design,dispersion=0.2)
mglmOneWay(d[1:10,],design,dispersion=0)

fit <- glmFit(d,design)
lrt <- glmLRT(d,fit,coef=2)
topTags(lrt)

example(glmFit)
fit <- glmFit(d,design,dispersion=dispersion.true,method="simple")
fit
fit <- glmFit(d,design,dispersion=dispersion.true,method="levenberg")
fit

# example using exactTest.matrix directly
y <- matrix(rnbinom(20,mu=10,size=1.5),nrow=5)
group <- factor(c(1,1,2,2))
y <- splitIntoGroupsPseudo(y,group,pair=c(1,2))
mus <- rep(10,5)
f <- exactTest.matrix(y$y1,y$y2,mus,r=1.5,all.zeros=rep(FALSE,length=nrow(y$y1)))
f
