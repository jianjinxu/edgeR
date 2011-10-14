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

d <- estimateTagwiseDisp(d,trend="none")
summary(d$tagwise.dispersion)
de <- exactTest(d,dispersion="common")
topTags(de)

de <- exactTest(d)
topTags(de)

d <- estimateTagwiseDisp(d,trend="movingave")
summary(d$tagwise.dispersion)
de <- exactTest(d)
topTags(de)

d2 <- estimateTagwiseDisp(d,trend="tricube")
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

# mglmOneWay
design <- model.matrix(~group,data=d$samples)
mglmOneWay(d[1:10,],design,dispersion=0.2)
mglmOneWay(d[1:10,],design,dispersion=0)

fit <- glmFit(d,design,dispersion=d$common.dispersion)
lrt <- glmLRT(d,fit,coef=2)
topTags(lrt)

fit <- glmFit(d,design)
lrt <- glmLRT(d,fit,coef=2)
topTags(lrt)

example(glmFit)
fit <- glmFit(d,design,dispersion=dispersion.true,method="simple")
fit
fit <- glmFit(d,design,dispersion=dispersion.true,method="levenberg")
fit

y <- matrix(rnbinom(20,mu=10,size=3/2),nrow=5)
group <- factor(c(1,1,2,2))
ys <- splitIntoGroupsPseudo(y,group,pair=c(1,2))
exactTestDoubleTail(ys$y1,ys$y2,dispersion=2/3)

y <- matrix(rnbinom(5*7,mu=10,size=3/2),nrow=5,ncol=7)
group <- factor(c(1,1,2,2,3,3,3))
ys <- splitIntoGroupsPseudo(y,group,pair=c(1,3))
exactTestDoubleTail(ys$y1,ys$y2,dispersion=2/3)

y[1,3:4] <- 0
design <- model.matrix(~group)
fit <- glmFit(y,design,dispersion=2/3)
summary(fit$coef)

