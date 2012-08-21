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

d2 <- estimateTagwiseDisp(d,trend="none")
summary(d2$tagwise.dispersion)
de <- exactTest(d2,dispersion="common")
topTags(de)

de <- exactTest(d2)
topTags(de)

d2 <- estimateTagwiseDisp(d,trend="movingave",span=0.4)
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

summary(exactTest(d2,rejection="smallp")$table$PValue)
summary(exactTest(d2,rejection="deviance")$table$PValue)

d2 <- estimateTagwiseDisp(d,trend="loess",span=0.8)
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

d2 <- estimateTagwiseDisp(d,trend="tricube",span=0.8)
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

# mglmOneWay
design <- model.matrix(~group,data=d$samples)
mglmOneWay(d[1:10,],design,dispersion=0.2)
mglmOneWay(d[1:10,],design,dispersion=0)

fit <- glmFit(d,design,dispersion=d$common.dispersion)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

fit <- glmFit(d,design,dispersion=d$common.dispersion,prior.count.total=2)
summary(fit$coef)

fit <- glmFit(d,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

dglm <- estimateGLMCommonDisp(d,design)
dglm$common.dispersion
dglm <- estimateGLMTagwiseDisp(dglm,design)
summary(dglm$tagwise.dispersion)
fit <- glmFit(dglm,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
dglm <- estimateGLMTrendedDisp(dglm,design)
summary(dglm$trended.dispersion)
dglm <- estimateGLMTrendedDisp(dglm,design,method="power")
summary(dglm$trended.dispersion)
dglm <- estimateGLMTrendedDisp(dglm,design,method="spline")
summary(dglm$trended.dispersion)
dglm <- estimateGLMTrendedDisp(dglm,design,method="bin.spline")
summary(dglm$trended.dispersion)
dglm <- estimateGLMTagwiseDisp(dglm,design)
summary(dglm$tagwise.dispersion)

# Continuous trend
nlibs <- 3
ntags <- 1000
dispersion.true <- 0.1
# Make first transcript respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ntags-1)))
mu.true <- 2^(beta.true %*% t(design))
# Generate count data
y <- rnbinom(ntags*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ntags,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("Gene",1:ntags,sep="")
d <- DGEList(y)
d <- calcNormFactors(d)
fit <- glmFit(d, design, dispersion=dispersion.true)
results <- glmLRT(fit, coef=2)
topTags(results)
d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
glmFit(d,design,dispersion=dispersion.true,method="simple")
glmFit(d,design,dispersion=dispersion.true,method="levenberg")

# Exact tests
y <- matrix(rnbinom(20,mu=10,size=3/2),nrow=5)
group <- factor(c(1,1,2,2))
ys <- splitIntoGroupsPseudo(y,group,pair=c(1,2))
exactTestDoubleTail(ys$y1,ys$y2,dispersion=2/3)

y <- matrix(rnbinom(5*7,mu=10,size=3/2),nrow=5,ncol=7)
group <- factor(c(1,1,2,2,3,3,3))
ys <- splitIntoGroupsPseudo(y,group,pair=c(1,3))
exactTestDoubleTail(ys$y1,ys$y2,dispersion=2/3)
exactTestBetaApprox(ys$y1,ys$y2,dispersion=2/3)

y[1,3:4] <- 0
design <- model.matrix(~group)
fit <- glmFit(y,design,dispersion=2/3)
summary(fit$coef)

# spliceVariants
z = matrix(c(2,0,4,6,4,3,7,1,1,0,1,1,0,3,1,2,0,1,2,1,0,3,1,0), 8, 3)
gz = c(1,2,2,2,2,2,2,2)
spliceVariants(DGEList(counts = z, group = c(1,2,2)), gz)
