
set.seed(0); u <- runif(100)

y <- matrix(rnbinom(2*4,mu=4,size=1.5),2,4)
lib.size <- rep(50000,4)
group <- c(1,1,2,2)

d<-list(data=y,lib.size=lib.size,group=group)

### msage
#msage<-msage(d,alpha=100)  # common likelihood

