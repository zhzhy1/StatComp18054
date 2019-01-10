## ------------------------------------------------------------------------
m1<-matrix(1,nr=2,nc=2)
m2<-matrix(2,nr=2,nc=2)
# Merge vector or matrix
rbind(m1,m2)
cbind(m1,m2)
# Matrix product
rbind(m1,m2)%*%cbind(m1,m2)
cbind(m1,m2)%*%rbind(m1,m2)
# Diagonal operation and diagonal matrix
diag(m1)
diag(rbind(m1,m2)%*%cbind(m1,m2))
diag(m1)<-10
m1
diag(3)
v<-c(10,20,30)
diag(v)
diag(2.1,nr=3,nc=5)
# Inverse matrix
solve(diag(v))
# Matrix QR decomposition
qr(diag(v))
# Eigenvalue and eigenvector
eigen(diag(v))
# Use table to output eigenvectors
knitr::kable(eigen(diag(v))$vectors,col.names=c("first eigenvector", "second eigenvector", "third eigenvector"), align = c('c','c','c'))
# Singular value decomposition
svd(diag(v))

## ------------------------------------------------------------------------
layout(matrix(1:4,2,2))
mat<-matrix(1:4,2,2)
mat
layout(mat)
layout.show(4)
layout(matrix(1:6,3,2))
layout.show(6)
layout(matrix(1:6,2,3))
layout.show(6)
m<-matrix(c(1:3,3),2,2)
layout(m)
layout.show(3)
m<-matrix(1:4,2,2)
layout(m,widths=c(1,3),heights=c(3,1))
layout.show(4)
m<-matrix(c(1,1,2,1),2,2)
layout(m,widths=c(2,1),heights=c(1,2))
layout.show(2)
m<-matrix(0:3,2,2)
layout(m,c(1,3),c(1,3))
layout.show(3)

## ---- echo=FALSE---------------------------------------------------------
library(knitr)
kable(matrix(c("p(x)",0.1,0.2,0.2,0.2,0.3),nrow=1),col.names = c("x","0","1","2","3","4"),align=c('l','c','c','c','c','c'))

## ------------------------------------------------------------------------
x<-c(0:4)
p<-c(0.1,0.2,0.2,0.2,0.3)
cp<-cumsum(p)
m<-1000 
#the random sample is in numeric r 
r<-numeric(m)
r<-x[findInterval(runif(m),cp)+1]
#the sample frequency is in ct
ct<-as.vector(table(r))
result<-rbind(x,ct/sum(ct)/p)
rownames(result)<-c("x","relative frequency")
result

## ------------------------------------------------------------------------
x<-c(0:4)
p<-c(0.1,0.2,0.2,0.2,0.3)
#the random sample is in numeric r1 
r1<-sample(x,size=1000,replace=TRUE,prob=p)
#the sample frequency is in ct1
ct1<-as.vector(table(r1))
result1<-rbind(x,ct1/sum(ct1)/p)
rownames(result1)<-c("x","relative frequency")
result1

## ------------------------------------------------------------------------
randomBeta<-function(n,a,b){
m=(1-a)/(2-a-b)
max=m^(a-1)*((1-m)^(b-1))/beta(a,b)#maximum of f(x)
c<-max+3 #take c bigger than the maximum
j<-k<-0
y<-numeric(n)
while (k < n) {
u <- runif(1)
j <- j + 1
x <- runif(1) #random variate from g
rho<-x^(a-1)*((1-x)^(b-1))/beta(a,b)/c
if (rho> u) {
#accept x
k <- k + 1
y[k] <- x
}
}
return(list("random sample"=y,"times"=j,"c"=c))
}
result<-randomBeta(n=1000,a=3,b=2)
y<-result[[1]]
hist(y,probability = TRUE,main="Beta(3,2) distribution",ylim = c(0,2))
x<-seq(0,1,0.01)
lines(x,x^2*(1-x)/beta(3,2)) #pdf of Beta(3,2)
times<-result[[2]] #times of acceptance-rejection
times 
c<-result[[3]] #constant c
c

## ------------------------------------------------------------------------
n=1000
r=4
beta=2
lambda=rgamma(n,r,beta)
x=rexp(n,lambda)
hist(x,prob=TRUE,main="Exponential-Gamma mixture")

## ------------------------------------------------------------------------
mypbeta<-function(a,b,x){ # antithetic variables to compute a Monte Carlo estimate of the Beta(3, 3) cdf
  m=1000000
  u=runif(m/2,min=0,max=x)
  u1=1-u
  g=1/beta(a,b)*u^(a-1)*(1-u)^(b-1)*x
  g1=1/beta(a,b)*u1^(a-1)*(1-u1)^(b-1)*x
  theta=(mean(g)+mean(g1))/2
  sd=sd(c(g,g1))/m
  return(list("theta"=theta,"sd"=sd))
}
x=seq(0.1,0.9,by=0.1)
MC=numeric(9)# Monte Carlo estimate
sd=numeric(9)
for(i in 1:9){
  result=mypbeta(3,3,x[i])
  MC[i]=result[[1]]
  sd[i]=result[[2]]
}
pBeta=pbeta(x,3,3)# the values returned by the pbeta function in R
cdf=rbind(MC,pBeta)
knitr::kable(cdf,col.names=x)
matplot(x,cbind(MC,pBeta),col=1:2,pch=1:2,xlim=c(0,1))# Compare using figure
legend("topleft", inset=.05,legend=c("MC","pbeta"),col=1:2,pch=1:2)

## ------------------------------------------------------------------------
myrRayleigh<-function(n,sigma,antithetic=TRUE){
  u=runif(n/2)
  if(!antithetic) 
    {v=runif(n/2)}else
    {v=1-u}
  u=c(u,v)
  x=sqrt(-2*sigma^2*log(1-u))
  return(x)
}
n=1000
set.seed(123)
MC1=myrRayleigh(n,2,antithetic = TRUE)
set.seed(123)
MC2=myrRayleigh(n,2,antithetic = FALSE)
qqplot(MC1,MC2) #Compare the order statistics of the two groups of samples
y=seq(0,10)
lines(y,y,col=2)
1-var(MC1)/var(MC2)#percent reduction in variance

## ------------------------------------------------------------------------
x=seq(1,5,0.1)
g=function(x)x^2/sqrt(2*pi)*exp(-(x^2)/2)
f1=function(x)1/((1-pnorm(-0.5))*sqrt(2*pi))*exp(-((x-1.5)^2)/2)
f2=function(x)x*exp(-((x^2)-1)/2)
plot(x,g(x),col=1,type="o",ylim=c(0,1),lty=1,ylab="y")
points(x,f1(x),col=2,type="o",lty=1)
lines(x,f2(x),col=3,type="o",lty=1)
legend("topright",inset=.05,legend=c("g(x)","f1(x)","f2(x)"),lty=1,col=1:3,horiz=FALSE)

n=1000
set.seed(123)
X=rnorm(1.5*n,mean=1.5,sd=1)
X=X[X>1]
X=X[1:n]
theta1=mean(g(X)/f1(X))
sd1=sd(g(X)/f1(X))/n
theta1
sd1
set.seed(123)
Z=rexp(n,rate=0.5)
Y=sqrt(Z+1)
theta2=mean(g(Y)/f2(Y))
sd2=sd(g(Y)/f2(Y))/n
theta2
sd2


## ------------------------------------------------------------------------
Giniratio=function(n,distribution){
  if(distribution=="lognormal") {
    x=exp(rnorm(n))
  }
  if(distribution=="uniform") {
    x=runif(n)
  }
  if(distribution=="Bernoulli") {
    x=rbinom(n,size=1000,0.1)
  }
  xorder=sort(x)
  G=0
  for(i in 1:n){
    G=(2*i-n-1)*xorder[i]
  }
  G=G/(n^2*mean(x))
  return(G)
}
n=1000
m=1000
Gini=matrix(0,nrow=m,ncol=3)
set.seed(123)
for(i in 1:m){
  Gini[i,1]=Giniratio(n=n,distribution="lognormal")
  Gini[i,2]=Giniratio(n=n,distribution="uniform")
  Gini[i,3]=Giniratio(n=n,distribution="Bernoulli")
}
colMeans(Gini)
decile=matrix(0,9,3)
med=numeric(3)
s=seq.int(0,m-1,by=m/10)[-1]
for(i in 1:3){
  decile[,i]=sort(Gini[,i])[s]
  med[i]=median(Gini[,i])
}
med
decile
hist(Gini[,1],freq=F)
hist(Gini[,2],freq=F)
hist(Gini[,3],freq=F)

## ------------------------------------------------------------------------
Giniratio=function(n,distribution){
  if(distribution=="lognormal") {
    x=exp(rnorm(n))
  }
  if(distribution=="uniform") {
    x=runif(n)
  }
  if(distribution=="Bernoulli") {
    x=rbinom(n,size=1000,0.1)
  }
  xorder=sort(x)
  G=0
  for(i in 1:n){
    G=(2*i-n-1)*xorder[i]
  }
  G=G/(n^2*mean(x))
  return(G)
}
n=1000
m=1000
Gini=numeric(m)
set.seed(123)
for(i in 1:m){
  Gini[i]=Giniratio(n=n,distribution="lognormal")
}
CI=c(mean(Gini)-qt(0.975,1e4-2)*sd(Gini),mean(Gini)+qt(0.975,1e4-2)*sd(Gini))
CI

## ------------------------------------------------------------------------
library(mvtnorm)
test=function(method,r){
n=100
m=1000
mean=c(0,0)
sigma=matrix(c(1,r,r,1),nrow=2)
p=mean(replicate(m,expr={
x=rmvnorm(n,mean,sigma)
cor.test(x[,1],x[,2],method=method)$p.value<=0.05
}))
return(p)
}
test1=function(method){             ##(X,Y)~(X,v*X),  X~N(0,1),P(V=1)=P(V=-1)=0.5
n=100
m=1000
p=mean(replicate(m,expr={
x=rnorm(n)
a=rbinom(n,1,0.5)
v=2*a-1
y=data.frame(x=x,y=v*x)
cor.test(y[,1],y[,2],method="pearson")$p.value<=0.05
}))
return(p)
}
test("pearson",0.2)
test("kendall",0.2)
test("spearman",0.2)
test1("pearson")
test1("kendall")
test1("spearman")

## ------------------------------------------------------------------------
library(bootstrap) #for the law data
law
theta.hat=cor(law$LSAT, law$GPA)
theta.hat
n=nrow(law)
j.cor=function(x,i)cor(x[i,1],x[i,2])
theta.jack=numeric(n)
for(i in 1:n){
  theta.jack[i]=j.cor(law,(1:n)[-i])
}
bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
bias.jack
se.jack=sd(theta.jack)*sqrt((n-1)^2/n)
se.jack

## ------------------------------------------------------------------------
library(boot)
aircondit
b.mean=function(x,i)mean(x[i])
mean.boot=boot(data=aircondit$hours,statistic=b.mean,R=1000)
mean.boot
CI=boot.ci(mean.boot,type=c("norm","basic","perc","bca"))
CI

## ------------------------------------------------------------------------
hist(aircondit$hours,breaks = 50,freq = F)

## ------------------------------------------------------------------------
library(bootstrap)
scor
n=nrow(scor)
sigma.hat=cov(scor)*(n-1)/n
eigenvalues.hat=eigen(sigma.hat)$values
theta.hat=eigenvalues.hat[1]/sum(eigenvalues.hat)
theta.jack=numeric(n)
for(i in 1:n){
  sigma.jack=cov(scor[-i,])*(n-1)/n
  eigenvalues.jack=eigen(sigma.jack)$values
  theta.jack[i]=eigenvalues.jack[1]/sum(eigenvalues.jack)
}
bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
bias.jack
se.jack=sd(theta.jack)*sqrt((n-1)^2/n)
se.jack

## ------------------------------------------------------------------------
library(DAAG) 
attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits

L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)

L4 <- lm(log(magnetic) ~ log(chemical))
plot(log(chemical), log(magnetic), main="Log-Log", pch=16)
logyhat4 <- L4$coef[1] + L4$coef[2] * log(a)
lines(log(a), logyhat4, lwd=2)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n/2) #store the mean of squared residual for leave-two-out samples point
# for n/2-fold(leave-two-out) cross validation
# fit models on leave-two-out samples
for (k in 1:(n/2)) {
index<-c(2*k-1,2*k)  #Subscript of leave-two-out samples point
y <- magnetic[-index]
x <- chemical[-index]

J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[index]
e1[k] <- mean((magnetic[index] - yhat1)^2)

J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[index] +
J2$coef[3] * chemical[index]^2
e2[k] <- mean((magnetic[index] - yhat2)^2)

J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[index]
yhat3 <- exp(logyhat3)
e3[k] <- mean((magnetic[index] - yhat3)^2)

J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[index])
yhat4 <- exp(logyhat4)
e4[k] <- mean((magnetic[index] - yhat4)^2)
}
c(mean(e1), mean(e2), mean(e3), mean(e4))
L2

## ------------------------------------------------------------------------
set.seed(123)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
x
y
CMtest<-function(S1, S2){
  xS1=sort(S1)
  M=length(xS1)
  xS2=sort(S2)
  N=length(xS2)
  
  a=data.frame(val=xS1,rang=seq(M),ens=rep(1,M))
  b=data.frame(val=xS2,rang=seq(N),ens=rep(2,N))
  c=rbind(a,b)
  c=c[order(c$val),]
  c=data.frame(c,rangTot=seq(M+N))
  dtfM=c[which(c$ens==1),]
  dtfN=c[which(c$ens==2),]
  somN = sum( (dtfN$rang - dtfN$rangTot)**2 )
  somM = sum( (dtfM$rang - dtfM$rangTot)**2 )
  U = N*somN + M*somM
  CvM = (  (U / (N*M)) / (N+M) ) - ((4*M*N - 1)/(6*(M+N)))
  return(CvM)
}
R <- 999
z <- c(x, y)
K <- 1:length(z)
n<-length(x)
reps <- numeric(R)
t0 <- CMtest(x, y)
for (i in 1:R) {
  kn <- sample(K, size = n, replace = FALSE)
  x1 <- z[kn]
  y1 <- z[-kn]
  reps[i] <- CMtest(x1, y1)
}
p <- mean(abs(c(t0, reps)) >= abs(t0))
p

## ------------------------------------------------------------------------
# library(RANN)
# library(boot)
# library(energy)
# library(Ball)
# 
# m <- 1e3; k<-3; p<-2; mu <- 0.5; su <- 0.75; set.seed(12)
# n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
# n3 <- 10; n4 <-100  #Unbalanced samples
# 
# Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
# }
# 
# eqdist.nn <- function(z,sizes,k){
#   boot.obj <- boot(data=z,statistic=Tn,R=R,
#   sim = "permutation", sizes = sizes,k=k)
#   ts <- c(boot.obj$t0,boot.obj$t)
#   p.value <- mean(ts>=ts[1])
#   list(statistic=ts[1],p.value=p.value)
# }


## ------------------------------------------------------------------------
#Unequal variances and equal expectations
# p.values1 <- matrix(NA,m,3)
# for(i in 1:m){
#   x <- matrix(rnorm(n1*p),ncol=p);
#   y <- cbind(rnorm(n2),rnorm(n2,sd=su));
#   z <- rbind(x,y)
#   p.values1[i,1] <- eqdist.nn(z,N,k)$p.value
#   p.values1[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#   p.values1[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12)$p.value
# }


## ------------------------------------------------------------------------
# #Unequal variances and unequal expectations
# p.values2 <- matrix(NA,m,3)
# for(i in 1:m){
#   x <- matrix(rnorm(n1*p),ncol=p);
#   y <- cbind(rnorm(n2),rnorm(n2,mean=mu,sd=su));
#   z <- rbind(x,y)
#   p.values2[i,1] <- eqdist.nn(z,N,k)$p.value
#   p.values2[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#   p.values2[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12)$p.value
# }


## ------------------------------------------------------------------------
# #t distribution with 1 df (heavy-tailed distribution)
# p.values3 <- matrix(NA,m,3)
# for(i in 1:m){
#   x <- matrix(rt(n1*p,df=1),ncol=p);
#   y <- cbind(rt(n2,df=1),rt(n2,df=1));
#   z <- rbind(x,y)
#   p.values3[i,1] <- eqdist.nn(z,N,k)$p.value
#   p.values3[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#   p.values3[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12)$p.value
# }

## ------------------------------------------------------------------------
# #Unbalanced samples (say, 1 case versus 10 controls)
# p.values4 <- matrix(NA,m,3)
# for(i in 1:m){
#   x <- matrix(rnorm(n3*p),ncol=p);
#   y <- cbind(rnorm(n4),rnorm(n4,mean = mu));
#   z <- rbind(x,y)
#   p.values4[i,1] <- eqdist.nn(z,N,k)$p.value
#   p.values4[i,2] <- eqdist.etest(z,sizes=c(n3,n4),R=R)$p.value
#   p.values4[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12)$p.value
# }
# alpha <- 0.1;
# pow1 <- colMeans(p.values1<alpha)
# pow2 <- colMeans(p.values2<alpha)
# pow3 <- colMeans(p.values3<alpha)
# pow4 <- colMeans(p.values4<alpha)
# pow1
# pow2
# pow3
# pow4

## ------------------------------------------------------------------------
set.seed(1)
rw.Metropolis <- function(n, sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (dt(y, n) / dt(x[i-1], n)))
                x[i] <- y  else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
        }

    n <- 1  #degrees of freedom for target Student t dist.
    N <- 5000
    sigma <- c(.05, .5, 2,  16)

    x0 <- 25
    rw1 <- rw.Metropolis(n, sigma[1], x0, N)
    rw2 <- rw.Metropolis(n, sigma[2], x0, N)
    rw3 <- rw.Metropolis(n, sigma[3], x0, N)
    rw4 <- rw.Metropolis(n, sigma[4], x0, N)

    #number of candidate points rejected
    print(c(rw1$k, rw2$k, rw3$k, rw4$k))
    
    # par(mfrow=c(2,2))  
     refline <- qt(c(.025, .975), df=n)
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    # for (j in 1:4) {
    #     plot(rw[,j], type="l",
    #          xlab=bquote(sigma == .(round(sigma[j],3))),
    #          ylab="X", ylim=range(rw[,j]))
    #     abline(h=refline)
    # }
    # par(mfrow=c(1,1)) #reset to default
    
    Nd=1000  #number of discarded sample
    a <- c(.05, seq(.1, .9, .1), .95)
    Q <- qt(a, n)
    mc <- rw[(Nd+1):N, ]
    Qrw <- apply(mc, 2, function(x) quantile(x, a))
    knitr::kable(cbind(Q, Qrw),col.names=c("True","sigma=0.05","sigma=0.5","sigma=2","sigma=16")) 


## ------------------------------------------------------------------------
set.seed(12)
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

size<-c(125,18,20,34)
prob <- function(y, size) {
        # computes (without the constant) the target density
        if (y < 0 || y >= 1)
            return (0)
        return((1/2+y/4)^size[1] *((1-y)/4)^size[2] * ((1-y)/4)^size[3] *(y/4)^size[4])
}

thetachain<-function(m){ #generate a Metropolis chain for theta
w<-0.5
u <- runif(m)  #for accept/reject step
v<-runif(m,-w,w) #proposal distribution
x<-rep(0,m)
x[1] <- 0.5
for (i in 2:m) {
    y <- x[i-1]+v[i] 
    if (u[i] <= prob(y, size) / prob(x[i-1], size))
        x[i] <- y  else
   x[i] <- x[i-1]
}
return(x)
}
#generate the chains
k=4      #number of chains to generate
n=15000  #length of chains
b=1000   #burn-in length
theta <- matrix(0, nrow=k, ncol=n)
for (i in 1:k) 
  theta[i, ] <- thetachain(m=n)
#compute diagnostic statistics
psi <- t(apply(theta, 1, cumsum))
for (i in 1:nrow(psi)) 
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
#par(mfrow=c(2,2))
#for (i in 1:k)
#plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))
#par(mfrow=c(1,1)) #restore default
#the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
#plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#abline(h=1.2, lty=2)


## ------------------------------------------------------------------------
#par(mfrow=c(2,2))
for (i in 1:k)
#hist(psi[i, 2001:n],xlab=i,freq=FALSE)
#par(mfrow=c(1,1)) #restore default
theta=rowMeans(psi[,2001:n]) 
#estimation of theta 
mean(theta)

## ------------------------------------------------------------------------
set.seed(12)
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

size<-c(125,18,20,34)
prob <- function(y, size) {
        # computes (without the constant) the target density
        if (y < 0 || y >= 1)
            return (0)
        return((1/2+y/4)^size[1] *((1-y)/4)^size[2] * ((1-y)/4)^size[3] *(y/4)^size[4])
}

thetachain<-function(m,x0){ #generate a Metropolis chain for theta
u <- runif(m)  #for accept/reject step
x<-rep(0,m)
x[1] <- x0
for (i in 2:m) {
    y <- runif(1) 
    if (u[i] <= prob(y, size) / prob(x[i-1], size))
        x[i] <- y  else
   x[i] <- x[i-1]
}
return(x)
}
#generate the chains
k=4      #number of chains to generate
n=15000  #length of chains
b=1000   #burn-in length
x0=c(0.3,0.4,0.5,0.6)
theta <- matrix(0, nrow=k, ncol=n)
for (i in 1:k) 
  theta[i, ] <- thetachain(m=n,x0=x0[i])
#compute diagnostic statistics
psi <- t(apply(theta, 1, cumsum))
for (i in 1:nrow(psi)) 
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
#par(mfrow=c(2,2))
#for (i in 1:k)
#plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))
#par(mfrow=c(1,1)) #restore default
#the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
#plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#abline(h=1.2, lty=2)


## ------------------------------------------------------------------------
nc=8000
thetahat=mean(rowMeans(psi[,nc:n])) 
thetahat

## ------------------------------------------------------------------------

mys<-function(a,k){
    M=sqrt(a^2*k/(k+1-a^2))
    return(pt(M,df=k,lower.tail = FALSE))
}
myf<-function(a,k){
 mys(a=a,k=(k-1))-mys(a=a,k=k)
}
kc=c(4:25,100,500,1000)
A=rep(0,length(kc))
for(i in 1:length(kc)){
  A[i]=uniroot(myf,c(0.5,sqrt(kc[i]-1)),k=kc[i])$root
}
cbind(kc,A)
plot(kc,A,type="o")


## ------------------------------------------------------------------------
cauchy<-function(x,eta,theta){
   1/(theta*pi*(1+((x-eta)/theta)^2))  
}
mypcauchy<-function(x,eta,theta){
  integrate(cauchy,lower=-Inf,upper=x,eta=eta,theta=theta)
}
x=seq(-2,2,by=0.5)
eta=1
theta=2
mycdf=rep(0,length(x))
for(i in 1:length(x)){
  mycdf[i]=mypcauchy(x[i],eta = eta,theta=theta)$value
}
Rcdf=pcauchy(x,location = eta,scale=theta)
cbind(x,mycdf,Rcdf)

## ------------------------------------------------------------------------
nA=28
nB=24
nOO=41
nAB=70
l<-function(theta,nAA,nBB){
  p=theta[1]
  q=theta[2]
  r=1-p-q
  return(-(2*nAA*log(p)+2*nBB*log(q)+(nA-nAA)*log(2*p*r)+(nB-nBB)*log(2*q*r)+nAB*log(2*p*q)+2*nOO*log(r)))
  }
MLE<-function(p0,q0){
  r0=1-p0-q0
  nAA0=nA*p0^2/(p0^2+2*p0*r0)
  nBB0=nB*q0^2/(q0^2+2*q0*r0)
  return(optim(c(p0,q0),l,lower=c(0.01,0.01),upper=c(0.49,0.49),method = "L-BFGS-B",nAA=nAA0,nBB=nBB0)$par)
}
EM<-function(p0,q0,l,maxiter){
  theta0=c(p0,q0)
  theta=MLE(p0,q0)
  record=rbind(theta0,theta)
  iter=1
  while((max(abs(theta-c(p0,q0)))>l)&(iter<=maxiter)){
    p0=theta[1]
    q0=theta[2]
    theta=MLE(p0,q0)
    record=rbind(record,theta)
    iter=iter+1
  }
  return(list(record=record,iter=iter))
}
result=EM(0.4,0.2,0.000001,1000)
result #the first column records MLE of p, the second column records MLE of q
matplot(result$record,type="o",pch=1:2,col=1:2,xlab="iteration time",ylab="MLE",main="EM-record")
legend("topright",inset=.05,legend=c("p","q"),
       pch=1:2,col=1:2,horiz=FALSE,cex=0.7)

## ------------------------------------------------------------------------
set.seed(12)
attach(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
result1=lapply(formulas,lm)
result1
result2=result1
for(i in 1:4){
  result2[[i]]<-lm(formulas[[i]])
}
result2

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]})
result3=lapply(bootstraps,function(a){lm(a$mpg~a$disp)})
result3
result4=result3
for(i in 1:10){
  result4[[i]]=lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp)
}
result4

## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
lapply(result1,rsq)
lapply(result2,rsq)
lapply(result3,rsq)
lapply(result4,rsq)

## ------------------------------------------------------------------------
set.seed(12)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(1:100, function(i){trials[[i]]$p.value})

## ------------------------------------------------------------------------
library(foreach)
library(datasets)

formulas <- list(
   mtcars$mpg ~ mtcars$disp,
   mtcars$mpg ~ I(1 / mtcars$disp),
   mtcars$mpg ~ mtcars$disp + mtcars$wt,
   mtcars$mpg ~ I(1 / mtcars$disp) + mtcars$wt
)

foreach(i=1:4) %do%
  lm(formulas[[i]])

## ------------------------------------------------------------------------
set.seed(12)
x=sample(c(1,2,3),1000,prob =c(0.2,0.5,0.3) ,replace = TRUE)
y=sample(c(4,5,6),1000,prob =c(0.4,0.4,0.2) ,replace = TRUE)
my_chisqtest<-function(x,y){
  A=table(x,y)
  n=nrow(A)
  p=ncol(A)
  N=sum(A)
  rowp=rowSums(A)/N
  colp=colSums(A)/N
  E=N*rowp%*%t(colp)
  chisq=sum((A-E)^2/E)
  return(chisq)
}
library(microbenchmark)
microbenchmark(chisq.test(x,y),times=1000)
microbenchmark(my_chisqtest(x,y),times=1000)
chisq.test(x,y)
my_chisqtest(x,y)

## ------------------------------------------------------------------------
my_table1<-function(x,y){
  x=factor(x)
  y=factor(y)
  xl=levels(x)
  yl=levels(y)
  n=length(x)
  nxl=length(xl)
  nyl=length(yl)
  table=matrix(0,nxl,nyl,dimnames = list(xl,yl))
  for(i in 1:n){
    for(j in 1:nxl){
      for(k in 1:nyl){
        if(x[i]==xl[j]&y[i]==yl[k])
          table[j,k]=table[j,k]+1
      }
    }
  }
  return(table)
}
library(microbenchmark)
microbenchmark(my_table1(x,y),times=10)
microbenchmark(table(x,y),times=10)
my_table1(x,y)
table(x,y)

## ------------------------------------------------------------------------
my_table2<-function(x,y){
  x=factor(x)
  y=factor(y)
  xl=levels(x)
  yl=levels(y)
  n=length(x)
  nxl=length(xl)
  nyl=length(yl)
  table=matrix(0,nxl,nyl,dimnames = list(xl,yl))
  for(i in 1:nxl){
    xi=factor(x,levels=xl[i])
    I=which(!is.na(xi))
    table[i,]=tabulate(y[I])
  }
  return(table)
}
library(microbenchmark)
microbenchmark(my_table2(x,y),times=10)
microbenchmark(table(x,y),times=10)
my_table2(x,y)
table(x,y)

## ------------------------------------------------------------------------
my_chisqtest1<-function(x,y){
  A=my_table2(x,y)
  n=nrow(A)
  p=ncol(A)
  N=sum(A)
  rowp=rowSums(A)/N
  colp=colSums(A)/N
  E=N*rowp%*%t(colp)
  chisq=sum((A-E)^2/E)
  return(chisq)
}
library(microbenchmark)
microbenchmark(chisq.test(x,y),times=1000)
microbenchmark(my_chisqtest1(x,y),times=1000)
chisq.test(x,y)
my_chisqtest1(x,y)

