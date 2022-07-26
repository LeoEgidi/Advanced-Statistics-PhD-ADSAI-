## ----prelim,echo=FALSE,include=FALSE,warning=FALSE--------------------------------------------------------------------------
library(knitr)
options(width=50)
mypref=function() opts_chunk$set(comment=NA, fig.width=6, fig.height=4,out.width='0.7\\textwidth',echo=FALSE,results='hide',fig.path="figure/050-mixedmodel-",tidy=FALSE)
mypref()
library(xtable)
library(splines)
library(mgcv)
library(coda)
library(rstan)


## ---------------------------------------------------------------------------------------------------------------------------
cmb=read.table("wmap.dat",header=TRUE)
lidar=read.table("lidar.dat",header=TRUE)
lidar=lidar[sort.list(lidar$range),]
bpd=read.table("bpd.dat",header=TRUE)
#cmb[1,]
#plot(cmb[,1],cmb[,2])


## ---------------------------------------------------------------------------------------------------------------------------
lidar$range=(lidar$range-mean(lidar$range))/sd(lidar$range)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
library(rstan)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
lidar=read.table("lidar.dat",header=TRUE)
x=(lidar$range-mean(lidar$range))/sd(lidar$range)
y=lidar$logratio


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
X=cbind(1,x,x^2,x^3)
dim(X)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
knots=seq(min(x),max(x),length=12)[2:11]


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
Z=matrix(NA,221,20)
for (i in 1:221) for (j in 1:20) 
  Z[i,j]=ifelse(x[i]-knots[j]>0,(x[i]-knots[j])^3,0)


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
Z=outer(x,knots,
        FUN=function(x,y) ifelse(x-y>0,(x-y)^3,0))


## ----echo=TRUE,tidy=FALSE,out.height="0.5\\textheight",out.width="0.7\\textwidth",fig.width=7,fig.height=5------------------
matplot(x,cbind(X,Z),type="l",col="black",lty=1)


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
dati1=list(N=length(y),
           p=3,
           K=length(knots),
           y=y,
           X=X,
           Z=Z)


## ----echo=TRUE,tidy=FALSE,size='tiny',eval=FALSE----------------------------------------------------------------------------
## data {
##   int<lower=0> N;
##   int<lower=0> p;
##   int<lower=0> K;
##   vector[N] y;
##   matrix[N, 4] X;
##   matrix[N, K] Z;
## }
## parameters {
##   vector[4] beta;
##   vector[K] b;
##   real<lower=0> sigma2;
##   real<lower=0> sigmab2;
## }
## model {
##   for (i in 1:N) {
##    y[i] ~ normal( X[i,]*beta[1:(p+1)] + Z[i,]*b[1:K],sqrt(sigma2));
##   }
##   for (k in 1:K){
##     b[k] ~ normal(0,sqrt(sigmab2));
##   }
##   for (j in 1:(p+1)) {
##     beta[j] ~ normal( 0  ,1000);
##   }
##   sigma2 ~ inv_gamma(0.001,0.001);
##   sigmab2 ~ inv_gamma(0.001, 0.001);
## }


## ----echo=TRUE,cache=TRUE,cache.path='cache/mcmc1',tidy=FALSE,size='tiny',eval=FALSE----------------------------------------
## fit_spline1 <- stan(file = 'spline1.stan', data = dati1, cores=parallel::detectCores())


## ----echo=FALSE,results='hide'----------------------------------------------------------------------------------------------
#save("fit_spline1",file="fit_spline1.Rdata")
load(file="fit_spline1.Rdata")


## ----echo=TRUE,size='tiny',results='markup'---------------------------------------------------------------------------------
print(fit_spline1)


## ----echo=TRUE,size='tiny',results='markup'---------------------------------------------------------------------------------
posterior=as.array(fit_spline1)
str(posterior)


## ----echo=TRUE,out.width="0.8\\textwidth", out.height="0.4\\textwidth",eval=TRUE,fig.width=8,fig.height=4-------------------
par(mfrow=c(1,2))
plot(posterior[,4,5],type="l")
plot(density(posterior[,4,5]))


## ----tpsfunction,out.width="0.8\\textwidth", out.height="0.8\\textheight",echo=FALSE----------------------------------------
tps=function(x,nodi=NULL){
  if (is.null(nodi)) nodi=quantile(x,seq(0.1,0.9,by=0.2))
  if (length(nodi)==1) nodi=seq(min(x),max(x),length=nodi+2)[2:(nodi+1)]
  K=length(nodi)
  X=matrix(NA,ncol=2+K,nrow=length(x))
  X[,1]=1
  X[,2]=x
  for (k in 1:K) X[,2+k]=abs(x-nodi[k])^3
  D=matrix(0,nrow=2+K,ncol=2+K)
  for (i in 1:K) for (j in 1:K) D[2+i,2+j]=abs(nodi[i]-nodi[j])
  return(list(X=X,D=D))
}
y.teo=function(y,x,l=0,nodi=NULL,p=3){
  bb=tpb(x,nodi,p)
  X=bb$X
  D=bb$D
  part=solve(t(X) %*% X + l*D ) %*% t(X)
  dof=sum(diag(part %*% X))
  yteo=X %*% part %*% y
  return(list(yteo=yteo,dof=dof))
}

par(mfrow=c(2,1),mar=c(2,3,0.2,0.2))
nodi=c(-1,0,0.5,1)
a=tps(x,nodi=nodi)
#plot(x,a$X[,]%*%c(1,1,1,1,1,2),type="l")
mm=matrix(NA,nrow=length(x),ncol=20)
for (i in 1:20){
  mm[,i]=a$X[,]%*%matrix(rnorm(ncol(a$X),0,10000),ncol=1)
}
matplot(x,mm,type="l")
matplot(x,a$X[,],type="l",xlab="x",ylab="")
points(nodi,rep(0,length(nodi)),pch=20)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
X=cbind(1,x)


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
Z_K<-(abs(outer(x,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% 
                    (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))


## ----out.width="0.8\\textwidth", out.height="0.8\\textheight",echo=FALSE----------------------------------------------------
par(mfrow=c(2,1),mar=c(3,4,0.3,0.3))
matplot(x,Z_K,type="l",col="black",lty=1)
matplot(x,Z,type="l",col="black",lty=1)


## ----echo=TRUE,tidy=FALSE,size='footnotesize',eval=FALSE--------------------------------------------------------------------
## model {
##   for (i in 1:N) {
##    y[i] ~ normal( X[i,]*beta[1:2] + Z[i,]*b[1:K],sigma);
##   }
##   for (k in 1:K){
##     b[k] ~ normal(0,sigmab);
##   }
##   for (j in 1:2) {
##     beta[j] ~ normal( 0  ,1000);
##   }
## }


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
dati2 <- list(N=length(y),p=3,
                         K=length(knots),
                         y=y,X=X,Z=Z)


## ----echo=FALSE,results='hide'----------------------------------------------------------------------------------------------
#save("fit_spline2",file="fit_spline2.Rdata")
load(file="fit_spline2.Rdata")
posterior=as.array(fit_spline2)


## ----echo=TRUE,cache=TRUE,cache.path='cache/mcmc2',tidy=FALSE,results='hide',eval=FALSE-------------------------------------
## fit_spline2 <- stan(file = 'spline2.stan', data = dati2, cores=parallel::detectCores())
## posterior=as.array(fit_spline2)


## ----echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------
## launch_shinystan(fit_spline2)


## ----echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------
## str(posterior)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
dimnames(posterior)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
posterior=extract(fit_spline2)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
str(posterior)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
beta=posterior$beta
dim(beta)


## ----grafposteriorsbeta,echo=FALSE,fig.width=9,fig.height=4,out.height="0.4\\textheight",out.width="0.9\\textwidth"---------
par(mfrow=c(1,2),mar=c(2,4,2,0.2))
hist(beta[,1],br=30,freq=FALSE,main=expression(paste(pi,group("(",beta[1],")"))))
lines(density(beta[,1]),col="red",lwd=2)
hist(beta[,2],br=30,freq=FALSE,main=expression(paste(pi,group("(",beta[2],")"))))
lines(density(beta[,2]),col="red",lwd=2)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
b=posterior$b
dim(b)


## ----grafposteriorsb,echo=FALSE,fig.width=9,fig.height=4,out.height="0.4\\textheight",out.width="0.9\\textwidth"------------
par(mfrow=c(1,1),mar=c(2,2,0.2,0.2))
boxplot(as.matrix(b))


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
f=X%*% t(beta) + Z %*% t(b)
dim(f)


## ----grafposteriorsplines,echo=-1,out.height="0.7\\textheight",out.width="0.8\\textwidth",fig.width=8,fig.height=7----------
par(mfrow=c(1,1),mar=c(5,4,0.3,0.3))
plot(x,y,pch=20)
lines(x,f[,1],col=gray(0.8)) 
for (i in 1:4000){ 
  lines(x,f[,i],col=gray(0.8))
}


## ----grafposteriors,echo=FALSE,out.height="0.6\\textheight",out.width="0.8\\textwidth",fig.width=8,fig.height=6-------------
par(mfrow=c(1,1),mar=c(5,4,0.3,0.3),xpd=FALSE)
plot(x,y,pch=20,col=gray(0.8))
for (pp in c(25,50,75,100,120,150:160,200)){ 
  points(x[pp],y[pp],col="red",pch=20)
  dd=density(f[pp,])
  segments(x[pp],min(y[pp],dd$x)-0.1,x[pp],max(y[pp],dd$x)+0.1,col=gray(0.5))
  lines(x[pp]-dd$y*0.01,dd$x)
}


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
pe=apply(f,1,median)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
ci=apply(f,1,quantile,probs=c(0.025,0.5,0.975))


## ----grafposteriorsconstime,echo=FALSE,out.height='0.6\\textheight',out.width='0.8\\textwidth',fig.width=8,fig.height=6-----
par(mfrow=c(1,1),mar=c(5,4,0.3,0.3),xpd=FALSE)
plot(x,y,pch=20,col=gray(0.8))
for (pp in c(25,50,75,100,120,150:160,200)){ 
  points(x[pp],y[pp],col="red",pch=20)
  dd=density(f[pp,])
  segments(x[pp],min(y[pp],dd$x)-0.1,x[pp],max(y[pp],dd$x)+0.1,col=gray(0.5))
  lines(x[pp]-dd$y*0.01,dd$x)
}
lines(x,pe,col="blue")
matlines(x,t(ci),lty=2,col="blue")


## ----echo=TRUE,eval=TRUE,out.height="0.6\\textheight",out.width="0.6\\textwidth"--------------------------------------------
plot(x,y,pch=20)
lines(x,pe,col="red")
matlines(x,t(ci),lwd=2,col="red",lty=c(2,1,2))


## ----echo=FALSE,out.width="0.8\\textwidth",fig.width=4,fig.height=2---------------------------------------------------------
intnodi=seq(0,1,length=7)
aggnodi=c(seq(-0.5,-0.1,length=4),seq(1.1,1.5,length=4))
par(mar=c(4,0,0,0),mfrow=c(1,1))
plot(c(aggnodi,intnodi),rep(0,length=15),xaxt="n",yaxt="n",xlab="",ylab="",bty="n",pch="|",xlim=c(-0.7,1.7))
abline(h=0)
text(intnodi,rep(-0.4,length(intnodi)),labels=c(expression(tau[1]),expression(tau[2]),"...","...","...",expression(tau[K-1]),expression(tau[K])),cex=0.7)
text(aggnodi,rep(-0.4,length(intnodi)),labels=c(expression(xi[1]),"...",expression(xi[M]),"a","b",expression(nu[1]),"...",expression(nu[M])),cex=0.7)
text(aggnodi,rep(-0.8,length(aggnodi)),labels=c(expression(kappa[1]),"...",expression(kappa[M]),"","",expression(kappa[K+M+1]),"...",expression(kappa[K+2*M])),cex=0.7)
text(intnodi,rep(-0.8,length(intnodi)),labels=c(expression(kappa[M+1]),expression(kappa[M+2]),"...","...","...",expression(kappa[M+K-1]),expression(kappa[M+K])),cex=0.7)


## ----bsplineA,out.width="0.8\\textwidth", out.height="0.8\\textheight",echo=FALSE-------------------------------------------
nodi=c(-1,0,0.5,1)
par(mfrow=c(3,1),mar=c(3,3,2,0.2))
lidar$range=x
base=bs(x,knots=nodi,degree=1,intercept=TRUE)
matplot(x,base,type="l",xlab="range (standardized)",main="1st degree B-splines")
points(nodi,rep(0,length(nodi)),pch=20)
base=bs(x,knots=nodi,degree=2,intercept=TRUE)
matplot(x,base,type="l",xlab="range (standardized)",main="2nd degree B-splines")
points(nodi,rep(0,length(nodi)),pch=20)
base=bs(x,knots=nodi,degree=3,intercept=TRUE)
matplot(x,base,type="l",xlab="range (standardized)",main="3rd degree B-splines")
points(nodi,rep(0,length(nodi)),pch=20)


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
nodi <- seq(min(x),max(x),length=22)[2:21]
B <- bs(x,knots=nodi,degree=3,
        Boundary.knots=c(min(x),max(x)),
        intercept=TRUE)


## ----echo=FALSE,out.width="0.8\\textwidth",out.height="0.8\\textheight"-----------------------------------------------------
par(mar=c(5,4,0.2,0.2))
matplot(x,B,type="l")


## ----createOmega,echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------
formOmega <- function(a,b,intKnots) {
  allKnots <- c(rep(a,4),intKnots,rep(b,4))
  K <- length(intKnots) ; L <- 3*(K+8)
  xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+
    rep(allKnots,each=3)[-c(1,2,L)])/2
  wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
  Bdd <- spline.des(allKnots,xtilde,
                    derivs=rep(2,length(xtilde)),
                    outer.ok=TRUE)$design
  Omega <- t(Bdd*wts) %*% Bdd
  return(Omega)
}


## ----compOmega,echo=TRUE----------------------------------------------------------------------------------------------------
Omega <- formOmega(min(x),max(x),nodi)
eigOmega <- eigen(Omega)
indsZ <- 1:(20+2)
UZ <- eigOmega$vectors[,indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
indsX <- (20+3):(20+4)
UX <- eigOmega$vectors[,indsX]
L <- cbind(UX,LZ)
stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))
if (sum(stabCheck^2) > 1.0001*(20+2)) print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
X <- cbind(rep(1,length(x)),x)
Z <- B%*%LZ


## ----echo=TRUE,tidy=FALSE,size='footnotesize'-------------------------------------------------------------------------------
model<- 'model {
  for (i in 1:N) {
   y[i] ~ dnorm(mu[i],tau)   
   mu[i] <- X[i,] %*% beta[1:2] + Z[i,] %*% b[1:K]
  }
  for (k in 1:K){
    b[k] ~ dnorm(0,taub)
  }
  taub ~ dgamma(1.0E-3,1.0E-3)
  tau ~ dgamma(1.0E-3,1.0E-3)
  for (j in 1:2) {
    beta[j] ~ dnorm( 0  ,1.0E-3)
  }
}'


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
dati3 <- list(N=length(y),
                        K=ncol(Z),
                         y=y,X=X,Z=Z)


## ----echo=TRUE,cache=TRUE,cache.path='cache/mcmc2B',tidy=FALSE,results='hide',eval=FALSE------------------------------------
## fit_spline3 <- stan(file = 'spline3.stan', data = dati3,
##                     cores=parallel::detectCores())


## ----echo=FALSE,results='hide'----------------------------------------------------------------------------------------------
#save("fit_spline3",file="fit_spline3.Rdata")
load(file="fit_spline3.Rdata")


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
posterior=extract(fit_spline3)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
beta=posterior$beta
dim(beta)
b=posterior$b
dim(b)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
f=X%*% t(beta) + Z %*% t(b)
dim(f)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
pe=apply(f,1,median)
ci=apply(f,1,quantile,probs=c(0.025,0.5,0.975))


## ----echo=TRUE,eval=TRUE,out.height='0.6\\textheight',out.width='0.6\\textwidth'--------------------------------------------
plot(x,y,pch=20)
lines(x,pe,col="red")
matlines(x,t(ci),lwd=2,col="red",lty=c(2,1,2))


## ----echo=FALSE,out.width="0.6\\textwidth", out.height="0.4\\textwidth",fig.width=6,fig.height=4----------------------------
tu=read.table('trade.union.txt',header=TRUE)
x=tu$wage
y=tu$union.member
plot(x,y,pch="|")


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
tu=read.table('trade.union.txt',header=TRUE)
x=tu$wage
y=tu$union.member


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
knots=seq(min(x),max(x),length=22)[2:21]
X=cbind(1,x)
Z_K<-(abs(outer(x,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% 
                    (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
dati4 <- list(N=length(y),K=length(knots),
                         y=y,X=X,Z=Z)


## ----echo=TRUE,tidy=FALSE,eval=FALSE----------------------------------------------------------------------------------------
## model {
##   for (i in 1:N) {
##    y[i] ~ bernoulli_logit(X[i,] * beta[1:2] + Z[i,] * b[1:K]);
## }
## for (k in 1:K){
##   b[k] ~ normal(0,sigmab);
## }
## sigmab ~ normal(0,1000);
## for (j in 1:2) {
##   beta[j] ~ normal( 0  ,1000);
## }
## }


## ----echo=FALSE,results='hide'----------------------------------------------------------------------------------------------
#save("fit_spline4",file="fit_spline4.Rdata")
load(file="fit_spline4.Rdata")
posterior=extract(fit_spline4)


## ----echo=TRUE,cache=TRUE,cache.path='cache/mcmc3',tidy=FALSE,results='hide',eval=FALSE-------------------------------------
## fit_spline4 <- stan(file = 'spline4.stan', data = dati4,
##                     cores=parallel::detectCores())
## posterior=extract(fit_spline4)


## ----echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------
## launch_shinystan("fit_spline4.Rdata")


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
logitp=X%*% t(posterior$beta) + Z %*% t(posterior$b)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
p=exp(logitp)/(1+exp(logitp))


## ----echo=TRUE,eval=TRUE,out.height="0.7\\textheight",out.width="0.8\\textwidth",tidy=FALSE,fig.width=8,fig.height=7--------
par(mfrow=c(1,1))
plot(x,y,pch="|")
sl=sort.list(x)
x.sorted=x[sl]
p=p[sl,]
ci=apply(p,1,quantile,
         probs=c(0.005,0.025,0.5,0.975,0.995))
matlines(x.sorted,t(ci),lwd=2,
         col="red",lty=c(2,2,1,2,2))


## ----echo=TRUE,eval=FALSE,tidy=FALSE,size='scriptsize'----------------------------------------------------------------------
## library(coda)
## hpd.int=HPDinterval(as.mcmc(t(p)),prob=c(0.95))
## matlines(x.sorted,hpd.int,col="green",lwd=2,lty=1)
## hpd.int=HPDinterval(as.mcmc(t(p)),prob=c(0.99))
## matlines(x.sorted,hpd.int,col="green",lwd=2,lty=1)


## ----echo=FALSE,eval=TRUE,out.height="0.7\\textheight",out.width="0.8\\textwidth",fig.width=8,fig.height=7------------------
par(mfrow=c(1,1))
plot(x,y,pch="|")
ci=apply(p,1,quantile,probs=c(0.005,0.025,0.5,0.975,0.995))
matlines(x.sorted,t(ci),lwd=2,col="red",lty=c(2,2,1,2,2))
hpd.int=HPDinterval(as.mcmc(t(p)),prob=c(0.95))
matlines(x.sorted,hpd.int,col="green",lwd=2,lty=1)
hpd.int=HPDinterval(as.mcmc(t(p)),prob=c(0.99))
matlines(x.sorted,hpd.int,col="green",lwd=2,lty=1)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
x.new=seq(0,50,length=50)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
X.new=cbind(1,x.new)
Z_K<-(abs(outer(x.new,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% 
                    (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Z.new<-t(solve(sqrt.OMEGA_all,t(Z_K)))


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
logitp=X.new%*% t(posterior$beta) + Z.new %*% t(posterior$b)
p=exp(logitp)/(1+exp(logitp))


## ----echo=TRUE,eval=FALSE,tidy=FALSE----------------------------------------------------------------------------------------
## plot(x,y,pch="|")
## matlines(x.new,t(apply(p,1,quantile,
##                        probs=c(0.025,0.5,0.975))),
##          lwd=2,col="blue",lty=c(2,1,2))


## ----scatteronions,echo=FALSE,out.height="0.6\\textheight",out.width="0.7\\textwidth",fig.width=8,fig.height=6--------------
onions=read.table("onions.txt",header=TRUE)
ind=onions$location+1
colori=c("red","blue")
plot(onions$dens,onions$yield,pch=c(20,20)[onions$location+1],xlab="Density",ylab="Yield",col=colori[ind])
legend(150,250,pch=c(20,20),col=colori,legend=c("Purnong Landing","Virginia"))


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
y=onions$yield
x=(onions$dens-mean(onions$dens))/sd(onions$dens)
X=cbind(1,x)
knots=seq(min(x),max(x),length=12)[-c(1,12)]
Z_K<-(abs(outer(x,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))
dati5 <- list(N=length(y),
                         K=length(knots),
                         y=y,X=X,Z=Z)


## ----echo=TRUE,tidy=FALSE,size='tiny',eval=FALSE----------------------------------------------------------------------------
## model {
##   for (i in 1:N) {
##    y[i] ~ normal(X[i,] * beta[1:2] + Z[i,] * mspl[1:K] +
##             X[i,] * suspbeta[ind[i],1:2]' +
##             ZB[i,] * susp[ind[i],1:KB]',sigma);
## }
## for (k in 1:K){
##   mspl[k] ~ normal(0,sigmaspl);
## }
## for (h in 1:2){
##    suspbeta[h,1] ~ normal( 0  ,sigmau);
##    suspbeta[h,2] ~ normal( 0  ,sigmau);
##    for (k in 1:KB){
##     susp[h,k] ~ normal(0,sigmasusp);
##    }}
## for (j in 1:2) {
##   beta[j] ~ normal( 0  ,1000);
## }
## }


## ----echo=FALSE,results='hide'----------------------------------------------------------------------------------------------
#save("fit_spline5",file="fit_spline5.Rdata")
load(file="fit_spline5.Rdata")
posterior=extract(fit_spline5)


## ----echo=TRUE,cache=TRUE,cache.path='cache/mcmconions',tidy=FALSE,results='hide',eval=FALSE--------------------------------
## fit_spline5 <- stan(file = 'spline5.stan', data = dati5,
##                     cores=parallel::detectCores())
## posterior=extract(fit_spline5)


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------

xnew=seq(min(x),max(x),length=30)
Xnew=cbind(1,xnew)
Z_K<-(abs(outer(xnew,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% 
                    (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Znew<-t(solve(sqrt.OMEGA_all,t(Z_K)))
f=Xnew%*% t(posterior$beta) +  Znew %*% t(posterior$mspl) 
plot(x,y,pch=c(20,20)[ind],col=c("red","blue")[ind])
matlines(xnew,t(apply(f,1,quantile,probs=c(0.025,0.5,0.975))),
         col="black",lwd=2,lty=c(2,1,2))


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
y=onions$yield
x=(onions$dens-mean(onions$dens))/sd(onions$dens)
ind=onions$location+1
X=cbind(1,x)
knots=seq(min(x),max(x),length=12)[-c(1,12)]
Z_K<-(abs(outer(x,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
knotsB=seq(min(x),max(x),length=7)[-c(1,7)]
Z_K<-(abs(outer(x,knotsB,"-")))^3 
OMEGA_all<-(abs(outer(knotsB,knotsB,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% 
                    (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
ZB<-t(solve(sqrt.OMEGA_all,t(Z_K)))


## ----echo=TRUE,tidy=FALSE---------------------------------------------------------------------------------------------------
dati6 <- list(N=length(y),K=length(knots),
                         KB=length(knotsB),y=y,X=X,
                         Z=Z,ZB=ZB,ind=ind)


## ----echo=TRUE,tidy=FALSE,eval=FALSE,size='tiny'----------------------------------------------------------------------------
## model {
##   for (i in 1:N) {
##    y[i] ~ normal(X[i,] * beta[1:2] + Z[i,] * mspl[1:K] +
##             X[i,] * suspbeta[ind[i],1:2]' +
##             ZB[i,] * susp[ind[i],1:KB]',sigma);
## }
## for (k in 1:K){
##   mspl[k] ~ normal(0,sigmaspl);
## }
## for (h in 1:2){
##    suspbeta[h,1] ~ normal( 0  ,sigmau);
##    suspbeta[h,2] ~ normal( 0  ,sigmau);
##    for (k in 1:KB){
##     susp[h,k] ~ normal(0,sigmasusp);
##    }}
## for (j in 1:2) {
##   beta[j] ~ normal( 0  ,1000);
## }
## }


## ----echo=FALSE,results='hide'----------------------------------------------------------------------------------------------
#save("fit_spline6",file="fit_spline6.Rdata")
load(file="fit_spline6.Rdata")
posterior=extract(fit_spline6)


## ----echo=TRUE,cache=TRUE,cache.path='cache/mcmcB',tidy=FALSE,results='hide',eval=FALSE-------------------------------------
## fit_spline6 <- stan(file = 'spline6.stan', data = dati6,
##                     cores=parallel::detectCores())
## posterior=extract(fit_spline6)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------
str(posterior)


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
suspA1=posterior$suspbeta[,1,]
susp1=posterior$susp[,1,]


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
suspA2=posterior$suspbeta[,2,]
susp2=posterior$susp[,2,]


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
xnew=seq(min(x),max(x),length=30)
beta=posterior$beta


## ----echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------
Xnew=cbind(1,xnew)
Z_K<-(abs(outer(xnew,knots,"-")))^3 
OMEGA_all<-(abs(outer(knots,knots,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
Znew<-t(solve(sqrt.OMEGA_all,t(Z_K)))

Z_K<-(abs(outer(xnew,knotsB,"-")))^3 
OMEGA_all<-(abs(outer(knotsB,knotsB,"-")))^3 
svd.OMEGA_all<-svd(OMEGA_all) 
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 
ZBnew<-t(solve(sqrt.OMEGA_all,t(Z_K)))


## ----echo=TRUE,eval=FALSE,tidy=FALSE----------------------------------------------------------------------------------------
## spl=posterior$mspl
## f1=Xnew%*% t(beta) + Xnew%*%t(suspA1[,1:2])+
##     Znew %*% t(spl) + ZBnew %*% t(susp1)
## f2=Xnew%*% t(beta) + Xnew%*%t(suspA2[,1:2])+
##     Znew %*% t(spl) + ZBnew %*% t(susp2)
## plot(x,y,pch=c(20,20)[ind],col=c("red","blue")[ind])
## matlines(xnew,t(apply(f1,1,quantile,
##                       probs=c(0.025,0.5,0.975))),
##          col="red",lwd=2,lty=c(2,1,2))
## matlines(xnew,t(apply(f2,1,quantile,
##                       probs=c(0.025,0.5,0.975))),
##          col="blue",lwd=2,lty=c(2,1,2))


## ----echo=FALSE,eval=TRUE,out.height="0.8\\textheight",out.width="0.8\\textwidth",tidy=FALSE--------------------------------
spl=posterior$mspl
f1=Xnew%*% t(beta) + Xnew%*%t(suspA1[,1:2])+ 
      Znew %*% t(spl) + ZBnew %*% t(susp1)
f2=Xnew%*% t(beta) + Xnew%*%t(suspA2[,1:2])+ 
      Znew %*% t(spl) + ZBnew %*% t(susp2)
plot(x,y,pch=c(20,20)[ind],col=c("red","blue")[ind])
matlines(xnew,t(apply(f1,1,quantile,probs=c(0.025,0.5,0.975))),
         col="red",lwd=2,lty=c(2,1,2))
matlines(xnew,t(apply(f2,1,quantile,probs=c(0.025,0.5,0.975))),
         col="blue",lwd=2,lty=c(2,1,2))


## ----echo=FALSE,out.height='0.6\\textwidth',out.width='0.7\\textheight'-----------------------------------------------------
pig=read.table("pig.weights.txt",header=TRUE)
pig2=reshape(pig,idvar="id.num",v.names="weight",timevar="num.weeks",direction="wide")
matplot(1:9,t(pig2[,-1]),pch=1,col="black",xlab="week",ylab="weight")
matlines(1:9,t(pig2[,-1]),col="black",lty=1)
dim(pig2)

