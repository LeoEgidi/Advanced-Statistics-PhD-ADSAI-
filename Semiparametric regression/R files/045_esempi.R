## 
##
##


library(knitr)
options(width=50)
mypref=function() opts_chunk$set(comment=NA, fig.width=6, fig.height=4,out.width='0.7\\textwidth',echo=FALSE,results='hide',fig.path="figure/035-esempi-",tidy=FALSE)
mypref()
library(xtable)
library(splines)
library(mgcv)
library(gamair)
library(SemiPar)
## ------------------------------------------------------------------------
cmb=read.table("wmap.dat",header=TRUE)
lidar=read.table("lidar.dat",header=TRUE)
lidar=lidar[sort.list(lidar$range),]
bpd=read.table("bpd.dat",header=TRUE)




############################################################################
############################################################################


## 
library(gamair)
data(hubble)
head(hubble)

## 
plot(hubble$x,hubble$y,pch=20)

## 
fit1=lm(y~x-1,data=hubble)
coef(fit1)
confint(fit1)

## 
plot(hubble$x,hubble$y,pch=20)
abline(fit1)

## 
pastmeas=data.frame(matrix(c(2016,73,1.75,2013,67.8,0.77,2012,69.32,0.80,2010,70.4,1.4,2010,71,2.5,2009,70.1,1.3,2009,71.9,2.7,2007,70.4,1.6,2006,77.6,14.9,2001,72,8,1958,75,NA),byrow=TRUE,ncol=3))
names(pastmeas)=c("year","const","error")

## 
fitGam=gam(y~s(x),data=hubble)
summary(fitGam)
plot(fitGam)

## 
plot(fitGam,all.terms=TRUE)

##
fitGamLin=gam(y~x,data=hubble)
AIC(fitGam,fitGamLin)
anova(fitGamLin,fitGam,test="F")

##
gam.check(fitGam)

## 
fitGam=gam(y~s(x),data=hubble,family=quasi(var=mu))
summary(fitGam)
plot(fitGam)

## 
plot(fitGam,all.terms=TRUE)

## 
gam.check(fitGam)

## 
library(sm)
sm.regression(hubble$x, hubble$y, 
              model = "linear")

## 
sm.regression(hubble$x, hubble$y, model = "linear")

############################################################################
############################################################################
##
##
## 
data(age.income)
age.income[1,]
plot(age.income$age,
     age.income$log.income)


##
data(age.income)
age.income[1,]
plot(age.income$age,
     age.income$log.income)
fit=gam(log.income~age,data=age.income)
curve(predict(fit,newdata=data.frame(age=x),type="response"),ad=TRUE,lwd=2)
fit=gam(log.income~s(age),data=age.income)
curve(predict(fit,newdata=data.frame(age=x),type="response"),ad=TRUE,lwd=2)

##
library(sm)
sm.regression(age.income$age, 
              age.income$log.income, 
              se = TRUE)

##
sm.regression(age.income$age, age.income$log.income, se = TRUE)

## 
library(sm)
sm.regression(age.income$age, 
              age.income$log.income, 
              model = "linear")

## 
sm.regression(age.income$age, age.income$log.income, model = "linear")

#######################################################################################
########################################################################################
### 
### Chicago pollution

data(chicago)
head(chicago)
chicago$data=seq(as.Date("1987/1/1"), as.Date("2000/12/31"), "days")
layout(matrix(c(1,1,1,1,2,3,4,5),byrow=TRUE,ncol=4))
par(mar=c(5,4,0,0))
plot(chicago$data,chicago$death,pch=".")
plot(chicago$pm10median,chicago$death,pch=".")
plot(chicago$so2median,chicago$death,pch=".")
plot(chicago$o3median,chicago$death,pch=".")
plot(chicago$tmpd,chicago$death,pch=".")

## 
options(width=90)
fit1=gam(death~s(time,bs="cr",k=200)+
           pm10median+so2median+o3median+
           tmpd,
         data=chicago,family=poisson)
par(mar=c(5,4,0,0))
plot(fit1)
summary(fit1)

## 
gam.check(fit1)

## 
fit1$model[sort.list(resid(fit1),decreasing=TRUE)[1:4],]
chicago[3110:3120,]
par(mfrow=c(1,2))
plot(fit1$model$tmpd,resid(fit1))
plot(fit1$model$o3median,resid(fit1))

chicago$death[3111:3125]
hist(chicago$death)
rug(chicago$death,col="red")

## try a more flexible model
fit1a=gam(death~s(time,bs="cr",k=200)+
           s(pm10median)+s(so2median)+s(o3median)+
           s(tmpd),
         data=chicago,family=poisson)
par(mar=c(5,4,0,0))
plot(fit1a)
summary(fit1a)
gam.check(fit1a)



## try lags
meanlag3=function(x){
  n=length(x)
  c(NA,NA,NA,(x[1:(n-3)]+x[2:(n-2)]+x[3:(n-1)]+x[4:n])/4)
}
chicago$pm10medianLag=meanlag3(chicago$pm10median)
chicago$pm25medianLag=meanlag3(chicago$pm25median)
chicago$o3medianLag=meanlag3(chicago$o3median)
chicago$so2medianLag=meanlag3(chicago$so2median)
chicago$tmpdLag=meanlag3(chicago$tmpd)

## 
layout(matrix(c(1,2,3,4,5,6,7,8),byrow=TRUE,ncol=4))
par(mar=c(5,4,0,0))
plot(chicago$pm10median,chicago$death,pch=".")
plot(chicago$so2median,chicago$death,pch=".")
plot(chicago$o3median,chicago$death,pch=".")
plot(chicago$tmpd,chicago$death,pch=".")
plot(chicago$pm10medianLag,chicago$death,pch=".")
plot(chicago$so2medianLag,chicago$death,pch=".")
plot(chicago$o3medianLag,chicago$death,pch=".")
plot(chicago$tmpdLag,chicago$death,pch=".")

##
chicago$data=seq(as.Date("1987/1/1"), as.Date("2000/12/31"), "days")
chicago$yearday=as.POSIXlt(chicago$data)$yday-366/2
chicago$weekday=as.POSIXlt(chicago$data)$wday
chicago$year=as.POSIXlt(chicago$data)$year

## 
layout(matrix(c(1,2,3,4,5,6,7,8),byrow=TRUE,ncol=4))
par(mar=c(5,4,0,0),mfrow=c(1,3))
plot(chicago$yearday,chicago$death,pch=".")
plot(as.factor(chicago$weekday),chicago$death,pch=".")
plot(as.factor(chicago$year),chicago$death,pch=".")

## some alternative models
fit2=gam(death~s(time,bs="cr",k=200)+
           pm10medianLag+so2medianLag+o3medianLag+
           tmpdLag,
         data=chicago,family=poisson)
fit3=gam(death~s(yearday,bs="cr")+
           pm10medianLag+so2medianLag+o3medianLag+
           s(tmpdLag),
         data=chicago,family=poisson)
fit4=gam(death~s(yearday,bs="cr")+
           so2medianLag+s(pm10medianLag)+
           s(o3medianLag,tmpdLag,k=40),
         data=chicago,family=poisson)
fit5=gam(death~s(yearday,bs="cr")+
           so2medianLag+s(pm10medianLag,bs="cr",k=6)+
           te(o3medianLag,tmpdLag,k=8) ,
         data=chicago,family=poisson)
fit6=gam(death~s(time,bs="cr",k=200) +
           te(o3medianLag,tmpdLag,pm10medianLag,k=c(8,8,6)),
         data=chicago,family=poisson)
fit7=gam(death~s(time,bs="cr",k=200) +
           te(o3medianLag,tmpdLag,k=8)+s(pm10medianLag,bs="cr",k=6),
         data=chicago,family=poisson)

## 
vis.gam(fit4,c("o3medianLag","tmpdLag"),theta=45,too.far=0.07)
vis.gam(fit4,c("o3medianLag","tmpdLag"),plot.type="contour",too.far=0.07)
points(fit4$model$o3medianLag,fit4$model$tmpdLag)


## 
plot(as.factor(na.omit(chicago[,-c(3,10)])$year),residuals(fit4))

##############################################################################
##############################################################################
##
##
##



## 
data(co2s)
plot(co2s$c.month,co2s$co2,type="l")


##
b=gam(co2~s(c.month,k=300,bs="cr"),data=co2s)
n=nrow(co2s)
pd <- data.frame(c.month=1:(n+36))
fv <- predict(b,pd,se=TRUE)
plot(pd$c.month,fv$fit,type="l")
points(co2s$c.month,co2s$co2,pch=".")
lines(pd$c.month,fv$fit+2*fv$se,col=2)
lines(pd$c.month,fv$fit-2*fv$se,col=2)


##
b2 <- gam(co2~s(month,bs="cc")+
            s(c.month,bs="cr",k=300),
          data=co2s,
          knots=list(month=seq(1,13,length=10)))

pd2 <- data.frame(c.month=1:(n+36),
                  month=rep(1:12,length.out=n+36))
fv <- predict(b2,pd2,se=TRUE)
plot(pd$c.month,fv$fit,type="l")
lines(pd$c.month,fv$fit+2*fv$se,col=2)
lines(pd$c.month,fv$fit-2*fv$se,col=2)


#############################################################################
#############################################################################
##
##
##
library(SemiPar)
library(maps)
data(ustemp)
attach(ustemp)
grey.levs <- min.temp+20
col.vec <- paste("grey",as.character(grey.levs),sep="")
map("usa")
points(-longitude,latitude,
     col=col.vec,pch=16,cex=3,xlim=c(-130,-60))
text(-longitude,latitude,as.character(city))
detach(ustemp)

##
fitSep=gam(min.temp~s(longitude)+s(latitude),data=ustemp)
par(mfrow=c(1,2),mar=c(5,4,0,0))
plot(fitSep)

##
par(mfrow=c(1,2),mar=c(5,4,0,0))
vis.gam(fitSep,type="response")
vis.gam(fitSep,type="response",
        plot.type="contour",color="topo")

##
fitJoint=gam(min.temp~s(longitude,latitude),data=ustemp)
par(mfrow=c(1,2),mar=c(5,4,0,0))
vis.gam(fitJoint,type="response")
vis.gam(fitJoint,type="response",
        plot.type="contour",color="topo")

##
par(mfrow=c(1,2),mar=c(5,4,0,0))
vis.gam(fitSep,type="response",
        plot.type="contour",color="topo")
vis.gam(fitJoint,type="response",
        plot.type="contour",color="topo")

##
par(mfrow=c(1,2),mar=c(5,4,0,0))
vis.gam(fitSep,type="response",
        color="topo",theta=30)
vis.gam(fitJoint,type="response",
        color="topo",theta=30)

## 
par(mfrow=c(1,2),mar=c(5,4,0,0))
fitK=sm.regression(cbind(ustemp$longitude,ustemp$latitude),
                   ustemp$min.temp)
fitK=sm.regression(cbind(ustemp$longitude,ustemp$latitude),
                   ustemp$min.temp,
                   display="image",ngrid=100)

## 
## library(mgcv)
## library(rworldmap)
## newmap <- getMap(resolution = "low")
## 
## lonseq=sort(-seq(min(ustemp$longitude),max(ustemp$longitude),length=100))
## latseq=seq(min(ustemp$latitude),max(ustemp$latitude),length=100)
## z=outer(lonseq,latseq,FUN=function(x,y) predict(fitSep,newdata=data.frame(longitude=-x,latitude=y,type="response")))
## plot(newmap,xlim = range(lonseq),ylim = range(latseq),asp = 1)
## contour(lonseq,latseq,z,add=TRUE,col="red")
## points(-ustemp$longitude,ustemp$latitude,
##      col=col.vec,pch=16,cex=3,xlim=c(-130,-60))
## text(-ustemp$longitude,ustemp$latitude,as.character(city))

############################################################################
############################################################################
##
##
##

data(mack)
data(coast)
plot(mack$lon,mack$lat,asp=1)
lines(coast$lon,coast$lat,
      col="blue")


## 
symbols(mack$lon,mack$lat,circles=mack$egg.count,inches=0.1,asp=1)
lines(coast$lon,coast$lat,col="blue")



## different bases and maximum smooth
x=sort(runif(100,0,1))
m=3*x+sin(2*pi*x)
y=m+rnorm(100,0,0.5*sd(m))
plot(x,y)

fit=gam(y~s(x,bs="tp"))
fit1=gam(y~s(x,bs="tp"),sp=10^10)
fit2=gam(y~s(x,bs="ts"))
fit3=gam(y~s(x,bs="ts"),sp=10^10)
plot(x,y)
curve(predict(fit,
              newdata=data.frame(x=x)),
      add=TRUE)
curve(predict(fit1,
              newdata=data.frame(x=x)),
      add=TRUE,col="red")
curve(predict(fit2,
              newdata=data.frame(x=x)),
      add=TRUE,col="blue")
curve(predict(fit3,newdata=data.frame(x=x)),
      add=TRUE,col="green")

## basic model
mack$log.net.area <- log(mack$net.area)
gm <- gam(egg.count ~ s(lon,lat,bs="ts")+
            s(I(b.depth^.5),bs="ts")+
            s(c.dist,bs="ts")+
            s(salinity,bs="ts")+
            s(temp.surf,bs="ts")+
            s(temp.20m,bs="ts")+
            offset(log.net.area),
          data=mack,
          family=poisson,
          scale=-1,gamma=1.4)

## 
par(mfrow=c(2,3),mar=c(5,4,0,0))
plot(gm,scale=0)

## 
par(mfrow=c(2,3))
gam.check(gm)

summary(gm)

## model suggested
gm2a <- gam(egg.count ~ s(lon,lat,bs="ts",k=100)+
             s(I(b.depth^.5),bs="ts")+
             s(c.dist,bs="ts")+
             s(temp.surf,bs="ts")+
             s(temp.20m,bs="ts")+
             offset(log.net.area),
           data=mack,
           family=poisson,
           scale=-1,gamma=1.4)

gm2 <- gam(egg.count ~ s(lon,lat,bs="ts")+
             s(I(b.depth^.5),bs="ts")+
             s(c.dist,bs="ts")+
             s(temp.surf,bs="ts")+
             s(temp.20m,bs="ts")+
             offset(log.net.area),
           data=mack,
           family=poisson,
           scale=-1,gamma=1.4)

## 
gm3 <- gam(egg.count ~ s(lon,lat,bs="ts")+
             s(I(b.depth^.5),bs="ts")+
             s(c.dist,bs="ts")+
             s(temp.20m,bs="ts")+
             offset(log.net.area),
           data=mack,
           family=poisson,
           scale=-1,gamma=1.4)

## 
summary(gm3)

## 
gm4<-gam(egg.count ~ s(lon,lat,bs="ts",k=40) +
           s(I(b.depth^.5),bs="ts") + 
           s(c.dist,bs="ts") +
           s(temp.20m,bs="ts") +
           offset(log.net.area),
         data=mack,
         family=negbin(1),
         control=gam.control(maxit=100),gamma=1.4)


model=gm2

## predictions
data(mackp)
mackp$log.net.area <- 0*mackp$lon # make offset column
lon<-seq(-15,-1,1/4);lat<-seq(44,58,1/4)
zz<-array(NA,57*57)
zz[mackp$area.index]<-predict(model,mackp)
par(mar=c(5,4,0,0))
image(lon,lat,matrix(zz,57,57),col=gray(0:32/32),
cex.lab=1.5,cex.axis=1.4)
contour(lon,lat,matrix(zz,57,57),add=TRUE)
lines(coast$lon,coast$lat,col=1)


library(MASS)
br1 <- mvrnorm(n=1000,coef(model),model$Vp)
Xp <- predict(model,newdata=mackp,type="lpmatrix")
mean.eggs1 <- colMeans(exp(Xp%*%t(br1)))
hist(mean.eggs1)

f<-fitted(model)
form<-egg.count~offset(log.net.area)+s(lon,lat,bs="ts",k=100)+
  s(I(b.depth^.5),bs="ts")+s(temp.20m,bs="ts")
mack.bs <- mack
n <- nrow(mack)
br <- matrix(0,0,length(coef(model)))
for (i in 1:19) { 
  e <- rpois(rep(1,n),f) - f
  y <- round(f+e*model$sig2^.5)
  y[y<0] <- 0
  mack.bs$egg.count <- y
  sp <- gam(form,data=mack.bs,family=poisson,scale=-1,gamma=1.4)$sp
  b <- gam(form,data=mack,family=poisson,sp=sp,scale=-1)
  br <- rbind(br,mvrnorm(n=100,coef(b),b$Vp))
}
br <- rbind(br,mvrnorm(n=100,coef(gm1a),gm1a$Vp))
mean.eggs <- colMeans(exp(Xp%*%t(br)))
hist(mean.eggs)


##########################################################################
##########################################################################
## 
## 
## 
data(brain)
brain=brain[brain$medFPQ>5e-3,]
brain$medFPQ.s=(brain$medFPQ-min(brain$medFPQ))/(max(brain$medFPQ)-min(brain$medFPQ))
brain$logmedFPQ=log(brain$medFPQ)
brain$logmedFPQ.s=(brain$logmedFPQ-min(brain$logmedFPQ))/(max(brain$logmedFPQ)-min(brain$logmedFPQ))
plot(brain$Y,brain$X,col=gray((brain$logmedFPQ.s)),pch=15,type="n")
rect(brain$Y,brain$X,brain$Y+1,brain$X+1,col=gray((brain$logmedFPQ.s)),border = gray((brain$logmedFPQ.s)))
#plot(brain$Y,brain$X,col=gray((brain$logmedFPQ.s)),pch=15,type="n")
#rect(brain$Y,brain$X,brain$Y+1,brain$X+1,col=gray((brain$medFPQ.s)),border = gray((brain$medFPQ.s)))

## clearly non gaussian
m0=gam(medFPQ~s(Y,X,k=100),data=brain)
par(mfrow=c(1,3),mar=c(5,4,0,0))
qqnorm(residuals(m0))
qqline(residuals(m0))
plot(fitted(m0),residuals(m0))
plot(log(fitted(m0)),log(residuals(m0)^2))

## study the relationship between variance and mean
lm(log(residuals(m0)^2)~log(fitted(m0)))

## alternative 1: transformation
m1=gam(medFPQ^.25~s(Y,X,k=100),data=brain)

## alternative 2: gamma model
m2=gam(medFPQ~s(Y,X,k=100),data=brain,
       family=Gamma(link=log))

## both imply bias on the response scale
mean(fitted(m1)^4);mean(fitted(m2));mean(brain$medFPQ)



## 
vis.gam(m2,plot.type="contour",too.far=0.02,n.grid=100)

## try a simpler additive structure
m3=gam(medFPQ~s(Y,k=30)+s(X,k=30),data=brain,
       family=Gamma(link=log))
m3$gcv.ubre
anova(m2,m3,test="F")

## 
vis.gam(m3,plot.type="contour",too.far=0.02,n.grid=100)

## testing more apprpriate using a nested structure
m4=gam(medFPQ~s(Y,k=30)+s(X,k=30)+s(Y,X,k=100),data=brain,
       family=Gamma(link=log))
m4$gcv.ubre
anova(m4,m3,test="F")

## 
vis.gam(m3,plot.type="contour",too.far=0.02,n.grid=100)

## checking for simmetry
brain$Xc=abs(brain$X-64.5)
brain$R=1*(brain$X<64.5)
ms1=gam(medFPQ~s(Y,Xc,k=100),data=brain,
       family=Gamma(link=log),method="gm")
ms2=gam(medFPQ~s(Y,Xc,k=100)+s(Y,Xc,k=100,by=R),data=brain,
       family=Gamma(link=log),method=gm)
anova(ms1,ms2,test="F")
anova(ms2)

## 
vis.gam(ms1,plot.type="contour",too.far=0.04,n.grid=100)
vis.gam(ms2,plot.type="contour",too.far=0.04,n.grid=100)

vis.gam(ms1,plot.type="contour",view=c("Xc","Y"),too.far=.03,
        color="gray",n.grid=60,zlim=c(-1,2),main="both sides")
vis.gam(ms2,plot.type="contour",view=c("Xc","Y"),
        cond=list(R=0),too.far=.03,color="gray",n.grid=60,
        zlim=c(-1,2),main="left side")
vis.gam(ms2,plot.type="contour",view=c("Xc","Y"),
        cond=list(R=1),too.far=.03,color="gray",n.grid=60,
        zlim=c(-1,2),main="right side")


## comparing two surfaces (5.2.5)






#############################################################################
#############################################################################
##
##
##
##
## 


library(gamair)
data(bpd)
bpd[1,]
plot(bpd$birthweight,bpd$BPD,pch="|")
fit=gam(BPD~birthweight,data=bpd,family=binomial)
fitS=gam(BPD~s(birthweight),data=bpd,family=binomial)
curve(predict(fit,newdata=data.frame(birthweight=x),
              type="response"),ad=TRUE,lwd=2)
curve(predict(fitS,
              newdata=data.frame(birthweight=x),
              type="response"),ad=TRUE,lwd=2,col="red")

## 
## 
library(sm)
M=mean(bpd$birthweight)
S=sd(bpd$birthweight)
bpd$birthweight.st=(bpd$birthweight-M)/S
sm.binomial(bpd$birthweight.st,
            bpd$BPD,
            h=0.5)
sm.binomial(bpd$birthweight.st,
            bpd$BPD,
            h=0.25,add=TRUE,col="red",pch=NA)

## 
par(mar=c(5,4,0.3,0.3))
nba=read.table("Esempi/TalentPerformance/NBA.csv",sep=",")
names(nba)=c("t","Performance")
plot(nba$t,nba$Performance)
vlim=0.5380863039399624
abline(v=vlim)
fit1=lm(Performance~t,data=subset(nba,t<vlim))
abline(fit1,lty=2)
lines(c(min(nba$t),vlim-0.01),coef(fit1)[1]+coef(fit1)[2]*c(min(nba$t),vlim-0.01),lwd=2,col="red")
fit2=lm(Performance~t,data=subset(nba,t>vlim))
abline(fit2,lty=2)
lines(c(max(nba$t),vlim+0.01),coef(fit2)[1]+coef(fit2)[2]*c(max(nba$t),vlim+0.01),lwd=2,col="red")
library(mgcv)
fitS=gam(Performance~s(t),data=nba)
#plot(fitS,add=TRUE)
curve(predict(fitS,data.frame(t=x)),from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2)
curve(predict(fitS,data.frame(t=x))
      -1.96*predict(fitS,data.frame(t=x),se.fit=TRUE)$se,
      from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2,lty=2)
curve(predict(fitS,data.frame(t=x))
      +1.96*predict(fitS,data.frame(t=x),se.fit=TRUE)$se,
      from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2,lty=2)

## 
par(mar=c(5,4,0.3,0.3))
nba=read.table("Esempi/TalentPerformance/NBA.csv",sep=",")
names(nba)=c("t","Performance")
plot(nba$t,nba$Performance)
vlim=0.5380863039399624
abline(v=vlim)
fit1=lm(Performance~t,data=subset(nba,t<vlim))
abline(fit1,lty=2)
lines(c(min(nba$t),vlim-0.01),coef(fit1)[1]+coef(fit1)[2]*c(min(nba$t),vlim-0.01),lwd=2,col="red")
fit2=lm(Performance~t,data=subset(nba,t>vlim))
abline(fit2,lty=2)
lines(c(max(nba$t),vlim+0.01),coef(fit2)[1]+coef(fit2)[2]*c(max(nba$t),vlim+0.01),lwd=2,col="red")
fitS=smooth.spline(nba$t,nba$Performance)
curve(predict(fitS,x)$y,from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2)

B=1000
xseq=seq(min(nba$t),max(nba$t),length=100)
yt=predict(fitS,nba$t)$y
res=nba$Performance-yt
fun=funD=matrix(NA,B,length(xseq))
for (b in 1:B){
  yb=yt+sample(res,length(res),replace=TRUE)
  fitb=smooth.spline(nba$t,yb)
  fun[b,]=predict(fitb,xseq)$y
  funD[b,]=predict(fitb,xseq,deriv=1)$y
}

par(mar=c(5,4,0.3,0.3),mfrow=c(2,1))
plot(nba$t,nba$Performance)
vlim=0.5380863039399624
abline(v=vlim)
fit1=lm(Performance~t,data=subset(nba,t<vlim))
abline(fit1,lty=2)
lines(c(min(nba$t),vlim-0.01),coef(fit1)[1]+coef(fit1)[2]*c(min(nba$t),vlim-0.01),lwd=2,col="red")
fit2=lm(Performance~t,data=subset(nba,t>vlim))
abline(fit2,lty=2)
lines(c(max(nba$t),vlim+0.01),coef(fit2)[1]+coef(fit2)[2]*c(max(nba$t),vlim+0.01),lwd=2,col="red")
fitS=smooth.spline(nba$t,nba$Performance)
curve(predict(fitS,x)$y,from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2)
for (i in 1:B) lines(xseq,fun[i,],col=gray(0.7))

plot(xseq,0*xseq,ylim=range(funD),type="n")
curve(predict(fitS,x,deriv=1)$y,from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2)
for (i in 1:B) lines(xseq,funD[i,],col=gray(0.7))
abline(h=0)

## 
par(mar=c(5,4,0.3,0.3))
nba=read.table("Esempi/TalentPerformance/NBA33.csv",sep=",")
names(nba)=c("t","Performance")
plot(nba$t,nba$Performance)
vlim=0.4404029692470837
abline(v=vlim)
fit1=lm(Performance~t,data=subset(nba,t<vlim))
abline(fit1,lty=2)
lines(c(min(nba$t),vlim-0.01),coef(fit1)[1]+coef(fit1)[2]*c(min(nba$t),vlim-0.01),lwd=2,col="red")
fit2=lm(Performance~t,data=subset(nba,t>vlim))
abline(fit2,lty=2)
lines(c(max(nba$t),vlim+0.01),coef(fit2)[1]+coef(fit2)[2]*c(max(nba$t),vlim+0.01),lwd=2,col="red")
fitS=gam(Performance~s(t),data=nba)
#plot(fitS,add=TRUE)
curve(predict(fitS,data.frame(t=x)),from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2)
curve(predict(fitS,data.frame(t=x))
      -1.96*predict(fitS,data.frame(t=x),se.fit=TRUE)$se,
      from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2,lty=2)
curve(predict(fitS,data.frame(t=x))
      +1.96*predict(fitS,data.frame(t=x),se.fit=TRUE)$se,
      from=min(nba$t),to=max(nba$t),add=TRUE,col="blue",lwd=2,lty=2)

## 
par(mar=c(5,4,0.3,0.3))
mlb=read.table("Esempi/TalentPerformance/MLB.csv",sep=",")
names(mlb)=c("t","Performance")
plot(mlb$t,mlb$Performance)
fit3=lm(Performance~t,data=subset(mlb))
abline(fit3,lty=2)
lines(range(mlb$t),coef(fit3)[1]+coef(fit3)[2]*range(mlb$t),lwd=2,col="red")
fitS=gam(Performance~s(t),data=mlb)
#plot(fitS,add=TRUE)
curve(predict(fitS,data.frame(t=x)),from=min(mlb$t),to=max(mlb$t),add=TRUE,col="blue",lwd=2)
curve(predict(fitS,data.frame(t=x))
      -1.96*predict(fitS,data.frame(t=x),se.fit=TRUE)$se,
      from=min(mlb$t),to=max(mlb$t),add=TRUE,col="blue",lwd=2,lty=2)
curve(predict(fitS,data.frame(t=x))
      +1.96*predict(fitS,data.frame(t=x),se.fit=TRUE)$se,
      from=min(mlb$t),to=max(mlb$t),add=TRUE,col="blue",lwd=2,lty=2)

