## Bayesian analysis to estimate D, a parameter describing NFDS
#module load R/4.2.2

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)
library(RColorBrewer)
library(lme4)

dat<-read.table("Tcris_master_32.csv",sep=",",header=TRUE)

hv<-which(dat$location=="HV")

hvs<-tapply(X=dat$striped[hv],INDEX=dat$year[hv],sum)
hvg<-tapply(X=dat$unstriped[hv],INDEX=dat$year[hv],sum)
hvp<-hvs/(hvg+hvs)
yr<-as.numeric(names(hvp))

pdf("HVStripe.pdf",width=5,height=4)
par(mar=c(4.5,4.5,.5,.5))
plot(yr[-c(1:2)],hvp[-c(1:2)],type='b',pch=19,xlab="Year",ylab="Stripe frequency",cex.lab=1.3)
dev.off()

## 2 or more samples
ids<-names(which(table(dat$location) >= 2))
N<-length(ids)
obs<-rep(0,N)
for(i in 1:N){
    yrs<-unique(dat$year[which(dat$location==ids[i])])
    J<-length(yrs)
    if(J >=2){
        obs[i]<-sum(yrs[-1]-yrs[-J]==1)
    }
}    

## focus on sites with 10 or more pairs
keep<-which(obs >= 10)
obs[keep]
# [1] 10 22 15 15 10 13 14 17 12 12
ids[keep]
# [1] "FH"   "HV"   "L"    "M"    "MBOX" "OG"   "OUT"  "PR"   "SC"   "VP"
kids<-ids[keep]
J<-length(keep)
Y0<-vector("list",J)
N0<-vector("list",J)
Y1<-vector("list",J)
N1<-vector("list",J)
Years<-vector("list",J)
Yr0<-vector("list",J)
P<-vector("list",J)
for(j in 1:J){
	site<-which(dat$location==kids[j])
	nS<-tapply(X=dat$striped[site],INDEX=dat$year[site],sum)
	nG<-tapply(X=dat$unstriped[site],INDEX=dat$year[site],sum)
	N<-nS+nG
	yrs<-as.numeric(names(nS))
	Years[[j]]<-yrs
	P[[j]]<-nS/N
	cnt<-1
	for(k in 1:(length(yrs)-1)){
		if((yrs[k+1]-yrs[k]) == 1){
			Yr0[[j]][cnt]<-yrs[k]
			Y0[[j]][cnt]<-nS[k]
			N0[[j]][cnt]<-N[k]
			Y1[[j]][cnt]<-nS[k+1]
			N1[[j]][cnt]<-N[k+1]
			cnt<-cnt+1
		}
	}
}

### plot of all time series
allYears<-1990:2023
Nyrs<-length(allYears)
Pmat<-matrix(NA,nrow=10,ncol=Nyrs)
for(j in 1:10){
	xx<-which(allYears %in% Years[[j]])
	Pmat[j,xx]<-P[[j]]
}

pdf("StripeTimeSeries.pdf",width=9,height=4.5)
par(mar=c(4.5,4.5,1,1))
cs<-alpha(brewer.pal(n=10,"Paired"),.7)
plot(allYears,Pmat[1,],ylim=c(0,1),col=cs[1],type='b',pch=19,xlab="Year",ylab="Stripe frequency",cex.lab=1.3,cex.axis=1.1)
for(j in 2:10){
	points(allYears,Pmat[j,],col=cs[j],type='b',pch=19)
}
legend(1990,.6,ids[keep],pch=19,col=cs,bty='n',ncol=3)
dev.off()

pdf("StripeTimeSeries2000.pdf",width=9,height=4.5)
layout(matrix(1:2,nrow=1,ncol=2),widths=c(7.5,1.5),heights=4.5)
par(mar=c(4.5,4.5,1,1))
cs<-alpha(brewer.pal(n=10,"Paired"),.7)
plot(allYears[11:Nyrs],Pmat[1,11:Nyrs],ylim=c(0,1),col=cs[1],type='b',pch=19,xlab="Year",ylab="Stripe frequency",cex.lab=1.3,cex.axis=1.1)
for(j in 2:10){
	points(allYears[11:Nyrs],Pmat[j,11:Nyrs],col=cs[j],type='b',pch=19)
}

par(mar=c(0,0,0,0))
plot(0:1,0:1,type='n',xlab="",ylab="",axes=FALSE)
legend(.1,.9,ids[keep],pch=19,col=cs,bty='n')
dev.off()

### bayes ###

blk<-rep(1:10,unlist(lapply(X=N1,length)))
Y0<-unlist(Y0)
Y1<-unlist(Y1)
N0<-unlist(N0)
N1<-unlist(N1)

dat<-list(Y0=Y0,Y1=Y1,N0=N0,N1=N1,blk=blk,nn=length(Y0),J=10)
fit<-stan("hmodel.stan",data=dat)
oo<-extract(fit,"beta")
best<-apply(oo[[1]],2,quantile,probs=c(0.5,0.05,0.95))
oo<-extract(fit,"alpha")
aest<-apply(oo[[1]],2,quantile,probs=c(0.5,0.05,0.95))
oo<-extract(fit,"dp")
dpest<-apply(oo[[1]],2,median)
oo<-extract(fit,"p0")
pest<-apply(oo[[1]],2,median)


pdf("StripeChangeNFDS.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
cs<-alpha(brewer.pal(n=10,"Paired"),.7)
plot(pest,dpest,pch=19,col=cs[blk],xlab="Stripe frequency",ylab="Change in frequency",cex.lab=1.2)
for(j in 1:10){
	abline(a=aest[1,j],b=best[1,j],col=cs[j])
}
legend(.335,.9,ids[keep],pch=19,col=cs,bty='n',ncol=3)
abline(h=0)
dev.off()

## compare simple lm
p0<-Y0/N0
p1<-Y1/N1
dp<-p1-p0

cor(dp,dpest)
#[1]  0.9404365

a<-rep(NA,10)
b<-rep(NA,10)
for(j in 1:10){
	pop<-which(blk==j)
	o<-lm(dp[pop] ~ p0[pop])
	oo<-summary(o)
	a[j]<-oo$coefficients[1,1]
	b[j]<-oo$coefficients[2,1]
}

pbar<-tapply(X=p0,INDEX=blk,mean)
eq<--a/b
pdf("nfdsEquil.pdf",width=4.5,height=4.5)
par(mar=c(4.5,4.5,1,1))
plot(pbar,eq,xlab="Mean stripe frequency",ylab="Equilibrium stripe frequency",pch=19)
abline(a=0,b=1)
dev.off()
cor.test(pbar,eq)

#	Pearson's product-moment correlation
#
#data:  pbar and eq
#t = 42.891, df = 8, p-value = 9.627e-11
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9904991 0.9995070
#sample estimates:
#      cor 
#0.9978327 

## first year 2000+ by block
fy<-c(1,12,34,48,63,74,86,101,117,130)


dpdf<-cbind(D=best[1,],eq,p0=pest[fy],yr=unlist(Yr0)[fy])
write.table(dpdf,file="dpdf.txt",row.names=FALSE,col.names=TRUE)

save(list=ls(),file="nfds.rdat")

## dampening analysis
a<-which(unlist(Yr0) > 1995) ## get rid of the few really old years
## first year = 2000, set to 0
D<-list(dp=abs(dpest)[a],Yr=unlist(Yr0)[a]-2000,bl=blk[a],N=length(a),J=10)

## simple linear
o<-lm(D$dp ~ D$Yr)
summary(o)

## linear mixed model
o<-lmer(D$dp ~ D$Yr + (1 | D$bl))
o
## hierarchical Bayesian
fitd<-stan("dampen.stan",data=D)
mub<-extract(fitd,"mub")
mean(mub[[1]] > 0)
#0.89075
quantile(mub[[1]],probs=c(.5,0.05,0.95))
#        50%          5%         95% 
# 0.001074846 -0.000440626  0.002571707
mua<-extract(fitd,"mua")
quantile(mua[[1]],probs=c(.5,0.05,0.95))
#        50%          5%         95% 
#0.03593379 0.02029517 0.05023829 

pdf("dpOverTime.pdf",width=5,height=4)
par(mar=c(4.5,4.5,.5,.5))
plot(D$Yr,D$dp,pch=19,xlab="Year",ylab="Change",cex.lab=1.3)
abline(a=0.048688498,b=0.005265305)
dev.off()

oo<-extract(fitd,"beta")
bestd<-apply(oo[[1]],2,quantile,probs=c(0.5,0.05,0.95))
oo<-extract(fitd,"alpha")
aestd<-apply(oo[[1]],2,quantile,probs=c(0.5,0.05,0.95))

\
## cartoon illustrations ##
simp<-function(D=NA,p0=.1,phat=.5,t=20){
        p<-rep(p0,t)
        for(j in 2:t){
                del<-p[j-1]-phat
                dp<-D*del
                p[j]<-p[j-1]+dp
                if(p[j] < 0){p[j]<-0}
                if(p[j] > 1){p[j]<-1}
        }
        return(p)
}

simd<-matrix(NA,nrow=3,ncol=24)
simd[1,]<-simp(D=-2,p0=.35,phat=.5,t=24)
simd[2,]<-simp(D=-2.05,p0=.35,phat=.5,t=24)
simd[3,]<-simp(D=-1.95,p0=.35,phat=.5,t=24)

cl<-1.3;ca<-1.0;cm<-1.3
pdf("SF_dpOverTime.pdf",width=5,height=11)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,1))
xy<-2000:2023
cs<-c("black","firebrick","cadetblue")
plot(xy,simd[1,],type='l',ylim=c(0,1),col=cs[1],xlab="Year",ylab="Frequency",cex.lab=cl,cex.axis=ca)
lines(xy,simd[2,],col=cs[2])
lines(xy,simd[3,],col=cs[3])
title(main="(a) Possible NFDS dynamics",cex.main=cm)

plot(xy[-24],abs(simd[1,-1]-simd[1,-24]),pch=19,ylim=c(0,1),col=cs[1],xlab="Year",ylab="Change infrequency",cex.lab=cl,cex.axis=ca)
points(xy[-24],abs(simd[2,-1]-simd[2,-24]),pch=19,col=cs[2])
points(xy[-24],abs(simd[3,-1]-simd[3,-24]),pch=19,col=cs[3])

o<-lm(abs(simd[2,-1]-simd[2,-24]) ~ xy[-24])
abline(o$coefficients,col=cs[2])
o<-lm(abs(simd[3,-1]-simd[3,-24]) ~ xy[-24])
abline(o$coefficients,col=cs[3])
abline(h=0.3)
legend(2000,1,c("constant","increasing","dampening"),fill=cs,bty='n')
title(main="(b) Patterns of change with time",cex.main=cm)

cs<-alpha(brewer.pal(n=10,"Paired"),.7)
plot(D$Yr,D$dp,pch=19,col=cs[D$bl],xlab="Year",ylab="Absolute change in frequency",cex.lab=1.2,
     axes=FALSE)
axis(2)
axis(1,at=seq(0,20,5),seq(2000,2020,5))
for(j in 1:10){
        abline(a=aestd[1,j],b=bestd[1,j],col=cs[j])
}
legend(0,.22,ids[keep],pch=19,col=cs,bty='n',ncol=3)
title(main="(c) Empirical change with time",cex.main=cm)
box()
dev.off()
sum(bestd[2,] > 0)
#[1] 7
## 7 CIs exclude 0 of ten

save(list=ls(),file="nfds.rdat")

### main text figus

cl<-1.4;ca<-1.05;cm<-1.4
pdf("F2_StripeChangeNFDS.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
cs<-alpha(brewer.pal(n=10,"Paired"),.7)

pp<-seq(.05,.95,.05)
eqd<-rep(1,19)

pds<-rep(0.1,10)
for(j in 2:10){
	pds[j]<-pds[j-1] + pds[j-1] * (1-pds[j-1]) * .4 + rnorm(1,0,.02)
}
pnfds<-rep(0.4,10)
for(j in 2:10){
	pnfds[j]<-pnfds[j-1] + .6 + pnfds[j-1] * -2 +  rnorm(1,0,.02)
}

plot(1:10,pds,pch=19,type='b',xlab="Generation",ylim=c(0,1),ylab="Frequency",cex.lab=cl,cex.main=cm,col="cadetblue")
title(main="(a) Directional selection",cex.main=cm)


plot(1:10,pnfds,pch=19,type='b',xlab="Generation",ylim=c(0,1),ylab="Frequency",cex.lab=cl,cex.main=cm,col="firebrick")
title(main="(b) Negative frequency-dependent selection",cex.main=cm)

plot(pnfds[-10],pnfds[-1] - pnfds[-10],xlim=c(0,1),xlab="Frequency",ylab="Change in frequency",cex.lab=cl,cex.axis=ca,col="firebrick",pch=19)
points(pds[-10],pds[-1] - pds[-10],col="cadetblue",pch=19)
abline(h=0)
o<-lm(pnfds[-1] - pnfds[-10] ~ pnfds[-10])
abline(o$coefficients,col="firebrick")
o<-lm(pds[-1] - pds[-10] ~ pds[-10])
abline(o$coefficients,col="cadetblue")
legend(.25,.45,c("directional","negative frequency dependent"),fill=c("cadetblue","firebrick"),title="Form of selection")
title(main="(c) Frequency vs change",cex.main=cm)

plot(pest,dpest,pch=19,col=cs[blk],xlab="Stripe frequency",ylab="Change in frequency",cex.lab=cl,cex.axis=ca)
for(j in 1:10){
	abline(a=aest[1,j],b=best[1,j],col=cs[j])
}
legend(.55,.25,ids[keep],pch=19,col=cs,bty='n',ncol=2)
abline(h=0)
title(main="(d) Time series D estimates",cex.main=cm)

dev.off()

## supp figure
cl<-1.3;ca<-1;cm<-1.3
pdf("SF_StripeEquil.pdf",width=10,height=5)
par(mfrow=c(1,2))
par(mar=c(4.5,5.5,2.5,1.5))
cs<-alpha(brewer.pal(n=10,"Paired"),.7)
plot(jitter(pp,2.5),jitter(pp,2.5),xlim=c(0,1),ylim=c(0,1),col="firebrick",xlab="Mean frequency",ylab="Equilibrium frequency",pch=19,cex.lab=cl,cex.axis=ca)
points(jitter((pp+eqd)/2,2.5),eqd,pch=19,col="cadetblue")
abline(a=0,b=1,lwd=1.3,col="firebrick")
xx<-pp+eqd/2
o<-lm(eqd ~ xx)
abline(o$coefficients,lwd=1.3,col="cadetblue")
legend(0,0.95,c("directional","negative frequency","dependence"),fill=c("cadetblue","firebrick",NA),title="Form of selection")
title(main="(a) Alternative equilibrium patterns",cex.main=cm)

plot(eq,pbar,xlim=c(0,1),ylim=c(0,1),col=cs,xlab="Mean stripe frequency",ylab="Equilibrium stripe frequency",pch=19,cex.lab=cl,cex.axis=ca)
abline(a=0,b=1)
mtext("r = 0.99",3,line=-2,adj=.1,cex=1.3)
title(main="(b) Equilibrium predictions",cex.main=cm)

dev.off()

