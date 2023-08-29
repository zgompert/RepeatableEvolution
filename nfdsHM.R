## NFDS analysis with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## read data
dat<-read.csv("TcristinaeGradient.csv",header=TRUE)
xx<-c(2:20) ## can't estimate relative fitness when fixed
yy<-10

D<-list(N=length(xx),S0=dat$Str_0[xx],G0=dat$Uns_0[xx],S1=dat$Str_1[xx],G1=dat$Uns_1[xx],p0=dat$Str_0[xx]/20,N_test=1,S0_test=dat$Str_0[yy],G0_test=dat$Uns_0[yy],S1_test=dat$Str_1[yy],G1_test=dat$Uns_1[yy])
fit<-stan("nfds_mod.stan",data=D,control=list(adapt_delta=0.9),iter=8000)
o<-summary(fit)

loo(fit)
#Computed from 16000 by 19 log-likelihood matrix

#         Estimate   SE
#elpd_loo    -70.9  6.7
#p_loo         5.6  1.1
#looic       141.8 13.3
#------
#Monte Carlo SE of elpd_loo is 0.1.

#All Pareto k estimates are good (k < 0.5)

sams<-vector("list",19)
for(i in 2:18){
	D<-list(N=length(xx[-i]),S0=dat$Str_0[xx][-i],G0=dat$Uns_0[xx][-i],S1=dat$Str_1[xx][-i],G1=dat$Uns_1[xx][-i],p0=dat$Str_0[xx][-i]/20,N_test=1,S0_test=dat$Str_0[xx][i],G0_test=dat$Uns_0[xx][i],S1_test=dat$Str_1[xx][i],G1_test=dat$Uns_1[xx][i])
	fit_cv<-stan("nfds_mod.stan",data=D,control=list(adapt_delta=0.9),iter=8000)
	sams[[i]]<-extract(fit_cv,pars=c("s1_err","g1_err","s1_test_est","g1_test_est"))
}

s_est<-matrix(NA,nrow=16000,ncol=17)
g_est<-matrix(NA,nrow=16000,ncol=17)
s_err<-matrix(NA,nrow=16000,ncol=17)
g_err<-matrix(NA,nrow=16000,ncol=17)
for(i in 1:17){
	s_err[,i]<-sams[[i+1]][[1]]
	g_err[,i]<-sams[[i+1]][[2]]
	s_est[,i]<-sams[[i+1]][[3]]
	g_est[,i]<-sams[[i+1]][[4]]
}


save(list=ls(),file="nfdsfit.rdat")

## sigmoid

dat<-read.csv("TcristinaeGradient.csv",header=TRUE)
xx<-c(2:20) ## can't estimate relative fitness when fixed
yy<-10

D<-list(N=length(xx),S0=dat$Str_0[xx],G0=dat$Uns_0[xx],S1=dat$Str_1[xx],G1=dat$Uns_1[xx],p0=dat$Str_0[xx]/20,N_test=1,S0_test=dat$Str_0[yy],G0_test=dat$Uns_0[yy],S1_test=dat$Str_1[yy],G1_test=dat$Uns_1[yy])
fit_sig<-stan("nfds_mod_sig.stan",data=D,control=list(adapt_delta=0.95),iter=20000)
o_sig<-summary(fit_sig)

loo(fit_sig)

#Computed from 40000 by 19 log-likelihood matrix

#         Estimate   SE
#elpd_loo    -71.6  7.0
#p_loo         6.4  1.3
#looic       143.2 14.1
#------
#Monte Carlo SE of elpd_loo is 0.2.

#All Pareto k estimates are good (k < 0.5).


sams_sig<-vector("list",19)
for(i in 2:18){
	D<-list(N=length(xx[-i]),S0=dat$Str_0[xx][-i],G0=dat$Uns_0[xx][-i],S1=dat$Str_1[xx][-i],G1=dat$Uns_1[xx][-i],p0=dat$Str_0[xx][-i]/20,N_test=1,S0_test=dat$Str_0[xx][i],G0_test=dat$Uns_0[xx][i],S1_test=dat$Str_1[xx][i],G1_test=dat$Uns_1[xx][i])
	fit_sig_cv<-stan("nfds_mod_sig.stan",data=D,control=list(adapt_delta=0.95),iter=20000)
	sams_sig[[i]]<-extract(fit_sig_cv,pars=c("s1_err","g1_err","s1_test_est","g1_test_est"))
}

s_sig_est<-matrix(NA,nrow=40000,ncol=17)
g_sig_est<-matrix(NA,nrow=40000,ncol=17)
s_sig_err<-matrix(NA,nrow=40000,ncol=17)
g_sig_err<-matrix(NA,nrow=40000,ncol=17)
for(i in 1:17){
	s_sig_err[,i]<-sams_sig[[i+1]][[1]]
	g_sig_err[,i]<-sams_sig[[i+1]][[2]]
	s_sig_est[,i]<-sams_sig[[i+1]][[3]]
	g_sig_est[,i]<-sams_sig[[i+1]][[4]]
}

### plots
hm_rw<-extract(fit,pars="rws")
hm_wbar<-extract(fit,pars="wbar")
hm_ws<-extract(fit,pars="ws")
hm_wg<-extract(fit,pars="wg")

est_ws<-apply(hm_ws[[1]],2,quantile,probs=c(.5,0.05,.95))
est_wg<-apply(hm_wg[[1]],2,quantile,probs=c(.5,0.05,.95))
est_rw<-apply(hm_rw[[1]],2,quantile,probs=c(.5,0.05,.95))
est_wbar<-apply(hm_wbar[[1]],2,quantile,probs=c(.5,0.05,.95))

xp<-seq(0.05,0.95,0.05)

library(scales)
plotCIs<-function(X=NA,Y=NA,lb=NA,ub=NA,lims=c(0,1),xl="Stripe frequency",yl="Value",cc="black",cl=1.3,ca=1.0){
	cs1<-alpha(cc,.5)
	cs2<-cc
	plot(X,Y,ylim=lims,type='n',xlab=xl,ylab=yl,cex.lab=cl,cex.axis=ca)
	polygon(c(X,rev(X)),c(lb,rev(ub)),col=cs1,border=NA)
	lines(X,Y,col=cs2)
}

pdf("nfdsBayes.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(5,5,1.5,1.5))
plotCIs(X=xp,Y=est_ws[1,],lb=est_ws[2,],ub=est_ws[3,],yl="Fitness (stripe)",cc="cadetblue")
plotCIs(X=xp,Y=est_wg[1,],lb=est_wg[2,],ub=est_wg[3,],yl="Fitness (green)",cc="forestgreen")
plotCIs(X=xp,Y=est_rw[1,],lb=est_rw[2,],ub=est_rw[3,],lims=c(0.5,2.1),yl="Relative fitness (stripe)",cc="violet")
abline(h=1,lty=2,col="red")
plotCIs(X=xp,Y=est_wbar[1,],lb=est_wbar[2,],ub=est_wbar[3,],yl="Mean fitness",cc="black")
dev.off()


pdf("nfdsBayesData.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(5,5,1.5,1.5))
plotCIs(X=xp,Y=est_ws[1,],lb=est_ws[2,],ub=est_ws[3,],yl="Fitness (stripe)",cc="cadetblue")
points(xp,D$S1/D$S0,pch=19)
plotCIs(X=xp,Y=est_wg[1,],lb=est_wg[2,],ub=est_wg[3,],yl="Fitness (green)",cc="forestgreen")
points(xp,D$G1/D$G0,pch=19)
plotCIs(X=xp,Y=est_rw[1,],lb=est_rw[2,],ub=est_rw[3,],lims=c(0.5,2.1),yl="Relative fitness (stripe)",cc="violet")
wbar<-(D$S1/D$S0+D$G1/D$G0)/2
points(xp,(D$S1/D$S0)/wbar,pch=19)
abline(h=1,lty=2,col="red")
plotCIs(X=xp,Y=est_wbar[1,],lb=est_wbar[2,],ub=est_wbar[3,],yl="Mean fitness",cc="black")
points(xp,wbar,pch=19)
dev.off()

## figure for repeatability

## read data
dat<-read.csv("TcristinaeGradient.csv",header=TRUE)
xx<-c(2:20) ## can't estimate relative fitness when fixed
yy<-10

D<-list(N=length(xx),S0=dat$Str_0[xx],G0=dat$Uns_0[xx],S1=dat$Str_1[xx],G1=dat$Uns_1[xx],p0=dat$Str_0[xx]/20,N_test=1,S0_test=dat$Str_0[yy],G0_test=dat$Uns_0[yy],S1_test=dat$Str_1[yy],G1_test=dat$Uns_1[yy])


cl<-1.3;ca<-1.0;cm<-1.3
pdf("F3_nfdsExperiment.pdf",width=8,height=8)
par(mfrow=c(2,2))

par(mar=c(4.5,5.5,2.5,1.5))

plot(xp,1-(.5*xp),ylim=c(0,1),type='l',ylab="Fitness",cex.lab=cl,cex.axis=ca,lwd=1.5,xlab="Stripe frequency")
ys<-rep(.8,length(xp))
ys[xp > .7]<-.3
lines(xp,ys,ylim=c(0,1),type='l',ylab="Fitness",cex.lab=cl,cex.axis=ca,lwd=1.5)
title(main="(a) Alternative fitness functions",cex.main=cm)

plotCIs(X=xp,Y=est_ws[1,],lb=est_ws[2,],ub=est_ws[3,],yl="Fitness (stripe)",cc="cadetblue")
points(xp,D$S1/D$S0,pch=19)
title(main="(b) Stripe fitness",cex.main=cm)

plotCIs(X=xp,Y=est_wg[1,],lb=est_wg[2,],ub=est_wg[3,],yl="Fitness (green)",cc="forestgreen")
points(xp,D$G1/D$G0,pch=19)
title(main="(c) Green fitness",cex.main=cm)

plotCIs(X=xp,Y=est_rw[1,],lb=est_rw[2,],ub=est_rw[3,],lims=c(0.5,2.1),yl="Relative fitness (stripe)",cc="firebrick")
wbar<-(D$S1/D$S0+D$G1/D$G0)/2
points(xp,(D$S1/D$S0)/wbar,pch=19)
abline(h=1,lty=2,col="black")
title(main="(d) Relative fitness",cex.main=cm)
dev.off()

## for sigmoid model
### plots
sig_rw<-extract(fit_sig,pars="rws")

est_sig_rw<-apply(sig_rw[[1]],2,quantile,probs=c(.5,0.05,.95))

xp<-seq(0.05,0.95,0.05)

library(scales)
plotCIs<-function(X=NA,Y=NA,lb=NA,ub=NA,lims=c(0,1),xl="Stripe frequency",yl="Value",cc="black",cl=1.3,ca=1.0){
        cs1<-alpha(cc,.5)
        cs2<-cc
        plot(X,Y,ylim=lims,type='n',xlab=xl,ylab=yl,cex.lab=cl,cex.axis=ca)
        polygon(c(X,rev(X)),c(lb,rev(ub)),col=cs1,border=NA)
        lines(X,Y,col=cs2)
}


cl<-1.3;ca<-1.0;cm<-1.3
pdf("SF_sigNfdsExperiment.pdf",width=6,height=6)

par(mar=c(4.5,5.5,2.5,1.5))
plotCIs(X=xp,Y=est_sig_rw[1,],lb=est_sig_rw[2,],ub=est_sig_rw[3,],lims=c(0.5,2.1),yl="Relative fitness (stripe)",cc="firebrick")
wbar<-(D$S1/D$S0+D$G1/D$G0)/2
points(xp,(D$S1/D$S0)/wbar,pch=19)
abline(h=1,lty=2,col="black")
dev.off()


## predicted equilibrium from experiment with same model as natural populations

D<-list(nn=length(xx),Y0=dat$Str_0[xx],N0=dat$Uns_0[xx]+dat$Str_0[xx],Y1=dat$Str_1[xx],N1=dat$Uns_1[xx]+dat$Str_1[xx])

fit_eq<-stan("simpleD.stan",data=D,control=list(adapt_delta=0.95),warmup=3000,iter=20000)
o_eq<-summary(fit_eq)
## D
dd<-extract(fit_eq,"beta")
quantile(dd[[1]],probs=c(.5,0.025,.975))
#       50%       2.5%      97.5%
#-0.7595888 -1.2140989 -0.3305669

## equli
-0.57/-0.76
[1] 0.75



## old panel A
#par(mar=c(4.5,3,2.5,3))

#plot(c(0,1),c(0,1),type='n',xlab="Percent stripe",ylab="",axes=FALSE,cex.lab=cl)
#lines(c(0,1),c(.2,.2),lwd=2)
#x<-(1:19)/20
#segments(x,rep(.1,19),x,rep(.3,19))
#text(x,rep(.39,19),paste(x*100,"%",sep=""),srt=90)
#text(x[1],.8,"1 stripe")
#text(x[1],.73,"19 green")
#text(x[10],.8,"10 stripe")
#text(x[10],.73,"10 green")
#text(x[19],.8,"19 stripe")
#text(x[19],.73,"1 green")
#lines(c(x[1],x[1]),c(.47,.7),lty=3,col="gray20")
#lines(c(x[10],x[10]),c(.47,.7),lty=3,col="gray20")
#lines(c(x[19],x[19]),c(.47,.7),lty=3,col="gray20")

#title(main="(a) Experimental design",cex.main=cm)
