## can drift explain the lack of dampening oscillations
library(scales)
library(RColorBrewer)

load(file="nfds.rdat")
## p0 = initial p, phat = equilibrium
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

## p0 = initial p, phat = equilibrium
simpDrift<-function(D=NA,p0=.1,phat=.5,t=20,Ne=100){
	p<-rep(p0,t)
	for(j in 2:t){
		del<-p[j-1]-phat
		dp<-D*del
		p[j]<-p[j-1]+dp
		p[j]<-rbinom(n=1,size=2*Ne,prob=p[j])/(2*Ne)
		if(p[j] < 0){p[j]<-0}
		if(p[j] > 1){p[j]<-1}
	}
	return(p)
}

## real data
ids<-c("FH","HV","L","M","MBOX","OG","OUT","PR","SC","VP")
cs<-alpha(brewer.pal(n=10,"Paired"),.7)
dat<-read.table("dpdf.txt",header=TRUE)

allYears<-1990:2023
Nyrs<-length(allYears)
Pmat<-matrix(NA,nrow=10,ncol=Nyrs)
for(j in 1:10){
        xx<-which(allYears %in% Years[[j]])
        Pmat[j,xx]<-P[[j]]
}

k2<-which(allYears >=2000)

pdf("OscillationsDrift.pdf",width=8,height=7)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,1))
for(i in 1:10){
	tt<-2023-dat$yr[i]+1
	p<-c(rep(NA,24-tt),simp(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],t=tt))
	pmat<-matrix(NA,nrow=20,ncol=24)
	for(j in 1:20){
		pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=110,t=tt))
	}
	plot(2000:2023,p,type='n',ylim=c(0,1),xlab="Year",ylab="Frequency")
	for(j in 1:20){
		lines(2000:2023,pmat[j,],col=alpha("darkgray",.7))
	}
	lines(2000:2023,p,col="black",lwd=1.5)
	points(allYears[k2],Pmat[i,k2],col=cs[i],type='b',pch=19,lwd=2)
	title(main=ids[i])
}
dev.off()

pdf("OscillationsDriftNe20.pdf",width=8,height=7)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,1))
for(i in 1:10){
	tt<-2023-dat$yr[i]+1
	p<-c(rep(NA,24-tt),simp(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],t=tt))
	pmat<-matrix(NA,nrow=20,ncol=24)
	for(j in 1:20){
		pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=20,t=tt))
	}
	plot(2000:2023,p,type='n',ylim=c(0,1),xlab="Year",ylab="Frequency")
	for(j in 1:20){
		lines(2000:2023,pmat[j,],col=alpha("darkgray",.7))
	}
	lines(2000:2023,p,col="black",lwd=1.5)
	points(allYears[k2],Pmat[i,k2],col=cs[i],type='b',pch=19,lwd=2)
	title(main=ids[i])
}
dev.off()

## stats
null<-matrix(NA,nrow=10,ncol=100)
for(i in 1:10){
	tt<-2023-dat$yr[i]+1
	pmat<-matrix(NA,nrow=100,ncol=24)
	for(j in 1:100){
		pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=110,t=tt))
	}
	dpmat<-abs(pmat[,-24]-pmat[,-1])
	null[i,]<-apply(dpmat,1,mean,na.rm=TRUE)
}

null20<-matrix(NA,nrow=10,ncol=100)
for(i in 1:10){
	tt<-2023-dat$yr[i]+1
	pmat<-matrix(NA,nrow=100,ncol=24)
	for(j in 1:100){
		pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=20,t=tt))
	}
	dpmat<-abs(pmat[,-24]-pmat[,-1])
	null20[i,]<-apply(dpmat,1,mean,na.rm=TRUE)
}

mns<-tapply(X=abs(dpest),INDEX=blk,mean)

ps<-rep(NA,10)
xf<-rep(NA,10)
for(i in 1:10){
	ps[i]<-mean(null[i,] >= mns[i])
	xf[i]<-mns[i]/mean(null[i,])
}
ps
#  [1] 0.03 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00

xf
#  [1] 1.32 1.62 1.60 1.69 2.65 1.70 2.33 2.48 1.61 2.25

ps20<-rep(NA,10)
xf20<-rep(NA,10)
for(i in 1:10){
	ps20[i]<-mean(null20[i,] >= mns[i])
	xf20[i]<-mns[i]/mean(null20[i,])
}
ps20## not much signal here with new Ds
xf20

## sensitivity
Ne<-rev(seq(20,110,10))
Nne<-length(Ne)
xfs<-matrix(NA,nrow=Nne,ncol=10)
pss<-matrix(NA,nrow=Nne,ncol=10)
for(k in 1:Nne){
	null<-matrix(NA,nrow=10,ncol=100)
	for(i in 1:10){
		tt<-2023-dat$yr[i]+1
		pmat<-matrix(NA,nrow=100,ncol=24)
		for(j in 1:100){
			pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=Ne[k],t=tt))
		}
		dpmat<-abs(pmat[,-24]-pmat[,-1])
		null[i,]<-apply(dpmat,1,mean,na.rm=TRUE)
		xfs[k,i]<-mns[i]/mean(null[i,])
		pss[k,i]<-mean(null[i,] >= mns[i])
	}
}	

save(list=ls(),file="drift.rdat")

sym<-(pss < 0.05) + 1

Ps<-pss
Ps[Ps == 0]<-0.01
chi<-apply(log(Ps),1,sum) * -2
pchisq(q=chi,df=20,lower.tail=FALSE)
round(pchisq(q=chi,df=20,lower.tail=FALSE),3)
# [1] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.001 0.129 0.925

## for Ne 110 for main text, with percision
#7.698246e-11


cl<-1.3;cm<-1.3;ca<-1.0
pdf("SF_Sensitivity.pdf",width=7,height=5.5)
par(mar=c(5,5,.5,.5))
plot(xfs[,1],ylim=c(0,2.75),type='n',axes=FALSE,xlab="Effective population size",ylab="X-fold mean change",cex.lab=cl)
polygon(c(-.5,8.5,8.5,-.5),c(-.5,-.5,3,3),border=NA,col=alpha("gray",.3))
abline(h=1,lty=3)
axis(1,at=1:10,Ne,cex.axis=ca)
axis(2,cex.axis=ca)
box()
for(j in 1:10){
lines(xfs[,j],pch=c(21,19)[sym[,j]],type='b',col=cs[j])
}
dev.off()


cl<-1.3;cm<-1.3;ca<-1.0
pdf("F4_OscillationsDrift.pdf",width=7,height=7)
par(mfrow=c(2,1))
par(mar=c(4.5,5,2.5,1))
plot(2000:2023,p,type='n',ylim=c(0,1),xlab="Year",ylab="Frequency",cex.lab=cl,cex.axis=ca)
for(i in c(2,7,10)){
	tt<-2023-dat$yr[i]+1
	p<-c(rep(NA,24-tt),simp(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],t=tt))
	pmat<-matrix(NA,nrow=20,ncol=24)
	for(j in 1:20){
		pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=110,t=tt))
	}
	for(j in 1:20){
		lines(2000:2023,pmat[j,],col=alpha("darkgray",.7))
	}
	lines(2000:2023,p,col="black",lwd=1.5)
	points(allYears[k2],Pmat[i,k2],col=cs[i],type='b',pch=19,lwd=2)
}
title(main="(a) Time series with drift",cex.main=cm)

boxplot(t(null),ylim=c(0,.085),pch=NA,xlab="Population",ylab="Mean change",col="white",names=ids,cex.lab=cl,cex.axis=ca)
for(k in 1:10){
	xx<-rnorm(100,k,.1)
	points(xx,null[k,],pch=19,col=alpha("gray",.1),cex=.5)
}
points(1:10,mns,col=cs,pch=19,cex=1.2)
title(main="(b) Change by drift",cex.main=cm)


dev.off()

cl<-1.3;cm<-1.3;ca<-1.0
pdf("F4_OscillationsDrift.pdf",width=11.3,height=7)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,1))
plot(2000:2023,p,type='n',ylim=c(0,1),xlab="Year",ylab="Frequency",cex.lab=cl,cex.axis=ca)
for(i in c(2,7,10)){
	tt<-2023-dat$yr[i]+1
	p<-c(rep(NA,24-tt),simp(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],t=tt))
	pmat<-matrix(NA,nrow=20,ncol=24)
	for(j in 1:20){
		pmat[j,]<-c(rep(NA,24-tt),simpDrift(D=dat$D[i],p0=dat$p0[i],phat=dat$eq[i],Ne=110,t=tt))
	}
	for(j in 1:20){
		lines(2000:2023,pmat[j,],col=alpha("darkgray",.7))
	}
	lines(2000:2023,p,col="black",lwd=1.5)
	points(allYears[k2],Pmat[i,k2],col=cs[i],type='b',pch=19,lwd=2)
}
title(main="(a) Effect of drift",cex.main=cm)

i<-2
oyr<-which((Years[[2]][-1] - Years[[2]][-25])==1)
oyr<-Years[[2]][oyr]
xyr<-2000:2022
dpd<-matrix(NA,nrow=23,ncol=2)
dpd[,1]<-xyr
dpd[which(xyr %in% oyr),2]<-dpest[blk==2][-1]

plot(dpd[,1],dpd[,2],col=cs[i],type='n',xlab="Year",ylab="Change in frequency",cex.lab=cl,cex.axis=ca)
points(dpd[,1],dpd[,2],col=cs[i],type='b',pch=19,lwd=2)
polygon(c(2020.5,2021.5,2021.5,2020.5),c(-.5,-.5,.5,.5),border=NA,col=alpha("red",.4))
title(main="(b) Effect of fire at HV",cex.main=cm)

boxplot(t(null),ylim=c(0,.2),pch=NA,xlab="Population",ylab="Mean change",col="white",names=ids,cex.lab=cl,cex.axis=ca)
for(k in 1:10){
	xx<-rnorm(100,k,.1)
	points(xx,null[k,],pch=19,col=alpha("gray",.1),cex=.5)
}
points(1:10,mns,col=cs,pch=19,cex=1.2)
title(main="(c) Change by drift (Ne = 110)",cex.main=cm)

boxplot(t(null20),ylim=c(0,.2),pch=NA,xlab="Population",ylab="Mean change",col="white",names=ids,cex.lab=cl,cex.axis=ca)
for(k in 1:10){
	xx<-rnorm(100,k,.1)
	points(xx,null20[k,],pch=19,col=alpha("gray",.1),cex=.5)
}
points(1:10,mns,col=cs,pch=19,cex=1.2)
title(main="(d) Change by drift (Ne = 110)",cex.main=cm)

dev.off()
