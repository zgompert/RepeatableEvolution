## analyze the melanic stripe selection experiment
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## read data
dat<-read.table("SelectionMelanicStripe.txt",header=TRUE)

df<-cbind(dat$stripe,dat$melanic)

## plot
pdf("StripeMelanRecap.pdf",width=7,height=8)
par(mfrow=c(1,2))
par(mar=c(4,4.5,2.5,1))
cl<-1.35;ca<-1.1;cm<-1.3
cs<-c("forestgreen","brown")
barplot(c(80,80),col=cs,ylab="Number released",names.arg=c("Stripe","Melanic"),cex.lab=cl,cex.axis=ca)
title(main="(A) Release",cex.main=cm)
for(i in 1:7){
	barplot(c(80-10*i,80-10*i),col=cs,add=TRUE,axes=FALSE)
}
barplot(apply(df,2,sum),ylim=c(0,80),col=cs,ylab="Number recaptured",names.arg=c("Stripe","Melanic"),cex.lab=cl,cex.axis=ca)
title(main="(B) Recapture",cex.main=cm)
for(i in 1:6){
	barplot(apply(df[1:(8-i),],2,sum),col=cs,add=TRUE,axes=FALSE)
}
barplot(df[1,],col=cs,add=TRUE,axes=FALSE)
dev.off()

pdf("StripeMelanRecap.pdf",width=7,height=4.5)
cl<-1.35;ca<-1.1;cm<-1.3
cs<-c("forestgreen","brown")
par(mar=c(5,5,1,1))
barplot(t(df),beside=TRUE,col=cs,names.arg=1:8,ylab="Number recaptured",cex.lab=cl,cex.axis=ca,ylim=c(0,10))
dev.off()

## simple binomial GLM


o<-glm(df ~ 1,family=binomial)
## int > 1, thus p > .5
#Coefficients:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept)   1.3291     0.3749   3.546 0.000392 ***


D<-list(y=dat$stripe, n=dat$stripe+dat$melanic,N=8)
o_glm<-stan("glm_simple.stan",data=D)
#Inference for Stan model: glm_simple.
#4 chains, each with iter=2000; warmup=1000; thin=1; 
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#       mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#beta   1.36    0.01 0.39   0.65   1.09   1.35   1.62   2.18  1393    1
#p      0.79    0.00 0.06   0.66   0.75   0.79   0.83   0.90  1444    1
#lp__ -13.12    0.02 0.75 -15.18 -13.27 -12.84 -12.65 -12.59  1796    1

#Samples were drawn using NUTS(diag_e) at Thu May 18 10:47:07 2023.
#For each parameter, n_eff is a crude measure of effective sample size,
#and Rhat is the potential scale reduction factor on split chains (at 
#convergence, Rhat=1).

o_hglm<-stan("glm_hm.stan",data=D,warmup=3000,iter=6000)
#Inference for Stan model: glm_hm.
#4 chains, each with iter=6000; warmup=3000; thin=1; 
#post-warmup draws per chain=3000, total post-warmup draws=12000.

#           mean se_mean   sd   2.5%    25%    50%    75% 97.5% n_eff Rhat
#beta       1.54    0.01 0.55   0.58   1.16   1.49   1.86  2.79  3502 1.00
#alpha[1]  -0.17    0.01 0.63  -1.66  -0.47  -0.09   0.18  0.99  5117 1.00
#alpha[2]  -0.62    0.02 0.78  -2.54  -1.03  -0.44  -0.06  0.48  2042 1.00
#alpha[3]   0.31    0.01 0.81  -1.04  -0.12   0.16   0.63  2.34  5233 1.00
#alpha[4]  -0.16    0.01 0.63  -1.62  -0.46  -0.09   0.19  1.05  5965 1.00
#alpha[5]   0.37    0.01 0.82  -0.92  -0.09   0.20   0.71  2.44  4171 1.00
#alpha[6]   0.00    0.01 0.89  -1.89  -0.36   0.00   0.37  1.92  8722 1.00
#alpha[7]   0.51    0.02 0.83  -0.64  -0.01   0.31   0.87  2.62  2658 1.00
#alpha[8]  -0.23    0.01 0.66  -1.80  -0.55  -0.14   0.14  0.94  5130 1.00
#sig        0.74    0.02 0.49   0.12   0.36   0.64   1.02  1.91  1025 1.01
#repp[1]    0.78    0.00 0.10   0.55   0.73   0.79   0.85  0.94 11292 1.00
#repp[2]    0.70    0.00 0.14   0.36   0.63   0.73   0.80  0.90  3306 1.00
#repp[3]    0.83    0.00 0.10   0.61   0.77   0.84   0.91  0.99  4780 1.00
#repp[4]    0.78    0.00 0.10   0.55   0.73   0.80   0.85  0.94 10599 1.00
#repp[5]    0.84    0.00 0.10   0.62   0.78   0.85   0.91  0.99  3393 1.00
#repp[6]    0.79    0.00 0.14   0.40   0.72   0.81   0.89  0.98  7253 1.00
#repp[7]    0.86    0.00 0.09   0.67   0.80   0.87   0.92  0.99  1964 1.00
#repp[8]    0.77    0.00 0.11   0.51   0.71   0.78   0.85  0.93 11604 1.00
#p          0.81    0.00 0.08   0.64   0.76   0.82   0.86  0.94  4046 1.00
#lp__     -20.88    0.21 5.24 -30.27 -24.50 -21.33 -17.52 -9.85   637 1.01

#Samples were drawn using NUTS(diag_e) at Thu May 18 11:09:45 2023.
#For each parameter, n_eff is a crude measure of effective sample size,
#and Rhat is the potential scale reduction factor on split chains (at 
#convergence, Rhat=1).


## check sum to zero of random effect
hmc<-extract(o_hglm)
summary(apply(hmc$alpha,1,mean))
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-2.0739248 -0.1337783  0.0003447  0.0017306  0.1346667  2.2440652 

save(list=ls(),file="stripeMelanic.rdat")
