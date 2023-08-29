data{
	int<lower=1> N; /* # of observations */
	int<lower=1> J; /* # of populations */
	real dp[N]; /* abs change in stripe */
	int<lower=0> Yr[N]; /* year */
	int<lower=1> bl[N]; /* block */
}

parameters{
	real alpha[J]; /* intercepts */
	real beta[J]; /* intercepts */
	real mua; /* hierarchical mean int */
	real<lower=0> sigmaa; /* hierarchical sd int */
	real mub; /* hierarchical mean slope */
	real<lower=0> sigmab; /* hierarchical sd slope */
	real<lower=0> sig; /* residual SD */
}

model{
	target += normal_lpdf(mua | 0, 1);
	target += gamma_lpdf(sigmaa | 2, 0.1);	
	target += normal_lpdf(mub | 0, 1);
	target += gamma_lpdf(sigmab | 2, 0.1);	
	target += gamma_lpdf(sig | 2, 0.1);	
	for(j in 1:J){
		target += normal_lpdf(alpha[j] | mua, sigmaa);
		target += normal_lpdf(beta[j] | mub, sigmab);
	}	
	for(n in 1:N){
		target += normal_lpdf(dp[n] | alpha[bl[n]] + Yr[n] * beta[bl[n]], sig);
	}
}

