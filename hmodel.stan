data{
	int<lower=1> nn; /* # of observations */
	int<lower=1> J; /* # of populations */
	int<lower=0> Y0[nn]; /* stripe count time 0 */
	int<lower=0> Y1[nn]; /* stripe count time 1 */
	int<lower=0> N0[nn]; /* total count time 0 */
	int<lower=0> N1[nn]; /* total count time 1 */
	int<lower=1> blk[nn]; /* block */
}

parameters{
	real alpha[J]; /* intercepts */
	real beta[J]; /* intercepts */
	real mua; /* hierarchical mean int */
	real<lower=0> sigmaa; /* hierarchical sd int */
	real mub; /* hierarchical mean slope */
	real<lower=0> sigmab; /* hierarchical sd slope */
	real<lower=0, upper=1> p0[nn]; /* time 0 stripe freq */
	real<lower=0, upper=1> p1[nn]; /* time 1 stripe freq */
	real<lower=0> sig; /* residual SD */
}

transformed parameters{
	real<lower=-2, upper=2> dp[nn]; /* change in stripe frequency */
	for(n in 1:nn){
		dp[n] = p1[n] - p0[n];
	}
}

model{
	target += normal_lpdf(mua | 0, 20);
	target += gamma_lpdf(sigmaa | 2, 0.1);	
	target += normal_lpdf(mub | 0, 20);
	target += gamma_lpdf(sigmab | 2, 0.1);	
	target += gamma_lpdf(sig | 2, 0.1);	
	for(j in 1:J){
		target += normal_lpdf(alpha[j] | mua, sigmaa);
		target += normal_lpdf(beta[j] | mub, sigmab);
	}	
	for(n in 1:nn){
		target += binomial_lpmf(Y0[n] | N0[n], p0[n]);
		target += binomial_lpmf(Y1[n] | N1[n], p1[n]);
		target += normal_lpdf(dp[n] | alpha[blk[n]] + p0[n] * beta[blk[n]], sig);
	}
}

