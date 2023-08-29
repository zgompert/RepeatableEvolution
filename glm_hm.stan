data{
	int<lower=1> N; /* sample size */
	int<lower=0> y[N]; /* Stripe recapture */
	int<lower=0> n[N]; /* Total recapature */
}

parameters{
	real beta; /* intercept term */
	real alpha[N]; /* random effects */
	real <lower=0> sig; /* SD for random effects */
}

transformed parameters{
	real<lower=0, upper=1> repp[N]; /* probability of stripe recap, rep */
	real<lower=0, upper=1> p; /* probability of stripe recap */
	p = inv_logit(beta);
	for(i in 1:N){
		repp[i] = inv_logit(beta + alpha[i]);
	}
}

model{
	/* prior on beta */
	target += normal_lpdf(beta | 0, 10);
	target += normal_lpdf(sig | 0, 1);
	for(i in 1:N){
		target += normal_lpdf(alpha[i] | 0, sig);
	}

	/* linear model and likelihood */
	for(i in 1:N){
		target += binomial_lpmf(y[i] | n[i], repp[i]);
	}
}

