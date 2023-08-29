data{
	int<lower=1> N; /* sample size */
	int<lower=0> y[N]; /* Stripe recapture */
	int<lower=0> n[N]; /* Total recapature */
}

parameters{
	real beta; /* intercept term */
}

transformed parameters{
	real<lower=0, upper=1> p; /* probability of stripe recap */
	p = inv_logit(beta);
}

model{
	/* prior on beta */
	target += normal_lpdf(beta | 0, 10);

	/* linear model and likelihood */
	for(i in 1:N){
		target += binomial_lpmf(y[i] | n[i], p);
	}
}

