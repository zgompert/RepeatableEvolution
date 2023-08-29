data{
	int<lower=0> N; /* # bushes */
	int S0[N]; /* vector of initial stripe count */
	int G0[N]; /* vector of initial green count */
	int S1[N]; /* vector of final stripe count */
	int G1[N]; /* vector of final green count */
	vector<lower=0.0,upper=1.0>[N] p0; /* initial stripe frequency */
	/* for cross validation */
	int<lower=0> N_test;
	int S0_test;
	int G0_test;
	int S1_test;
	int G1_test;
}


parameters{
	real mu;
	real beta;
	real<lower=0.0> L; /* max of sigmoid */
	real k; /* steepness of curve */
	real<lower=0.0,upper=1.0> x0; /* x-value for midpoint */
}

transformed parameters {
	vector<lower=0.0,upper=1.0>[N] wbar; /* mean fitness of morphs */
	vector<lower=0.0>[N] rws; /* relative fitness of stripe */
	vector<lower=0.0,upper=1.0>[N] ws; /* absolute fitness of stripe */
	vector<lower=0.0,upper=1.0>[N] wg; /* absolute fitness of green */
	/* linear models */
	wbar = inv_logit(mu + beta * p0);
	for(i in 1:N){
		rws[i] = L/(1 + exp(-1 * k * (p0[i]-x0)));
	}
	for(i in 1:N){
		ws[i] = rws[i] * wbar[i];
		wg[i] = 2*wbar[i] - rws[i] * wbar[i];
	} 
}

model{

	/* likelihood */
	for(i in 1:N){
		/* increment binomial based on survival */
		target += binomial_lpmf(S1[i] | S0[i], ws[i]);
		target += binomial_lpmf(G1[i] | G0[i], wg[i]);
	}
	/* priors for terms in model */
	target += normal_lpdf(mu | 0, 10);
	target += normal_lpdf(beta | 0, 10);
	target += beta_lpdf(x0 | 1, 1);
	target += normal_lpdf(k | 0, 10);
	target += normal_lpdf(L | 0, 10);
}

generated quantities {
	/* linear models */
        vector[N] log_lik;
	real<lower=0.0,upper=1.0> p0_test = S0_test/20.0;
	real<lower=0.0,upper=1.0> wbar_test = inv_logit(mu + beta * p0_test);
	real<lower=0.0> rws_test = L/(1 + exp(-1 * k * (p0_test-x0)));
	real<lower=0.0,upper=1.0> ws_test = rws_test * wbar_test;
	real<lower=0.0,upper=1.0> wg_test = 2*wbar_test - rws_test * wbar_test;
	int s1_test_est = binomial_rng(S0_test, ws_test);
	int s1_err = S1_test - s1_test_est;
	int g1_test_est = binomial_rng(G0_test, wg_test);
	int g1_err = G1_test - g1_test_est;

        for (i in 1:N){
                 log_lik[i] = binomial_lpmf(S1[i] | S0[i], ws[i]) + binomial_lpmf(G1[i] | G0[i], wg[i]);
        } 
}

