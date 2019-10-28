
functions {
 // log of Gaussian density with mean zero
  real dlnorm(real x, real sigma) {
    return -log(2 * pi()) / 2 - log(sigma) - square(x) / (2 * square(sigma));
  }
  // spike and slab prior. w is an estimated parameter akin to probability of inclusion in the model
  real spikeSlab_lpdf(vector beta, real spike_scale, real slab_scale, vector w) {
    vector[num_elements(beta)] lprob;
    for (i in 1:num_elements(beta)) {
      lprob[i] = (1 - w[i]) * dlnorm(beta[i], spike_scale) + w[i] * dlnorm(beta[i], slab_scale);
     }
    return log_sum_exp(lprob);
    }
  // generalized double pareto shrinkage prior with parameters: scale (xi), shape (alpha)
  real gdp_lpdf(vector theta, real xi, real alpha) {
    vector[num_elements(theta)] lprob;
    for (i in 1:num_elements(theta)) {
      lprob[i] = -log(2 * xi) - (1 + alpha) * log( 1 + fabs(theta[i]) / (alpha*xi) );
	    }
	return log_sum_exp(lprob);
  }
  // regularized horseshoe
  //  vector horseshoe(vector z, real sigma, real aux1_global, real aux2_global,
  //		   real aux1_local, real aux2_local, real scale_global,
  //		   real slab_scale, real caux) {
  //  real<lower = 0> tau;
  //  vector<lower=0>[num_elements(z)] lambda;
  // vector<lower=0>[num_elements(z)] lambda_tilde;
  //  real <lower=0> c;
  //  tau = aux1_global * sqrt(aux2_global) * scale_global * sigma;
  //  lambda = aux1_local .* sqrt(aux2_local);
  //  c = slab_scale * sqrt(caux);
  //  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));
  //  return z .* lambda_tilde * tau; 
}

 









