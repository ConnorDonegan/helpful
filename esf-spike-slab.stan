
#include sparsity-priors.stan

data {
  //number of observations
  int n;
  //number of eigenvectors
  int n_ev;
  //number of covariates
  int k;
  //outcome variable
  vector[n] y;
  //matrix of covariates
  matrix[n, k] x;
  //matrix of eigenvectors for spatial filtering
  matrix[n, n_ev] ev;
  // scale parameters for the two normal distributions, spike at 0 and a wider slab
  real spike_scale;
  real slab_scale;
}

transformed data {
    // use the QR decomposition on the matrix of covariates
  matrix[n, k] Q_ast;
  matrix[k, k] R_ast;
  matrix[k, k] R_inverse;
  //hyperparameters for prior on weights w
  real shape1 = 1;
  real shape2 = 1;
  Q_ast = qr_Q(x)[, 1:k] * sqrt(n - 1);
  R_ast = qr_R(x)[1:k, ] / sqrt(n - 1);
  R_inverse = inverse(R_ast);
}

parameters {
  //weights for spike and slab style mixture of normals
  vector<lower=0, upper=1>[n_ev] w;
  //parameters for the eigenvectors
  vector[n_ev] beta_ev;
  //parameters for the covariates (with qr decomposition)
  vector[k] beta_tilde;
  //intercept and standard deviation of the model
  real alpha;
  real<lower = 0> sigma;
  // real<lower=0> spike_scale;
  // real<lower=0> slab_scale;
}

transformed parameters {
  vector[k] beta;
  vector[n] fitted;
  fitted = alpha + ev * beta_ev + Q_ast * beta_tilde;
  beta = R_inverse * beta_tilde;
}

model {
  w ~ beta(shape1, shape2);
  beta_ev ~ spikeSlab(spike_scale, slab_scale, w);
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
  sigma ~ cauchy(0, 2);
  target +=  normal_lpdf(y | fitted, sigma);
}

generated quantities{
  vector[n] sf;
  vector[n] residual;
  vector[n] log_lik;
  //spatial filter
  sf = ev * beta_ev;
  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma);
    residual[i] = fitted[i] - y[i];
  }
}

