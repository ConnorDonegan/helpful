#include sparsity-priors.stan

data {
  int n;
  int n_ev;
  int k;
  vector[n] y;
  matrix[n, k] x;
  matrix[n, n_ev] ev;
 real gdp_scale;
 real gdp_shape;
}

transformed data {
    // use the QR decomposition on the matrix of covariates
  matrix[n, k] Q_ast;
  matrix[k, k] R_ast;
  matrix[k, k] R_inverse;
  Q_ast = qr_Q(x)[, 1:k] * sqrt(n - 1);
  R_ast = qr_R(x)[1:k, ] / sqrt(n - 1);
  R_inverse = inverse(R_ast);
}

parameters {
  vector[n_ev] beta_ev;
  vector[k] beta_tilde;
  real alpha;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[k] beta;
  vector[n] fitted;
  fitted = alpha + ev * beta_ev + Q_ast * beta_tilde;
  beta = R_inverse * beta_tilde;
}

model {
  beta_ev ~ gdp(gdp_scale, gdp_shape);
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
  sigma ~ cauchy(0, 2);
  target +=  normal_lpdf(y | fitted, sigma);
}

generated quantities{
  vector[n] sf;
  vector[n] residual;
  vector[n] log_lik;
  sf = ev * beta_ev;
  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma);
    residual[i] = fitted[i] - y[i];
  }
}

