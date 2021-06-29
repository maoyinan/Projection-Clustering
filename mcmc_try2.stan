data {
  int<lower=0> N; // total number of observations
  vector[N] y; // response
  int<lower=0> nBasis;
  row_vector[nBasis] Z[N];
  int<lower=0> nSub; // subjects
  int<lower=1,upper=nSub> sub[N];
}
parameters {
  simplex[nBasis] pai; // covariance trace proportion
  real<lower=0> tau;      // prior scale
  cholesky_factor_corr[nBasis] LOmega;  // correlation LKJ distribution
  vector[nBasis] b[nSub]; // standard normal for coefficients
  real<lower=0> sigmaepsilon;
  real beta0;
}
transformed parameters {
  vector[N] yhat;
  vector[nBasis] beta[nSub]; // random effects of part ZA
  vector<lower=0>[nBasis] sigma= nBasis * tau^2 * pai;
  matrix[nBasis,nBasis] SigmaBeta= diag_pre_multiply(sigma, LOmega);
  matrix[nBasis,nBasis] VarBeta= SigmaBeta*(SigmaBeta');

  for(i in 1:nSub)
    beta[i] = SigmaBeta * b[i];
  for (i in 1:N){
    yhat[i] = beta0 + Z[i]  * beta[sub[i]];
    // print("i: ", i, " yhat[i] ", yhat[i]);
  }
}
model {
  vector[nBasis] alpha= rep_vector(1,nBasis); // concentration
  pai ~ dirichlet(alpha);
  tau ~ gamma(1, 1); // shape, scale
  LOmega ~ lkj_corr_cholesky(1); // regularization

  for( i in 1:nSub){
    for(j in 1:nBasis)
      b[i,j] ~ normal(0,1);
  }

  sigmaepsilon ~ exponential(1); // prior of sigma of error term
  y ~ normal(yhat, sigmaepsilon);

}


