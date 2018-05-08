functions{
  real pn_sphere_lpdf(row_vector theta,vector mu,matrix sigma){
    vector[3] u;
    real A;
    real B;
    real C;
    real D;
    real tmp;
    real p;
    u[1] = cos(theta[1]) * sin(theta[2]);
    u[2] = sin(theta[1]) * sin(theta[2]);
    u[3] = cos(theta[2]);
    A = u' * inverse(sigma) * u;
    B = u' * inverse(sigma) * mu;
    C = (-0.5) * (mu' * inverse(sigma) * mu);
    D = B/sqrt(A);
    tmp = normal_cdf(D,0,1)/(exp(-0.5*(D^2))/sqrt(2*pi()));
    p = - 1.5 * log(A) - 0.5 * log(determinant(sigma)) + C 
        + log( (1 + (D * tmp))*D + tmp);
    return p;
  }
}

data{
  int N; //sample size
  int M; //mixture num
  vector<lower=0,upper=2*pi()>[N] theta1; //data
  vector<lower=0,upper=pi()>[N] theta2; //data
}

transformed data{
  matrix[N,2] theta;
  theta[1:N,1] = theta1;
  theta[1:N,2] = theta2;
}
parameters{
  vector[3] mu[M];
  // unit_vector[3] mu[M];
  vector[2] gamma[M];
  cov_matrix[2] star[M];
  simplex[M] alpha;
}

transformed parameters{
  cov_matrix[3] sigma[M];
  for(m in 1:M){
    sigma[m,1:2,1:2] = star[m] + (gamma[m] * gamma[m]');
    sigma[m,3,1:2] = gamma[m]';
    sigma[m,1:2,3] = gamma[m];
    sigma[m,3,3] = 1.0;
  }
}

model{
  mu ~ multi_normal(rep_vector(0,3),diag_matrix(rep_vector(10^5,3)));
  gamma ~ multi_normal(rep_vector(0,2),diag_matrix(rep_vector(10^5,2)));
  alpha ~ dirichlet(rep_vector(2.0,M));
  for(m in 1:M){
    star[m] ~ inv_wishart(4,diag_matrix(rep_vector(1,2)));
  }
  for(n in 1:N){
    vector[M] ps;
    for(m in 1:M){
      ps[m] = log(alpha[m]) + pn_sphere_lpdf(theta[n]|mu[m],sigma[m]);
    }
    target += log_sum_exp(ps);
  }

}
generated quantities{
  vector[N] log_likelihood;
  for(n in 1:N){
    vector[M] ps;
    for(m in 1:M){
      ps[m] = log(alpha[m]) + pn_sphere_lpdf(theta[n]|mu[m],sigma[m]);
    }
    log_likelihood[n] =  log_sum_exp(ps);
  }
}

