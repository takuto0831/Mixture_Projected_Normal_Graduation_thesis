functions{
  real pn_circle_lpdf(real theta,vector mu,matrix sigma){
    vector[2] u;
    real A; real B; real C; real D; real p;
    u[1] = cos(theta); u[2] = sin(theta);
    A = quad_form(inverse_spd(sigma), u); B = u' * inverse(sigma) * mu;
    C = (-0.5) * quad_form(inverse_spd(sigma), mu); D = B/sqrt(A);

    p = -log(A)- 0.5*log(determinant(sigma)) + C
    + log(1+(D * normal_cdf(D,0,1)/exp(normal_lpdf(D|0,1))));    
    return p;
  }
}

data{
  int N; //sample size
  int M; //mixture num
  real<lower=0,upper=2*pi()> theta[N]; //data
}

parameters{
  vector[2] mu[M];
  //ordered[M] mu_tmp;
  //unit_vector[2] mu[M]; 
  real<lower=0.001,upper=100> tau[M];
  real rho[M];
  simplex[M] alpha;
}

transformed parameters{
  cov_matrix[2] sigma[M];
  //vector[2] mu[M];
  //real tmp;
  for(m in 1:M){
    sigma[m,1,1] = tau[m]; sigma[m,1,2] = sqrt(tau[m])*rho[m];
    sigma[m,2,1] = sqrt(tau[m])*rho[m]; sigma[m,2,2] = 1.0;
    //tmp = inv_logit(mu_tmp[m])*2*pi();
    //mu[m,1] = cos(tmp); mu[m,2] = sin(tmp);
  }
}
  
model{
  // priors
  mu ~ multi_normal(rep_vector(0,2),diag_matrix(rep_vector(10^5,2)));
  //mu_tmp ~ normal(0,10^5);
  tau ~ cauchy(0,5);
  //tau ~ inv_gamma(0.01,0.01);
  rho ~ uniform(-1.0,1.0);
  alpha ~ dirichlet(rep_vector(2.0,M));
  for(n in 1:N){
    vector[M] ps;
    for(m in 1:M){
      ps[m] = log(alpha[m]) + pn_circle_lpdf(theta[n]|mu[m],sigma[m]);
    }
    target += log_sum_exp(ps);
  }
}
