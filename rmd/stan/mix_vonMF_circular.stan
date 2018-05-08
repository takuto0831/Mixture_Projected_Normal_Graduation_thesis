data{
  int<lower=0> N;
  int<lower=0> M;
  real theta[N];
}

parameters{
  // ordered[M] mu_tmp;
  vector<lower=0,upper=2*pi()>[M] mu;
  vector<lower=0,upper=100>[M] kappa;
  simplex[M] alpha;
}

//transformed parameters{
//  vector<lower=0,upper=2*pi()>[M] mu;
  //inv_logit で0~1にして2piをかける
//  for(m in 1:M){
//    mu[m] = inv_logit(mu_tmp[m])*2*pi();
//  }
//}

model{
  mu ~ normal(0,10^5);
  //mu_tmp ~ normal(0,10^5);
  kappa ~ exponential(0.05);
  alpha ~ dirichlet(rep_vector(2.0,M));
  //alpha ~ beta(1,1); //混合比率の事前分布
  for(n in 1:N){
    vector[M] ps;
    for(m in 1:M){
      ps[m] = log(alpha[m]) + von_mises_lpdf(theta[n]|mu[m],kappa[m]);
    }
    target += log_sum_exp(ps);
  }
}
