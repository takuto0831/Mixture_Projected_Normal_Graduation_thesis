functions{
  real pn_circle_lpdf(real theta,vector mu,matrix sigma){
    vector[2] u;
    real A;
    real B;
    real C;
    real D;
    real p;
    u[1] = cos(theta);
    u[2] = sin(theta);
    A = u' * inverse(sigma) * u;
    B = u' * inverse(sigma) * mu;
    C = (-0.5) * (mu' * inverse(sigma) * mu);
    D= B/sqrt(A);
    p = -log(A)-log(sqrt(determinant(sigma))) + C
    + log(1+(D * normal_cdf(D,0,1)/exp(normal_lpdf(D|0,1))));    
    return p;
  }
}

data{
  int N; //sample size
  real<lower=0,upper=2*pi()> theta[N]; //data
}

parameters{
  vector[2] mu;
  real<lower = 0> tau;
  real rho;
}

transformed parameters{
  matrix[2,2] sigma;
  sigma[1,1] = tau; sigma[1,2] = sqrt(tau)*rho;
  sigma[2,1] = sqrt(tau)*rho; sigma[2,2] = 1.0;
}

model{
  mu ~ multi_normal(rep_vector(0,2),diag_matrix(rep_vector(10^5,2)));
  tau ~ cauchy(0,5);
  rho ~ uniform(-1.0,1.0);
  for(n in 1:N){
    theta[n] ~ pn_circle(mu,sigma);
  }
}
