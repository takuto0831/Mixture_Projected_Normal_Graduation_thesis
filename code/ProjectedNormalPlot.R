library(dplyr)
library(tidyverse)
library(foreach)
library(ggthemes)
library(circular)

# project normal func
projected_normal_circular <- function(theta,mu,Sigma){
  u = matrix(c(cos(theta),sin(theta)),ncol=1)
  A = t(u) %*% solve(Sigma) %*% u
  B = t(u) %*% solve(Sigma) %*% mu
  C = (-1/2) * (t(mu) %*% solve(Sigma) %*% mu)
  tmp = B/sqrt(A)
  p = (1/(2*pi*A*sqrt(det(Sigma)))) * exp(C) * 
    (1 + tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1))
  return(p)
}

mu = matrix(c(-0.19,2.09),ncol=1); sigma1 = 1.58; sigma2 = 1.4; rho = -0.84
Sigma = matrix(c(sigma1^2,rho*sigma1*sigma2,
                 rho*sigma1*sigma2,sigma2^2),ncol=2)

# 予測値をプロット
theta <- seq(0, 2 * pi ,0.1)
x <- foreach(i=theta, .combine = c) %do% projected_normal_circular(i,mu,Sigma) 

data.frame(theta = theta,prob = x) %>% 
  ggplot(aes(x = theta, y = prob)) + 
  geom_line() 

# project normal 2峰性

mu = matrix(c(-0.24,0.15),ncol=1); sigma1 = 0.458; sigma2 = 1; rho = 0.15
Sigma = matrix(c(sigma1^2,rho*sigma1*sigma2,
                 rho*sigma1*sigma2,sigma2^2),ncol=2)
# 予測値をプロット
theta <- seq(0, 2 * pi ,0.1)
x <- foreach(i=theta, .combine = c) %do% projected_normal_circular(i,mu,Sigma) 

data.frame(theta = theta,prob = x) %>% 
  ggplot(aes(x = theta, y = prob)) + 
  geom_line() 

# von Mises sample
x <- foreach(i=theta, .combine = c) %do% dvonmises(i, mu=circular(3/4), kappa=5)
data.frame(theta = theta,prob = x) %>% 
  ggplot(aes(x = theta, y = prob)) + 
  geom_line() 

