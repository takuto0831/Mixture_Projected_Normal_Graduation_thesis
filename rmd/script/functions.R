# projcted normal (式(1))
# サンプルデータを作成する
# 円の場合
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

# mcmcで予測したパラメータを元に分布を生成する関数
projected_normal_dist <- function(fit){
  #結果の確認
  print(fit)
  
  # 値の抽出
  fit_ext <- extract(fit,permuted=T)
  # シュミレーション値1
  mu = matrix(c(fit_ext$mu[,1] %>% mean(),fit_ext$mu[,2] %>% mean()),ncol=1)
  Sigma = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),
                   fit_ext$sigma[,2,1] %>% mean(),1),ncol=2)
  # 予測値をプロット
  theta = seq(0,2*pi,0.01)
  x <- foreach(i=theta, .combine = c) %do% projected_normal_circular(i,mu,Sigma) 
  data.frame(theta = theta,pred = x) %>%
    ggplot(aes(x=theta,y=pred)) +
    geom_line() + 
    theme_bw()
  }
  
# 分布を元にthetaの乱数を発生させる
PN_sample <- function(x,theta,iter){
  cl <- makeCluster(4)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  # 累積確率
  y <- foreach(i=1:length(x), .combine = c) %dopar% {
    sum(x[1:i])
  }
  count_func <- function(y,theta){
    tmp <- runif(1,min = min(y),max = max(y))
    for(i in 2:length(y)) 
      if(tmp > y[i-1] && tmp < y[i]) 
        return(theta[i])
  }
  
  pn_sample <- foreach(i=1:iter, .combine = c) %dopar% {
    count_func(y,theta)
  }
  return(pn_sample)
}

# 複数の方法でパラメータを予測したときそれらの結果の比較を行う関数
projected_normal_dist_comp <- function(fit, fit_){
  #結果の確認
  print(fit)
  print(fit_)
 
  # 値の抽出
  fit_ext <- extract(fit,permuted=T)
  fit_ext_ <- extract(fit_,permuted=T)
  
  # シュミレーション値1
  mu = matrix(c(fit_ext$mu[,1] %>% mean(),fit_ext$mu[,2] %>% mean()),ncol=1)
  Sigma = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),
                    fit_ext$sigma[,2,1] %>% mean(),1),ncol=2)
  # シュミレーション値2
  mu_ = matrix(c(fit_ext_$mu[,1] %>% mean(),fit_ext_$mu[,2] %>% mean()),ncol=1)
  Sigma_ = matrix(c(fit_ext_$sigma[,1,1] %>% mean(),fit_ext_$sigma[,1,2] %>% mean(),
                   fit_ext_$sigma[,2,1] %>% mean(),1),ncol=2)
  
  # 予測値をプロット
  theta = seq(0,2*pi,0.01)
  x <- foreach(i=theta, .combine = c) %do% projected_normal_circular(i,mu,Sigma) #gamma
  x_ <- foreach(i=theta, .combine = c) %do% projected_normal_circular(i,mu_,Sigma_) #half cauchy
  data.frame(theta = theta,pred_gamma = x, pred_cauchy = x_) %>%
    gather(label,pred,-theta) %>% 
    ggplot(aes(x=theta,y=pred)) +
    geom_line(aes(color = label))
}

############### 球の場合 ####################
projected_normal_sphere <- function(theta1,theta2,mu,Sigma){
  u = matrix(c(cos(theta1)*sin(theta2),sin(theta1)*sin(theta2),cos(theta2)),ncol=1)
  A = t(u) %*% solve(Sigma) %*% u
  B = t(u) %*% solve(Sigma) %*% mu
  C = (-1/2) * (t(mu) %*% solve(Sigma) %*% mu)
  D = B/sqrt(A)
  p = (1/(2*pi*A))^(1.5) * (exp(C)/sqrt(det(Sigma))) *
    ( (D*pnorm(D,0,1)/dnorm(D,0,1) + 1)*D + pnorm(D,0,1)/dnorm(D,0,1) )
  return(p)
}

# 分布を元にtheta1,theta2の乱数を発生させる
PN_sample_sphere <- function(x,theta1,theta2,iter){
  len <- length(theta1) 
  ans <- c() # 累積確率を格納すす
  
  count_func <- function(data,theta1,theta2){
    options(digits=10)
    tmp <- runif(1,min = min(data$z),max = max(data$z))
    for(i in 2:length(data$z))  
      if(tmp > data$z[i-1] && tmp <= data$z[i]) 
        return(c(theta1[data$x[i]],theta2[data$y[i]]))
  }
  
  # 並列処理
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  # 累積確率
  for(i in 1:len){
    foreach(j=1:len, .combine = c) %dopar% {
      options(digits=10)
      sum(x[1:i,1:j])
    } -> tmp
    rbind(ans,tmp) -> ans
  }
  data.frame(x = rep(seq(1,len,1),len),
             y = rep(1:len, each=len),
             z = ans %>% as.vector()) -> data
  data <- data[order(data$z),]
  
  # theta1,theta2の行列を作る
  pn_sample <- foreach(i=1:iter, .combine = "rbind") %dopar% {
    count_func(data,theta1,theta2)
  }
  stopCluster(cl)
  return(pn_sample %>% as.matrix())
}

# mcmcで予測したパラメータを元に分布を生成する関数
projected_normal_sphere_dist <- function(fit,theta1,theta2,func){
  #結果の確認
  print(fit)
  # 値の抽出
  fit_ext <- extract(fit,permuted=T)
  # シュミレーション値1
  mu = matrix(c(fit_ext$mu[,1] %>% mean(),fit_ext$mu[,2] %>% mean(),fit_ext$mu[,3] %>% mean()),ncol=1)
  ###### ここから作る #####
  Sigma = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),fit_ext$sigma[,1,3] %>% mean(),
                   fit_ext$sigma[,2,1] %>% mean(),fit_ext$sigma[,2,2] %>% mean(),fit_ext$sigma[,2,3] %>% mean(),
                   fit_ext$sigma[,3,1] %>% mean(),fit_ext$sigma[,3,2] %>% mean(),fit_ext$sigma[,3,3] %>% mean()),ncol = 3)
  #star = matrix(c(fit_ext$star[,1,1] %>% mean(),fit_ext$star[,1,2] %>% mean(),
  #               fit_ext$star[,2,1] %>% mean(),fit_ext$star[,2,2] %>% mean()),ncol = 2)

  # 並列処理
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  x <- c()
  for(i in theta1){
    foreach(j=theta2, .combine = cbind) %dopar%{
      func(i,j,mu,Sigma)
    } -> tmp 
    rbind(x,tmp) -> x
  } 
  stopCluster(cl)
  return(x)
}