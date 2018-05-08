########### 必要なパッケージ ###########
rm(list=ls())
library(mvtnorm)
library(rstan)
library(tidyverse)
library(ggplot2)
library(MCMCpack)
library(loo)
library(rgl)
########### 関数 #############
ClusterProjectedNormal <- function(data,fit,clust){
  fit_pn <- rstan::extract(fit,permuted=T); para <- list(); data_len <- dim(data)[1]; ans <- c();
  # 平均と分散抽出
  for(i in 1:clust){
    para[[i]] = matrix(c(fit_pn$mu[,i,1] %>% mean(),fit_pn$mu[,i,2] %>% mean(),fit_pn$mu[,i,3] %>% mean()), ncol=1)
    para[[i+clust]] = matrix(c(fit_pn$sigma[,i,1,1] %>% mean(),fit_pn$sigma[,i,1,2] %>% mean(),fit_pn$sigma[,i,1,3] %>% mean(),
                               fit_pn$sigma[,i,2,1] %>% mean(),fit_pn$sigma[,i,2,2] %>% mean(),fit_pn$sigma[,i,2,3] %>% mean(),
                               fit_pn$sigma[,i,3,1] %>% mean(),fit_pn$sigma[,i,3,2] %>% mean(),1),ncol=3)}
  for(i in 1:data_len){
    tmp <- c() # 予測ラベルを格納する
    for(k in 1:clust){tmp[k] <-projected_normal_sphere(as.numeric(data[i,1]),as.numeric(data[i,2]),para[[k]],para[[k+clust]])}
    ans[i] <- which.max(tmp)
  }
  return(ans)}

MakeData <- function(iter,mu,Sigma,number){
  set.seed(831)
  rmvnorm(iter,mu,Sigma) %>% 
    as.data.frame() %>% as.tbl() %>% 
    mutate(r = sqrt(V1^2 + V2^2 + V3^2)) %>% 
    mutate(x = V1/r, y = V2/r, z = V3/r) %>% 
    mutate(theta2 = acos(z)) %>% 
    mutate(theta1 = ifelse(y/sin(theta2) > 0, acos(x/sin(theta2)), - acos(x/sin(theta2)))) %>% 
    mutate(theta1 = ifelse(theta1 < 0, theta1 + 2*pi, theta1)) %>% 
    mutate(label = number) %>% 
    dplyr::select(theta1,theta2,label,x,y,z) -> data
  return(data)
}
MakeParameterList <- function(fit,clust){
  # data : 列 theta1 , theta2 の (? * 2) matrix
  fit_pn <- rstan::extract(fit,permuted=T); para <- c(); 
  # 平均と分散抽出
  for(i in 1:clust){
    tmp = matrix(c(fit_pn$mu[,i,1] %>% mean(),fit_pn$mu[,i,2] %>% mean(),fit_pn$mu[,i,3] %>% mean()), nrow=1)
    #para[[i+clust]] = matrix(c(fit_pn$sigma[,i,1,1] %>% mean(),fit_pn$sigma[,i,1,2] %>% mean(),fit_pn$sigma[,i,1,3] %>% mean(),
    #                           fit_pn$sigma[,i,2,1] %>% mean(),fit_pn$sigma[,i,2,2] %>% mean(),fit_pn$sigma[,i,2,3] %>% mean(),
    #                           fit_pn$sigma[,i,3,1] %>% mean(),fit_pn$sigma[,i,3,2] %>% mean(),1),ncol=3)
    para <- rbind(para,tmp)
  }
  return(as.data.frame(para))
}

########## 乱数生成 ############
iter=500; 
set.seed(1)
mu1 = matrix(c(1,0,0),ncol=1); Sigma1 = riwish(4,matrix(c(1,0,0,0,1,0,0,0,1),ncol = 3));
set.seed(2)
mu2 = matrix(c(0,0.5,2),ncol=1); Sigma2 = riwish(4,matrix(c(1,0,0,0,1,0,0,0,1),ncol = 3));
set.seed(5)
mu3 = matrix(c(1.5,1,0),ncol=1); Sigma3 = riwish(4,matrix(c(1,0,0,0,1,0,0,0,1),ncol = 3));
set.seed(4)
mu4 = matrix(c(1,1,1),ncol=1); Sigma4 = riwish(4,matrix(c(1,0,0,0,1,0,0,0,1),ncol = 3));

########## 角度データ, 正解ラベル取得 #############
data1 <- MakeData(iter,mu1,Sigma1,1)
data2 <- MakeData(iter,mu2,Sigma2,2)
data3 <- MakeData(iter,mu3,Sigma3,3)
data4 <- MakeData(iter,mu4,Sigma4,4)
all_data <- rbind(data1,data2,data3,data4) 

########## 正解ラベルの確認 ##############
all_data%>% 
  ggplot(aes(x= theta2,y = theta1)) + 
  geom_point(aes(colour = factor(label))) +
  labs(shape="label",colour="label") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )
################################## 並列処理用 ##############################
#rstan_options(auto_write=TRUE)
#options(mc.cores=parallel::detectCores())
fit3 <- stan("stan/mix_pn_sphere.stan", 
            data=list(N=dim(all_data)[1],M=3,theta1=all_data$theta1,theta2=all_data$theta2),
            iter = 20000,chains = 1,open_progress = FALSE)
fit4 <- stan("stan/mix_pn_sphere.stan", 
             data=list(N=dim(all_data)[1],M=4,theta1=all_data$theta1,theta2=all_data$theta2),
             iter = 10000,chains = 1,open_progress = FALSE)
fit5 <- stan("stan/mix_pn_sphere.stan", 
             data=list(N=dim(all_data)[1],M=5,theta1=all_data$theta1,theta2=all_data$theta2),
             iter = 20000,chains = 1,open_progress = FALSE)

############ データ保存 ##############
save(fit3,file="Rdata/sim3.Rdata")
save(fit4,file="Rdata/sim4.Rdata")
save(fit5,file="Rdata/sim5.Rdata")

####################### 結果の確認 ##############################
all(summary(fit3)$summary[,"Rhat"] <= 1.10, na.rm=T)
all(summary(fit4)$summary[,"Rhat"] <= 1.10, na.rm=T)
all(summary(fit5)$summary[,"Rhat"] <= 1.10, na.rm=T)

############  mcmc 更新 ##########
dat <- MakeParameterList(fit4,4)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

fit4_ <- stan("stan/mix_pn_sphere_first_point.stan", 
                 data=list(N=dim(all_data)[1],M=4,theta1=all_data$theta1,
                           theta2=all_data$theta2,mu1=dat),
                 iter = 10000,
                 chains = 4,
                 open_progress = FALSE)

####################### WAICの確認 ##############################
log_lik3 <- extract_log_lik(fit3,"log_likelihood")
loo::waic(log_lik3)
log_lik4 <- extract_log_lik(fit4,"log_likelihood")
loo::waic(log_lik4)
log_lik5 <- extract_log_lik(fit5,"log_likelihood")
loo::waic(log_lik5)

########################  クラスター分析 ##########################
all_data$clus3 <- ClusterProjectedNormal(all_data,fit3,3)
all_data$clus4 <- ClusterProjectedNormal(all_data,fit4,4)
all_data$clus5 <- ClusterProjectedNormal(all_data,fit5,5)

#### plot meikakuni
all_data %>% 
  mutate(clus4_ = ifelse(clus4==1,3,
                        ifelse(clus4==2,1,
                               ifelse(clus4==3,4,2)))) -> all_data



all_data%>% 
  ggplot(aes(x= theta2,y = theta1)) + 
  geom_point(aes(colour = factor(clus3))) +
  labs(shape="label",colour="label") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

all_data%>% 
  ggplot(aes(x= theta2,y = theta1)) + 
  geom_point(aes(colour = factor(clus4_))) +
  labs(shape="label",colour="label") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

all_data%>% 
  ggplot(aes(x= theta2,y = theta1)) + 
  geom_point(aes(colour = factor(clus5)))

################### 正解ラベルとの比較 ########################
table(all_data$label,all_data$clus4_)

#################### 真の平均方向 ######################
a1 <- matrix(c(0,0,0),ncol=3)
a2 <- matrix(c(1,0,0),ncol=3)
a3 <- matrix(c(0,0.5,2),ncol=3)
a4 <- matrix(c(1.5,1,0),ncol=3)
a5 <- matrix(c(1,1,1),ncol=3)

plot3d(all_data$x,all_data$y,all_data$z,col = rainbow(4)[factor(all_data$label)],
       size = 3,xlab = "x",ylab = "y",zlab = "z")

arrow3d(a1,a2/(a2^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)
arrow3d(a1,a3/(a3^2 %>% sum() %>% sqrt()),type = "fslat",thickness = 1)
arrow3d(a1,a4/(a4^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)
arrow3d(a1,a5/(a5^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)

#################### 予測の平均方向 ######################
a1 <- matrix(c(0,0,0),ncol=3)
a2 <- matrix(c(1.52,1.65,-0.04),ncol=3)
a3 <- matrix(c(1.98,-0.15,-0.44),ncol=3)
a4 <- matrix(c(0.66,0.58,0.55),ncol=3)
a5 <- matrix(c(0.00,0.10,0.44),ncol=3)

plot3d(all_data$x,all_data$y,all_data$z,col = rainbow(4)[factor(all_data$clus4_)],
       size = 3,xlab = "x",ylab = "y",zlab = "z")

arrow3d(a1,a2/(a2^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)
arrow3d(a1,a3/(a3^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)
arrow3d(a1,a4/(a4^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)
arrow3d(a1,a5/(a5^2 %>% sum() %>% sqrt()),type = "flat",thickness = 1)

