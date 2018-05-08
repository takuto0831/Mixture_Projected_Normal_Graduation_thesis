ChangeRange <- function(data,num){
  #num = 1のとき(-pi~pi を 0~2*pi に変更)
  if(num == 1){
    data[data<0] + 2*pi -> data[data<0]
    return(data)
  }
  #num = 2のとき(0~2*pi を -pi~pi に変更)
  if(num == 2){
    data[data>pi] - 2*pi -> data[data>pi]
    return(data)
  }
} 

CheckData <- function(data,label){
  # それぞれの分布
  data.frame(data,z = as.factor(label)) %>% 
    ggplot(aes(x=data,fill=z)) +
    geom_density(alpha = 0.3) +
    theme_economist() -> p
  print(p)
  
  # 混合分布
  data.frame(data) %>% 
    ggplot(aes(x=data)) +
    geom_density(alpha = 0.3) +
    theme_economist() -> p
  print(p)
}

# 混合分布のそれぞれを分布する関数

VonMisesPlot <- function(fit,clust){
  # パラメータ抽出
  fit_von <- extract(fit,permuted=T); para <- list(); theta = seq(0,2*pi,0.01); plots <- list(); data <-c()
  # 平均と分散抽出
  for(i in 1:clust){
    para[[i]] = fit_von$mu[,i] %>% mean() 
    para[[i+clust]] = fit_von$kappa[,i] %>% mean() 
  }
  for(i in 1:clust){
    x <- foreach(j=theta, .combine = c) %do% dvonmises(circular(j),circular(para[[i]]),para[[i+clust]])
    # 推定した混合比率を用いて分布作成
    data <- rbind(data,c(x) * (fit_von$alpha[,i] %>% mean()))
    data.frame(theta = theta,pred = x) %>%
      ggplot(aes(x=theta,y=pred)) +
      geom_line() +  labs(title = i) -> plots[[i]]
  }  
  multiplot(plotlist = plots,cols = 3)
  
  # 混合分布
  data.frame(density = apply(data,2,sum),theta) %>% 
    ggplot(aes(x=theta,y=density)) +
    geom_line() +
    theme_economist() -> p
  print(p)
  }

PNPlot <- function(fit,clust){
  # パラメータ抽出
  fit_pn <- extract(fit,permuted=T); para <- list(); theta = seq(0,2*pi,0.01); plots <- list();data <-c()
  # 平均と分散抽出
  for(i in 1:clust){
    para[[i]] = matrix(c(fit_pn$mu[,i,1] %>% mean(),fit_pn$mu[,i,2] %>% mean()),ncol=1)
    para[[i+clust]] = matrix(c(fit_pn$sigma[,i,1,1] %>% mean(),fit_pn$sigma[,i,1,2] %>% mean(),
                               fit_pn$sigma[,i,2,1] %>% mean(),1),ncol=2)
  }
  for(i in 1:clust){
    x <- foreach(j=theta, .combine = c) %do% projected_normal_circular(j,para[[i]],para[[i+clust]]) 
    # 推定した混合比率を用いて分布作成
    data <- rbind(data,c(x) * (fit_pn$alpha[,i] %>% mean()))
    data.frame(theta = theta,pred = x) %>%
      ggplot(aes(x=theta,y=pred)) +
      geom_line() + labs(title = i) -> plots[[i]]
  }
  multiplot(plotlist = plots,cols = 3)
  # 混合分布
  data.frame(density = apply(data,2,sum),theta) %>% 
    ggplot(aes(x=theta,y=density)) +
    geom_line() +
    theme_economist() -> p
  print(p)
}

# 混合分布を用いてクラスタリングする.
# あらゆるクラスター数に対応させる
ClusterVonMises <- function(data,fit,clust){
  # パラメータ抽出
  fit_von <- extract(fit,permuted=T); para <- list(); 
  # データ定義
  data_len <- length(data); theta = seq(0,2*pi,0.01); ANS <- matrix(0, nrow=data_len, ncol=clust);
  # 平均と分散抽出
  for(i in 1:clust){
    para[[i]] = fit_von$mu[,i] %>% mean() 
    para[[i+clust]] = fit_von$kappa[,i] %>% mean() 
  }
  for(i in 1:data_len){
    # 予測ラベルを格納する
    ans <- c(); tmp <- c();
    for(j in 1:length(data[[i]])){
      for(k in 1:clust){
        tmp[k] <- dvonmises(circular(data[[i]][j]),circular(para[[k]]),para[[k+clust]])
      }
      ans[j] <- which.max(tmp %>% as.vector())
    }
    table(ans) %>% as.data.frame()-> tmp
    for(m in 1:dim(tmp)[1]){
      ANS[i,tmp[m,1] %>% as.character() %>% as.numeric()] <- tmp[m,2]
    }
    #which.max(tmp$Freq) %>% tmp$Freq[.]/sum(tmp$Freq) -> ANS[i]
  }
  return(ANS)
}

ClusterProjectedNormal <- function(data,fit,clust){
  # パラメータ抽出
  fit_pn <- extract(fit,permuted=T); para <- list(); 
  # データ定義
  data_len <- length(data); theta = seq(0,2*pi,0.01);ANS <- matrix(0, nrow=data_len, ncol=clust);
  # 平均と分散抽出
  for(i in 1:clust){
    para[[i]] = matrix(c(fit_pn$mu[,i,1] %>% mean(),fit_pn$mu[,i,2] %>% mean()),ncol=1)
    para[[i+clust]] = matrix(c(fit_pn$sigma[,i,1,1] %>% mean(),fit_pn$sigma[,i,1,2] %>% mean(),
                                 fit_pn$sigma[,i,2,1] %>% mean(),1),ncol=2)
  }
  for(i in 1:data_len){
    # 予測ラベルを格納する
    ans <- c(); tmp <- c();
    for(j in 1:length(data[[i]])){
      for(k in 1:clust){
        tmp[k] <-projected_normal_circular(data[[i]][j],para[[k]],para[[k+clust]])
      }
      ans[j] <- which.max(tmp)
    }
    table(ans) %>% as.data.frame() -> tmp
    for(m in 1:dim(tmp)[1]){
      ANS[i,tmp[m,1] %>% as.character() %>% as.numeric()] <- tmp[m,2]
    }
    #which.max(tmp$Freq) %>% tmp$Freq[.]/sum(tmp$Freq) -> ANS[i]
  }
  return(ANS)
}
