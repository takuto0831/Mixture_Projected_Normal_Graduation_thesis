setwd("C:/Users/SHIO-160412-4/Desktop/projected normal")
source('script/functions_mix_dist.R', encoding = 'UTF-8')
source('script/functions.R', encoding = 'UTF-8')

ans_box <- c()
for(i in 1:20){
  #真のクラスターが3つの場合に予測クラスターを3つとして調べる
  rstan_options(auto_write=TRUE)
  #options(mc.cores=parallel::detectCores())
  fit3 <- stan("stan/mix_pn_circular.stan", 
               data=list(N=dat,M=3,theta=data),
               iter = 20000,
               chains = 1,
               open_progress = FALSE)

  # 結果の確認
  if(all(summary(fit3)$summary[,"Rhat"] <= 1.10, na.rm=T)){
    stan_trace(fit3)
    stan_ac(fit3)
    # stan file の保存 読み込み
    filename <- sprintf("stan_fit/pn_model3_%02d.rda",i)
    save(fit3,  file=filename)
    # 結果の確認
    options(max.print=1000)
    print(fit3)
    PNPlot(fit3,3) %>% 
      multiplot(plotlist = .,cols = 3)
    # クラスタリング分析
    ClusterProjectedNormal(data_list,fit3,3) %>% 
      rbind(ans_box,.) -> ans_box
  }
}
