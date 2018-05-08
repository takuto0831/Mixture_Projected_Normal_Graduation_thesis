setwd("C:/Users/SHIO-160412-4/Desktop/projected normal")
source('script/functions_mix_dist.R', encoding = 'UTF-8')
source('script/functions.R', encoding = 'UTF-8')

ans_box <- c()
for(i in 1:20){
  #真のクラスターが3つの場合に予測クラスターを2つとして調べる
  rstan_options(auto_write=TRUE)
  #options(mc.cores=parallel::detectCores())

  fit2 <- stan("stan/mix_pn_circular.stan", 
             data=list(N=dat,M=2,theta=data),
             iter = 20000,
             chains = 1,
             open_progress = FALSE)

  # 結果の確認
  if(all(summary(fit2)$summary[,"Rhat"] <= 1.10, na.rm=T)){
    stan_trace(fit2)
    stan_ac(fit2)

    # stan file の保存 読み込み
    filename <- sprintf("stan_fit/pn_model2_%02d.rda",i)
    save(fit2, file=filename)

    # 結果の確認
    options(max.print=1000)
    print(fit2)
    PNPlot(fit2,2) %>% 
      multiplot(plotlist = .,cols = 3)

    # クラスタリング分析
    ClusterProjectedNormal(data_list,fit2,2) %>% 
      rbind(ans_box,.) -> ans_box
  }
}