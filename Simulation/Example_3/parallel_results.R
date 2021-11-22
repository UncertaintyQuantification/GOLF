library(ggplot2)

kernel_type = 'exp' #c('matern_5_2', 'exp')
missing_per = 0.2 #c(0.5,0.2)
missing_pattern = "disk" #c("random","disk")

m = 100
figure = FALSE

results_5000 = data.frame( random_seed = rep(0,m),delta_RMSE = rep(0,m), delta_L= rep(0,m), delta_U= rep(0,m), fast_RMSE= rep(0,m), slow_RMSE= rep(0,m), fast_CVG= rep(0,m), slow_CVG= rep(0,m), fast_INT= rep(0,m), slow_INT= rep(0,m))



for (index in 1:m){
  filename_full_GP = paste0("Parallel/" ,"size_25_",kernel_type,"_", missing_pattern ,"_",missing_per,"/", "Results_full_GP_seed_",index,".RData")
  filename_fast = paste0("Parallel/" ,"size_25_",kernel_type,"_", missing_pattern ,"_",missing_per,"/", "Results_fast_seed_",index,".RData")
  
  load(filename_fast)
  load(filename_full_GP)
  results_5000$random_seed[index] = index
  
  # 1: delta_rmse
  
  results_5000$delta_RMSE[index] = sqrt(mean((Y_avg_sample_mean[missing_index]-Y_sample_mean_avg_full_GP[missing_index] )^2 ))
  
  # 2: delta_L
  # lower_interval_95_values_full_GP=rowMaxs(lower_interval_95_full_GP,value=T)
  # upper_interval_95_values_full_GP=rowMins(upper_interval_95_full_GP,value=T)
  
  results_5000$delta_L[index] = sqrt(mean((lower_interval_95_values-lower_interval_95_values_full_GP)^2))
  
  # 3: delta_U
  
  results_5000$delta_U[index] = sqrt(mean((upper_interval_95_values-upper_interval_95_values_full_GP))^2)
  
  # 4: RMSE, CVG , INT
  
  results_5000$fast_RMSE[index] = sqrt(mean((Y_avg_sample_mean[missing_index]-Y_full[missing_index] )^2 ))
  results_5000$slow_RMSE[index] = sqrt(mean((Y_sample_mean_avg_full_GP[missing_index]-Y_full[missing_index] )^2 ))
  
  results_5000$fast_CVG[index] = sum( (Y_full[missing_index]<=upper_interval_95_values)& (Y_full[missing_index]>=lower_interval_95_values))/length(missing_index)
  
  results_5000$slow_CVG[index] = sum( (Y_full[missing_index]<=upper_interval_95_values_full_GP)& (Y_full[missing_index]>=lower_interval_95_values_full_GP))/length(missing_index)
  
  results_5000$fast_INT[index] = mean(upper_interval_95_values-lower_interval_95_values)
  
  results_5000$slow_INT[index] = mean(upper_interval_95_values_full_GP-lower_interval_95_values_full_GP)

  
  # 5: histograme of parameters

  if (figure){
    param_all = as.data.frame(rbind(param_record[,1:4],param_record_full_GP))
    colnames(param_all) = c('beta',  'eta', 'beta1', 'sigma_0_2')
    param_all$method = rep(c("fast", "full GP"),each=5000 )
    
    
    plot_1 = ggplot(data = param_all, aes(x = log(beta),group=method, color = method)) + 
      geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
      geom_density(alpha=.2) + 
      xlab( expression(log(beta)))
    
    plot_2 = ggplot(data = param_all, aes(x = log(eta),group=method, color = method)) + 
      geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
      geom_density(alpha=.2)+ 
      xlab( expression(log(eta)))
    
    plot_3 = ggplot(data = param_all, aes(x = log(beta1),group=method, color = method)) + 
      geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
      geom_density(alpha=.2)+ 
      xlab( expression(log(beta[1])))
    
    plot_4 = ggplot(data = param_all, aes(x = log(sigma_0_2),group=method, color = method)) + 
      geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
      geom_density(alpha=.2)+ 
      xlab( expression(log(sigma[0]^2)))

    grid.arrange(plot_1, plot_2,plot_3, plot_4, ncol=2, nrow = 2)
  }
  
  
}

round(results_5000,3)
