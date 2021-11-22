library(ggplot2)
library(gridExtra)
kernal = c('matern_5_2', 'exp')
missing_per = c(0.5,0.2)
missing_pattern = c("random","disk")

results_5000_average = data.frame(kernal = rep(kernal,each=3),
                                  missing_per = rep(c(0.5,0.2,0.2),2),
                                  missing_pattern = rep(c("random","random","disk"),2),
                                  delta_RMSE = rep(0,6),
                                  delta_L= rep(0,6),
                                  delta_U= rep(0,6),
                                  GOLF_RMSE= rep(0,6),
                                  GOLF_CVG= rep(0,6),
                                  GOLF_INT= rep(0,6),
                                  Exact_GP_RMSE= rep(0,6),
                                  Exact_GP_CVG= rep(0,6),
                                  Exact_GP_INT= rep(0,6))

results_5000_average$kernal = as.character(results_5000_average$kernal)

for (index_avg in 1:6){
  kernel_type = results_5000_average$kernal[index_avg] #c('matern_5_2', 'exp')
  missing_per = results_5000_average$missing_per[index_avg] #c(0.5,0.2)
  missing_pattern = results_5000_average$missing_pattern[index_avg] #c("random","disk")
  
  m = 100
  
  results_5000_each = data.frame( random_seed = rep(0,m),delta_RMSE = rep(0,m), delta_L= rep(0,m), delta_U= rep(0,m), GOLF_RMSE= rep(0,m), GOLF_CVG= rep(0,m), GOLF_INT= rep(0,m), Exact_GP_RMSE= rep(0,m),  Exact_GP_CVG= rep(0,m),  Exact_GP_INT= rep(0,m))
  
  
  
  for (index in 1:m){
    if(index_avg==3 & index %in% c(39,44,69,77)){
      next
    }
    filename_full_GP = paste0("Parallel/" ,"size_25_",kernel_type,"_", missing_pattern ,"_",missing_per,"/", "Results_full_GP_seed_",index,".RData")
    filename_GOLF = paste0("Parallel/" ,"size_25_",kernel_type,"_", missing_pattern ,"_",missing_per,"/", "Results_fast_seed_",index,".RData")
    
    load(filename_GOLF)
    load(filename_full_GP)
    results_5000_each$random_seed[index] = index
    
    # 1: delta_rmse
    
    results_5000_each$delta_RMSE[index] = sqrt(mean((Y_avg_sample_mean[missing_index]-Y_sample_mean_avg_full_GP[missing_index] )^2 ))
    
    # 2: delta_L
    # lower_interval_95_values_full_GP=rowMaxs(lower_interval_95_full_GP,value=T)
    # upper_interval_95_values_full_GP=rowMins(upper_interval_95_full_GP,value=T)
    
    results_5000_each$delta_L[index] = mean(abs((lower_interval_95_values-lower_interval_95_values_full_GP)))
    
    # 3: delta_U
    
    results_5000_each$delta_U[index] = mean(abs((upper_interval_95_values-upper_interval_95_values_full_GP)))
    
    # 4: RMSE, CVG , INT
    
    results_5000_each$GOLF_RMSE[index] = sqrt(mean((Y_avg_sample_mean[missing_index]-Y_full[missing_index] )^2 ))
    results_5000_each$Exact_GP_RMSE[index] = sqrt(mean((Y_sample_mean_avg_full_GP[missing_index]-Y_full[missing_index] )^2 ))
    
    results_5000_each$GOLF_CVG[index] = sum( (Y_full[missing_index]<=upper_interval_95_values)& (Y_full[missing_index]>=lower_interval_95_values))/length(missing_index)
    
    results_5000_each$Exact_GP_CVG[index] = sum( (Y_full[missing_index]<=upper_interval_95_values_full_GP)& (Y_full[missing_index]>=lower_interval_95_values_full_GP))/length(missing_index)
    
    results_5000_each$GOLF_INT[index] = mean(upper_interval_95_values-lower_interval_95_values)
    
    results_5000_each$Exact_GP_INT[index] = mean(upper_interval_95_values_full_GP-lower_interval_95_values_full_GP)
    

    
    
    # 5: histograme of parameters
    
    if ((index_avg==1) & (index==1)){
      param_all = as.data.frame(rbind(param_record[1001:5000,1:4],param_record_full_GP[1001:5000,]))
      colnames(param_all) = c('beta0',  'eta1', 'beta1', 'sigma_0_2')
      param_all$method = rep(c("GOLF", "Exact GP"),each=4000 )
      
      
      plot_1 = ggplot(data = param_all, aes(x = log(beta0),group=method, color = method)) + 
        geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
        geom_density(alpha=.2) + 
        xlab( expression(log(beta[0])))
      
      # plot_1_trace = ggplot(data = param_all, aes(y = log(beta0), x = rep(1:5000,2), group=method, color = method)) + 
      #   geom_line() +  
      #   ylab( expression(log(beta[0]))) +
      #   theme(axis.title.x=element_blank())
      
      plot_2 = ggplot(data = param_all, aes(x = log(eta1),group=method, color = method)) + 
        geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
        geom_density(alpha=.2)+ 
        xlab( expression(log(eta[1])))
      
      # plot_2_trace = ggplot(data = param_all, aes(y = log(eta1), x = rep(1:5000,2), group=method, color = method)) + 
      #   geom_line() +  
      #   ylab( expression(log(eta[1]))) +
      #   theme(axis.title.x=element_blank())
      
      plot_3 = ggplot(data = param_all, aes(x = log(beta1),group=method, color = method)) + 
        geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
        geom_density(alpha=.2)+ 
        xlab( expression(log(beta[1])))
      
      # plot_3_trace = ggplot(data = param_all, aes(y = log(beta1), x = rep(1:5000,2), group=method, color = method)) + 
      #   geom_line() +  
      #   ylab( expression(log(beta[1]))) +
      #   theme(axis.title.x=element_blank())
      
      plot_4 = ggplot(data = param_all, aes(x = log(sigma_0_2),group=method, color = method)) + 
        geom_histogram(aes(y=..density..),fill="white", bins = 100,position="identity") + 
        geom_density(alpha=.2)+ 
        xlab( expression(log(sigma[0]^2)))
      
      # figure_2 = grid.arrange(plot_1, plot_2,plot_3, ncol=3, nrow = 1)
      # grid.arrange(plot_1_trace, plot_2_trace, plot_3_trace, ncol=1, nrow = 3)
    }
    
    results_5000_average[index_avg,4:12] = round(apply(results_5000_each, 2,mean)[2:10],3)
  }
}

# Figure 2
grid.arrange(plot_1, plot_2,plot_3, ncol=3, nrow = 1)

# Table 2
results_5000_average

Y_missing_disk = Y_full
Y_missing_disk[missing_index] = NA

Y_missing_random = Y_full
Y_missing_random[missing_index] = NA

par(mfrow=c(1,1))
image2D(Y_full)
image2D(Y_missing_disk)
image2D(Y_missing_random)
