rm(list = ls(all.names = TRUE))
source('../../functions/functions.R')
rfun_names = ls()
library(foreach)
library(doParallel)
library(ggplot2)
library(plot3D)
library(RobustGaSP)
library(FastGaSP)
library(Rcpp)
library(RcppEigen)
sourceCpp(file='../../src/functions_Dec_2019.cpp')
library(maps)#Install maps package if not done before
library(rARPACK)
library( Rfast)
library(gridExtra)
package_names = .packages(F)

disk_missing = function(center, radius){
  m = matrix(0,n1,n2)
  g = expand.grid(1:nrow(m), 1:nrow(m))
  g$d2 = sqrt((g$Var1-center[1])^2 + (g$Var2-center[2])^2)
  g$inside = g$d2<=radius
  return(which(g$inside==T))
}

Accept_proposal = function(r){
  if(r>=1){
    return(TRUE)
  }else{
    u = runif(1,0,1)
    if (u<r){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}

# creat 6 folders to save parallel results for different parameter settings
creat_dir = FALSE # default setting is false. Switch it to TRUE if you want to creat new folders
if (creat_dir){
  dir.create("./result")
  dir.create("./result/Parallel")
  kernal = c('matern_5_2', 'exp')
  missing_per = c(0.5,0.2) 
  missing_pattern = c("random","disk")
  
  items = data.frame(kernal = rep(kernal,3), missing_per = c(rep(missing_per,each=2),0.2,0.2), missing_pattern = rep(missing_pattern,each=4)[1:6])
  
  for (i in 1:dim(items)[1]){
    filename = paste0("./result/Parallel/size_25_",as.character(items$kernal[i]),"_", items$missing_pattern[i] ,"_",items$missing_per[i],"/")
    dir.create(filename)
  }
}

# select kernels and missing values
kernel_type = 'exp' # 'matern_5_2', 'exp'
missing_per = 0.2 # 0.5, 0.2
missing_pattern = "random" # "random", "disk"

# generate m * 625 data matrix
m = 20 # number of simulations
M = 5000 # number of MCMC samples
M_0 = 1000 # burn-in samples
n1 = 25 # data size
n2 = 25
N = n1*n2
Y_truth_all = matrix(0, m, N) # m * 625 matrix
Y_full_all = matrix(0, m, N) # m * 625 matrix
Y_all = matrix(0, m, N) # m * 625 matrix
missing_index_all = list() # conatins all missing indices for different random seed

##let me do equally space first
input1=seq(0,1,1/(n1-1)) ##row
input2=seq(0,1,1/(n2-1)) ##column

input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)

input=cbind(rep(input1,n2),input_2_vec)

if(kernel_type=='matern_5_2'){
  beta=c(1,3)
  R1=matern_5_2_funct(abs(outer(input1,input1,'-')),beta[1])
  R2=matern_5_2_funct(abs(outer(input2,input2,'-')),beta[2])
}else if(kernel_type=='exp'){
  beta=c(1/4,1)
  R1=pow_exp_funct(abs(outer(input1,input1,'-')),beta[1],alpha_i=1)
  R2=pow_exp_funct(abs(outer(input2,input2,'-')),beta[2],alpha_i=1)
  
}
L1=t(chol(R1))
L2=t(chol(R2))

# generating m * 625 matrix, each row refers to a unique random seed from 1 to m
for (i in 1:m){
  set.seed(i)
  
  Z=matrix(rnorm(n1*n2),n1,n2)
  Z2=t(L2%*%t(Z))
  eta=0.1
  
  
  Y_truth_all[i,]=L1%*%Z2
  Y_full_all[i,]=L1%*%Z2+eta*matrix(rnorm(n1*n2),n1,n2)
  
  if (missing_pattern=="random"){
    missing_index=sample(1:(n1*n2),size=missing_per*n1*n2)
  } else if (missing_pattern=="disk"){
    if (missing_per==0.5){
      radius = 10
      missing_index = disk_missing(c(12,12),radius)
    }
    if (missing_per==0.2){
      radius = 7
      missing_index = disk_missing(c(12,12),radius)
    }
  }
  
  Y_all[i,]=Y_full_all[i,]
  Y_all[i,missing_index]=NA
  missing_index_all[[i]] = missing_index
}
#image2D(matrix(as.vector(Y_all[3,]),nrow=25,ncol=25))





## jointly robust prior
a=1/2-1

b=1

##number of d
d=n1

param_record=matrix(0,M,4) ## beta,  eta, beta1, sigma_2_0                                                                                                                                                                                                                                                                           


CL = rep(0,2)    ###CL is also used in the prior so I make it a model parameter


CL[1] = (max(input[,1])-min(input[,1]))/n1
CL[2] = (max(input[,2])-min(input[,2]))/n2


##step size for MCMC
sd_log_beta=0.05
sd_log_eta=0.2
sd_log_beta1=0.05

index_print=0

##################################################
# for loop for fast method
##################################################

for(index in 1:m){
  set.seed(index)
  Y_truth = matrix(Y_truth_all[index,],25,25)
  Y_full = matrix(Y_full_all[index,],25,25)
  Y = matrix(Y_all[index,],25,25)
  n_missing = sum(is.na(Y))
  missing_index = missing_index_all[[index]]
  
  Y_sample=Y
  
  Y_sample[missing_index]= mean(Y_sample[-missing_index])
  
  Y_sample_1=Y_sample
  
  ##initial value
  log_beta1_cur=0
  log_beta_eta_cur=c(0,-3)  ##this is for the second one
  
  
  delta_x_2=input2[2:n2]-input2[1:(n2-1)]
  index_obs= as.integer(rep(1,n2))
  
  R0_1=abs(outer(input1,input1,'-'))
  ##already defined
  #kernel_type='exp'
  
  eigen_R1_cur=Get_eigen_1(R0_1=R0_1,log_beta1=log_beta1_cur,kernel_type=kernel_type,d=d)
  ###this change a little as now it start from d
  tilde_eta_vec_cur=exp(log_beta_eta_cur[2])/eigen_R1_cur$values[1:d]
  
  Cov_2_cur=as.list(1:d)
  
  for(i_d in 1:d){
    log_beta_eta_tilde=c(log_beta_eta_cur[1],log(tilde_eta_vec_cur[i_d]))
    Cov_2_cur[[i_d]]=Construct_Cov_KF(log_beta_eta_tilde,delta_x_2,index_obs,kernel_type)
    
    
  }
  
  
  Cov_2_propose=as.list(1:d)
  
  C_R_K_Q_propose=list()
  
  # ######zero mean first #############
  # 
  # H2=cbind(rep(1,n2))
  # q2=dim(H2)[2]
  # 
  # B2_cur=matrix(0,q2,n1)
  # tilde_B2_sample_0=matrix(0,dim(B2_cur)[1],d)
  # tilde_B2_hat=tilde_B2_sample_0
  # 
  # 
  # B2_vec_hat_record=matrix(0,M,q2*n1)
  # B2_vec_cur_record=matrix(0,M,q2*n1)
  # 
  # B2_hat=matrix(0,q2,n1)
  # mean_cur2=t(B2_cur)%*%t(H2)
  # 
  # H2_t_H2_inv_H2_t=solve(t(H2)%*%H2)%*%t(H2)
  # H2_t_H2_inv=solve(t(H2)%*%H2)
  # L_H2_t_H2_inv=t(chol(H2_t_H2_inv))
  # I_H2_H2_t_H2_inv_H2_t=diag(n2)-H2%*%H2_t_H2_inv_H2_t
  #mean_cur=mean_cur2
  
  mean_cur=matrix(0,n1,n2)
  
  Uz_hat_sum=0
  index_print=0
  
  
  
  ###U1_cur is the first d eigenvector of R1
  ##I start with zero mean so the first sample has zero mean
  Y_sample_tilde_cur=t(eigen_R1_cur$vectors[,1:d])%*%(Y_sample-mean_cur)
  
  Trace_Y_sample_2=0
  
  for(i_row in 1:n1){
    Trace_Y_sample_2=Trace_Y_sample_2+sum((Y_sample[i_row,]-mean_cur[i_row,])*(Y_sample[i_row,]-mean_cur[i_row,]) )
  }
  
  
  
  Trace_Y_sample_tilde_2_cur=0
  
  
  for(i_d in 1:d){
    Trace_Y_sample_tilde_2_cur=Trace_Y_sample_tilde_2_cur+sum(Y_sample_tilde_cur[i_d,]*Y_sample_tilde_cur[i_d,])
    
  }
  
  
  ##starting value
  sigma_2_0_cur=0.1
  
  beta_post_cur=matrix(0,d,3)
  colnames(beta_post_cur)=c('log_lik','S_2','log_det')
  beta_post_propose=matrix(0,d,3)
  colnames(beta_post_propose)=c('log_lik','S_2','log_det')
  
  for(i_d in 1:d){
    
    log_beta_eta=c(log_beta_eta_cur[1],log_beta_eta_cur[2])
    
    #beta_post_S_2_log_det_propose=Get_beta_post_S_2_log_det(t(Y_sample_tilde_cur[i_row,] ),log_beta_eta_propose,tilde_eta_vec_propose[i_row],Cov_2_propose)
    
    beta_post_S_2_log_det=Get_beta_post_S_2_log_det(t(Y_sample_tilde_cur[i_d,] ),log_beta_eta,tilde_eta_vec_cur[i_d],Cov_2_cur[[i_d]],add_log_prior = F)
    beta_post_cur[i_d,]=beta_post_S_2_log_det
  }
  
  
  z_hat_matrix=matrix(0,d,n2)
  z_sample_matrix=matrix(0,d,n2)
  
  ###################################
  
  Uz_hat_sum=0
  Uz_sample=0
  index_print=0
  
  ##record some RMSE
  Y_pred_sum_ori=0
  Y_pred_sum_ori_sample_mean=0
  
  Y_sample_mean_sum=0
  
  record_RMSE_ori=rep(NA,M)
  record_RMSE_ind=rep(NA,M)
  
  record_RMSE_ori_sample_mean=rep(NA,M)
  record_RMSE_ind_sample_mean=rep(NA,M)
  
  record_log_post_beta1=rep(0,M)
  record_log_post_beta=matrix(0,M,d)
  
  #M_0=M/2
  M_0=1000
  
  record_RMSE_hat_mean=rep(NA,M)
  
  ##record Y_sample_mean_all
  Y_sample_mean_all_record=matrix(0,20,n1*n2)
  
  
  
  M_025=0.025*(M-M_0)
  
  lower_interval_95=matrix(0,length(missing_index),M_025)
  upper_interval_95=matrix(0,length(missing_index),M_025)
  
  sample_difference=rep(NA,M)
  
  
  
  system.time(
    for(i_M in 1:M){
      #print(i_M)
      #set.seed(i_M)
      if((i_M%%(M/20))==0){
        index_print=index_print+1
        print(paste0("finish: ", 5*index_print, " percent"))
        Y_sample_mean_all_record[index_print,]=as.vector(Y_sample_mean_all)
        
      }
      
      log_beta_propose=log_beta_eta_cur[1]+sd_log_beta*rnorm(1);
      log_eta_propose=log_beta_eta_cur[2]+sd_log_eta*rnorm(1); 
      log_beta1_propose=log_beta1_cur+sd_log_beta1*rnorm(1);
      
      log_beta_eta_propose=c(log_beta_propose,log_eta_propose)
      
      
      eigen_R1_propose=Get_eigen_1(R0_1=R0_1,log_beta1=log_beta1_propose,kernel_type=kernel_type,d=d)
      
      
      tilde_eta_vec_propose=exp(log_beta_eta_propose[2])/eigen_R1_propose$values[1:d]
      Y_sample_tilde_propose=t(eigen_R1_propose$vectors[,1:d])%*%(Y_sample-mean_cur)
      
      for(i_row in 1:d){
        
        ###tilde_eta_propose=exp(log_eta_propose)/eigen_R1_cur$values[i_row]
        #  eigen_R1_values_cur[i_row]
        
        log_beta_eta_tilde_propose=c(log_beta_propose,log(tilde_eta_vec_propose[i_row]))
        
        Cov_2_propose[[i_row]]=Construct_Cov_KF(log_beta_eta_tilde_propose,delta_x_2,index_obs,kernel_type)
        
        beta_post_S_2_log_det_propose=Get_beta_post_S_2_log_det(t(Y_sample_tilde_propose[i_row,] ),log_beta_eta_propose,tilde_eta_vec_propose[i_row],Cov_2_propose[[i_row]],add_log_prior = F)
        
        beta_post_propose[i_row,]=beta_post_S_2_log_det_propose
        
      }
      
      
      Trace_Y_sample_tilde_2_propose=sum(Y_sample_tilde_propose^2)
      
      
      
      log_JR_prior_cur=log_approx_ref_prior(log_beta_eta_cur,0,T,CL[2],a,b)+sum(log_beta_eta_cur)
      log_JR_prior_propose=log_approx_ref_prior(log_beta_eta_propose,0,T,CL[2],a,b)+sum(log_beta_eta_propose)
      
      
      log_JR_prior1_cur=log_approx_ref_prior(log_beta1_cur,0,F,CL[1],a,b)
      log_JR_prior1_propose=log_approx_ref_prior(log_beta1_propose,0,F,CL[1],a,b)
      
      log_JR_prior_all_cur=log_JR_prior_cur+log_JR_prior1_cur
      log_JR_prior_all_propose=log_JR_prior_propose+log_JR_prior1_propose
      
      
      
      ##note only 1:d
      
      log_post_all_cur=n2/2*sum(log(tilde_eta_vec_cur))-1/2*sum(beta_post_cur[1:d,3])-(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+ sum(beta_post_cur[1:d,2]) )/(2*sigma_2_0_cur)+log_JR_prior_all_cur+log_beta1_cur+sum(log_beta_eta_cur)
      
      log_post_all_propose=n2/2*sum(log(tilde_eta_vec_propose))-1/2*sum(beta_post_propose[1:d,3])-(Trace_Y_sample_2-Trace_Y_sample_tilde_2_propose+ sum(beta_post_propose[1:d,2]) )/(2*sigma_2_0_cur)+log_JR_prior_all_propose+log_beta1_propose+sum(log_beta_eta_propose)
      
      
      r_ratio=exp( sum(log_post_all_propose)-sum(log_post_all_cur));
      r_ratio
      
      decision=Accept_proposal(r_ratio);
      
      
      param_record[i_M,1:3]=exp(c(log_beta_eta_cur,log_beta1_cur))
      
      
      if(decision){
        log_beta_eta_cur=log_beta_eta_propose
        tilde_eta_vec_cur=tilde_eta_vec_propose
        param_record[i_M,1:2]=exp(log_beta_eta_propose)
        Cov_2_cur=Cov_2_propose
        log_beta1_cur=log_beta1_propose
        param_record[i_M,3]=exp(log_beta1_propose)
        
        eigen_R1_cur=eigen_R1_propose
        
        #U1_cur=U1_propose
        #eigen_R1_values_cur=eigen_R1_values_propose
        tilde_eta_vec_cur=tilde_eta_vec_propose
        
        beta_post_cur=beta_post_propose
        
        log_post_all_cur=log_post_all_propose
        
        Y_sample_tilde_cur=Y_sample_tilde_propose
        
        
        ###need to think why these two are  the same
        Trace_Y_sample_tilde_2_cur=Trace_Y_sample_tilde_2_propose
        
      }
      
      tau_sample=rgamma(1,shape=(n1)*n2/2,rate=(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+ sum(beta_post_cur[1:d,2]) )/2)
      
      sigma_2_0_cur=1/tau_sample
      
      param_record[i_M,4]=sigma_2_0_cur
      
      ##update sigma_2_cur along with the sigma_2_0_cur
      sigma_2_cur=sigma_2_0_cur/exp(log_beta_eta_cur[2])
      
      Y_sample_tilde_cur=t(eigen_R1_cur$vectors[,1:d])%*%(Y_sample-mean_cur)
      
      
      
      need_mean=TRUE;
      need_sample=TRUE;
      need_var=FALSE;
      for(i_row in 1:d){
        ###this needs to be changed
        if(kernel_type=='exp'){
          KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_exp(Cov_2_cur[[i_row]]$C_R_K_Q, Cov_2_cur[[i_row]]$W0,  Cov_2_cur[[i_row]]$GG,  Cov_2_cur[[i_row]]$W,tilde_eta_vec_cur[i_row],
                                                                                index_obs,t(Y_sample_tilde_cur[i_row,] )/sqrt(eigen_R1_cur$values[i_row]*sigma_2_cur),need_sample,need_mean,need_var)
        }else if(kernel_type=='matern_5_2'){
          KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_matern_5_2(Cov_2_cur[[i_row]]$C_R_K_Q, Cov_2_cur[[i_row]]$W0,  Cov_2_cur[[i_row]]$GG,  Cov_2_cur[[i_row]]$W,tilde_eta_vec_cur[i_row],
                                                                                       index_obs,t(Y_sample_tilde_cur[i_row,] )/sqrt(eigen_R1_cur$values[i_row]*sigma_2_cur),need_sample,need_mean,need_var)
        }
        z_hat_matrix[i_row,]=KF_smoother_sample_list[[1]] ##this is even worse
        z_sample_matrix[i_row,]=KF_smoother_sample_list[[2]] ##this is even worse
      }
      
      
      Uz_hat=eigen_R1_cur$vectors[,1:d]%*%((sqrt(eigen_R1_cur$values*sigma_2_cur))*z_hat_matrix)
      
      Y_sample_mean_all=mean_cur+eigen_R1_cur$vectors[,1:d]%*%((sqrt(eigen_R1_cur$values*sigma_2_cur))*z_sample_matrix)
      Y_sample_all=Y_sample_mean_all+sqrt(sigma_2_0_cur )*rnorm(n1*n2)
      
      Y_sample[missing_index]=Y_sample_all[missing_index]
      Y_sample_tilde_cur=t(eigen_R1_cur$vectors[,1:d])%*%(Y_sample-mean_cur) ### t(t(Y_sample)%*%U1_cur)
      
      #image2D(Y_sample_mean_all)
      #image2D(Uz_hat+t(H2%*%B2_cur))
      #image2D(Uz_hat+t(H2%*%B2_hat))
      
      if(i_M>M_0){
        Uz_hat_sum=Uz_hat_sum+Uz_hat
        #Y_sample_sum=Y_sample_sum+Y_sample
        Y_sample_mean_sum=Y_sample_mean_sum+Y_sample_mean_all
        
        
        Y_avg_sample_mean=Y_sample_mean_sum/(i_M-M_0)
        
        record_RMSE_hat_mean[i_M]=sqrt(mean((Y_avg_sample_mean[missing_index]-Y_full[missing_index] )^2 ))
        
        if( (i_M-M_0)<=M_025 ){
          lower_interval_95[,i_M-M_0]=Y_sample[missing_index]
          upper_interval_95[,i_M-M_0]=Y_sample[missing_index]
        }else{
          max_lower=rowMaxs(lower_interval_95)
          min_upper=rowMins(upper_interval_95)
          
          index_lower=cbind(1:length(missing_index),max_lower)
          index_upper=cbind(1:length(missing_index),min_upper)
          
          max_lower_values=lower_interval_95[index_lower]
          min_upper_values=upper_interval_95[index_upper]
          
          
          Y_here=Y_sample[missing_index]
          
          index_fill_lower=which(Y_here<max_lower_values)
          lower_interval_95[index_lower[index_fill_lower,]]=Y_here[index_fill_lower]
          
          index_fill_upper=which(Y_here> min_upper_values)
          
          upper_interval_95[index_upper[index_fill_upper,]]=Y_here[index_fill_upper]
          
          #index  Y_sample[missing_index]
        }
        
        
      }
      
      Y_pred_sum_ori_sample_mean=Y_pred_sum_ori_sample_mean+Y_sample_mean_all
      Y_sample_mean_avg=Y_pred_sum_ori_sample_mean/i_M
      
      record_RMSE_ori_sample_mean[i_M]=sqrt(mean( (Y_sample_mean_avg[missing_index]-Y_full[missing_index])^2))
      record_RMSE_ind_sample_mean[i_M]=sqrt(mean((Y_sample_mean_all[missing_index]-Y_full[missing_index] )^2 ))
      
      sample_difference[i_M]=mean(Y_sample_mean_all[missing_index]-Y_full[missing_index] )
      
      
      #then I need to recompute the likelihood 
      Trace_Y_sample_2=0
      
      for(i_row in 1:n1){
        Trace_Y_sample_2=Trace_Y_sample_2+sum((Y_sample[i_row,]-mean_cur[i_row,])*(Y_sample[i_row,]-mean_cur[i_row,]) )
      }
      
      Trace_Y_sample_tilde_2_cur=sum(Y_sample_tilde_cur^2)
      
      ##update beta_post_cur[1:d]
      
      for(i_d in 1:d){
        
        log_beta_eta=c(log_beta_eta_cur[1],log_beta_eta_cur[2])
        
        #beta_post_S_2_log_det_propose=Get_beta_post_S_2_log_det(t(Y_sample_tilde_cur[i_row,] ),log_beta_eta_propose,tilde_eta_vec_propose[i_row],Cov_2_propose)
        
        beta_post_S_2_log_det=Get_beta_post_S_2_log_det(t(Y_sample_tilde_cur[i_d,] ),log_beta_eta,tilde_eta_vec_cur[i_d],Cov_2_cur[[i_d]],add_log_prior = F)
        beta_post_cur[i_d,]=beta_post_S_2_log_det
      }
      #log_JR_prior=log_approx_ref_prior(log_beta_eta_cur[c(1,2)],0,T,CL[2],a,b)
      
      ##beta_post_cur[,1]= beta_post_cur[,1]+log_JR_prior/d
      
      ##record S_2
      param_record[i_M,1:(2)]=exp(log_beta_eta_cur)     
      param_record[i_M,(3)]=exp(log_beta1_cur)
      
      
      ##one can comment it
      #print(c(record_RMSE_ori[i_M],record_RMSE_ind[i_M]))
      #print( param_record[i_M,(2*d)+2])
      
      #print(c(record_RMSE_ori_sample_mean[i_M], record_RMSE_ind_sample_mean[i_M],sample_difference[i_M],exp(log_beta_eta_cur), exp(log_beta1_cur),sigma_2_0_cur))
      
    }
  )
  
  Uz_avg=Uz_hat_sum/(M-M_0)
  
  Y_pred=Uz_avg
  
  # sqrt(mean((Y_sample[missing_index]-Y_full[missing_index] )^2 ))
  # 
  # sqrt(mean((Y_pred[missing_index]-Y_full[missing_index] )^2 ))
  # sqrt(mean((Y_avg_sample_mean[missing_index]-Y_full[missing_index] )^2 ))
  # sqrt(mean((Y_avg_sample_mean[missing_index]-Y_truth[missing_index] )^2 ))
  # sqrt(mean((Y_pred[missing_index]-Y_truth[missing_index] )^2 ))
  
  
  ###confidence interval 
  lower_interval_95_values=rowMaxs(lower_interval_95,value=T)
  upper_interval_95_values=rowMins(upper_interval_95,value=T)
  # missing_index_record = rep(0,n1*n2)
  # missing_index_record[missing_index] = 1
  # results_fast = cbind(param_record, as.vector(Y_avg_sample_mean),as.vector(Y_full), as.vector(Y_truth),missing_index_record,lower_interval_95_values,upper_interval_95_values)
  # colnames(results_fast) = c('beta',  'eta', 'beta1', 'sigma_2_0', 'Y_avg_sample_mean', 'Y_full','Y_truth','missing_index_record','lower_interval_95_values', 'upper_interval_95_values')
  # 
  filename = paste0("results_parallel/" ,"size_25_",kernel_type,"_", missing_pattern ,"_",missing_per,"/", "Results_fast_seed_",index,".RData")
  # write.csv(results_fast, file = filename,row.names=F)
  
  save(param_record, Y_avg_sample_mean, Y_full, Y_truth, missing_index, lower_interval_95_values, upper_interval_95_values, file = filename)
}

# load(filename)






##########################################################
# parallel computation for full GP method
##########################################################

cl <- makeCluster(8)  # number of clusters                                                                                   
registerDoParallel(cl)

foreach(index = 1:m, .packages = package_names, .export = rfun_names) %dopar% {
  #sourceCpp(file='functions_Dec_2019.cpp')
  ###full Gaussian process
  set.seed(index)
  Y_truth = matrix(Y_truth_all[index,],25,25)
  Y_full = matrix(Y_full_all[index,],25,25)
  Y = matrix(Y_all[index,],25,25)
  n_missing = sum(is.na(Y))
  missing_index = missing_index_all[[index]]

  
  param_record_full_GP=matrix(0,M,4)
  pred_mean=matrix(0,n1,n2)
  
  log_beta_eta_beta1_cur=c(0,0,0)
  sigma_2_0_cur=0.1
  
  input_obs=input[-missing_index,]
  
  Y_obs=Y[-missing_index]
  n_obs = length(Y_obs)
  if(kernel_type=='matern_5_2'){
    R1=matern_5_2_funct(abs(outer(input_obs[,1],input_obs[,1],'-')),exp(log_beta_eta_beta1_cur[3]))
    R2=matern_5_2_funct(abs(outer(input_obs[,2],input_obs[,2],'-')),exp(log_beta_eta_beta1_cur[1]))
  }else if(kernel_type=='exp'){
    R1=pow_exp_funct(abs(outer(input_obs[,1],input_obs[,1],'-')),exp(log_beta_eta_beta1_cur[3]),alpha_i=1)
    R2=pow_exp_funct(abs(outer(input_obs[,2],input_obs[,2],'-')),exp(log_beta_eta_beta1_cur[1]),alpha_i=1)
    
  }
  R_cur=R1*R2
  tilde_R_cur=R_cur+exp(log_beta_eta_beta1_cur[2])*diag(n_obs)
  
  ##
  L_cur=t(chol(tilde_R_cur))
  R_tilde_inv_Y_cur=backsolve(t(L_cur),forwardsolve(L_cur,Y_obs))            
  
  #backsolve(t(L_cur),forwardsolve(L_cur,Y_obs))-solve(tilde_R_cur)%*%Y_obs            
  
  
  
  S_2_cur=t(Y_obs)%*%R_tilde_inv_Y_cur*exp(log_beta_eta_beta1_cur[2])
  log_JR_prior1_cur=log_approx_ref_prior(log_beta_eta_beta1_cur[3],0,F,CL[1],a,b)
  log_JR_prior_cur=log_approx_ref_prior(log_beta_eta_beta1_cur[c(1,2)],0,T,CL[2],a,b)
  
  log_JR_prior_all=log_JR_prior1_cur+log_JR_prior_cur
  log_post_cur=-sum(log(diag(L_cur)))+n_obs/2*log_beta_eta_beta1_cur[2]-S_2_cur/(2*sigma_2_0_cur)+log_JR_prior_all+sum(log_beta_eta_beta1_cur)
  
  
  
  Y_sample_mean_sum_missing=0
  
  
  M_025=0.025*(M-M_0)
  
  lower_interval_95_full_GP=matrix(0,length(missing_index),M_025)
  upper_interval_95_full_GP=matrix(0,length(missing_index),M_025)
  
  r0.list=as.list(1:2)
  r0.list[[1]]=abs(outer(input[missing_index,1],input[-missing_index,1],'-'))
  r0.list[[2]]=abs(outer(input[missing_index,2],input[-missing_index,2],'-'))
  
  r00.list=as.list(1:2)
  r00.list[[1]]=abs(outer(input[missing_index,1],input[missing_index,1],'-'))
  r00.list[[2]]=abs(outer(input[missing_index,2],input[missing_index,2],'-'))
  
  #-n_obs/2*log(sigma_2_0_cur)
  #1/2*determinant(tilde_R_cur)$modulus[1]
  
  Y_sample_mean_full_GP=matrix(0,n1,n2)
  Y_sample_mean_full_GP[-missing_index]=Y[-missing_index]
  
  Y_sample_mean_sum_full_GP=0;
  
  record_RMSE_mean_full_GP_record=rep(NA,M-M_0)
  
  sigma_2_pred=rep(0,n_missing)
  sigma_pred_sum=rep(0,n_missing)
  sigma_2_pred_sum=0
  
  lower_interval_95_full_GP=matrix(0,length(missing_index),M_025)
  upper_interval_95_full_GP=matrix(0,length(missing_index),M_025)
  
  system.time(
    for(i_M in 1:M){
      print(i_M)
      #set.seed(i_M)
      if((i_M%%(M/20))==0){
        index_print=index_print+1
        print(paste0("finish: ", 5*index_print, " percent"))
        #Y_sample_mean_all_record[index_print,]=as.vector(Y_sample_mean_all)
        
      }
      
      log_beta_propose=log_beta_eta_beta1_cur[1]+sd_log_beta*rnorm(1);
      log_eta_propose=log_beta_eta_beta1_cur[2]+sd_log_eta*rnorm(1); 
      log_beta1_propose=log_beta_eta_beta1_cur[3]+sd_log_beta1*rnorm(1); 
      
      
      if(kernel_type=='matern_5_2'){
        R1=matern_5_2_funct(abs(outer(input_obs[,1],input_obs[,1],'-')),exp(log_beta1_propose))
        R2=matern_5_2_funct(abs(outer(input_obs[,2],input_obs[,2],'-')),exp(log_beta_propose))
      }else if(kernel_type=='exp'){
        R1=pow_exp_funct(abs(outer(input_obs[,1],input_obs[,1],'-')),exp(log_beta1_propose),alpha_i=1)
        R2=pow_exp_funct(abs(outer(input_obs[,2],input_obs[,2],'-')),exp(log_beta_propose),alpha_i=1)
        
      }
      R_propose=R1*R2
      
      tilde_R_propose=R_propose+exp(log_eta_propose)*diag(n_obs)
      
      L_propose=t(chol(tilde_R_propose))
      R_tilde_inv_Y_propose=backsolve(t(L_propose),forwardsolve(L_propose,Y_obs))            
      
      #L_propose
      
      S_2_propose=t(Y_obs)%*%R_tilde_inv_Y_propose*exp(log_eta_propose)
      log_JR_prior1_propose=log_approx_ref_prior(log_beta1_propose,0,F,CL[1],a,b)
      log_JR_prior_propose=log_approx_ref_prior(c(log_beta_propose,log_eta_propose),0,T,CL[2],a,b)
      
      log_JR_prior_all=log_JR_prior1_propose+log_JR_prior_propose
      log_post_propose=-sum(log(diag(L_propose)))+n_obs/2*log_eta_propose-S_2_propose/(2*sigma_2_0_cur)+log_JR_prior_all+sum(log_beta1_propose+log_eta_propose+log_beta_propose)
      log_post_propose
      
      
      r_ratio=exp( log_post_propose-log_post_cur);
      r_ratio
      
      #-n_obs/2*log(sigma_2_0_cur)
      
      #sum(beta_post_propose[,1])
      decision=Accept_proposal(r_ratio);
      decision
      param_record_full_GP[i_M,1:3]=exp(log_beta_eta_beta1_cur)
      if(decision){
        log_beta_eta_beta1_cur[1]=log_beta_propose
        log_beta_eta_beta1_cur[2]=log_eta_propose
        log_beta_eta_beta1_cur[3]=log_beta1_propose
        
        param_record_full_GP[i_M,1:3]=exp(log_beta_eta_beta1_cur)
        log_post_cur=log_post_propose
        S_2_cur=S_2_propose
        R_tilde_inv_Y_cur=R_tilde_inv_Y_propose
        L_cur=L_propose
        tilde_R_cur=tilde_R_propose
        R_cur=R_propose
        
      }
      
      ##done for updating beta1
      
      ##update sigma_2 here 
      #tau_sample=rgamma(1,shape=(n1-d)*n2/2,rate=(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+ sum(S_2_cur) )/2)
      tau_sample=rgamma(1,shape=n_obs/2,rate=(S_2_cur)/2)
      
      sigma_2_0_cur=1/tau_sample
      
      param_record_full_GP[i_M,4]=sigma_2_0_cur
      
      ##update sigma_2_cur along with the sigma_2_0_cur
      sigma_2_cur=sigma_2_0_cur/exp(log_beta_eta_beta1_cur[2])
      
      
      #R1=matern_5_2_funct(abs(outer(input1,input1,'-')),beta[1])
      #R2=matern_5_2_funct(abs(outer(input2,input2,'-')),beta[2])
      #R=R1*R2
      ########################
      log_JR_prior1_cur=log_approx_ref_prior(log_beta_eta_beta1_cur[3],0,F,CL[1],a,b)
      log_JR_prior_cur=log_approx_ref_prior(log_beta_eta_beta1_cur[c(1,2)],0,T,CL[2],a,b)
      
      log_JR_prior_all=log_JR_prior1_cur+log_JR_prior_cur
      
      log_post_cur=-sum(log(diag(L_cur)))+n_obs/2*log_beta_eta_beta1_cur[2]-S_2_cur/(2*sigma_2_0_cur)+log_JR_prior_all+sum(log_beta_eta_beta1_cur)
      
      
      #  -n_obs/2*log(sigma_2_0_cur)
      
      #R1=matern_5_2_funct(abs(outer(input1,input1,'-')),log_beta_eta_beta1_cur[3])
      #R2=matern_5_2_funct(abs(outer(input2,input2,'-')),log_beta_eta_beta1_cur[1])
      #R_all=outer(R1,R2)
      ##It took some time
      if(kernel_type=='matern_5_2'){
        r=separable_kernel(r0.list,exp(log_beta_eta_beta1_cur[c(3,1)]),kernel_type,alpha=c(1,1))
        R_star=separable_kernel(r00.list,exp(log_beta_eta_beta1_cur[c(3,1)]),kernel_type,alpha=c(1,1))
      }else if(kernel_type=='exp'){
        r=separable_kernel(r0.list,exp(log_beta_eta_beta1_cur[c(3,1)]),'pow_exp',alpha=c(1,1))
        R_star=separable_kernel(r00.list,exp(log_beta_eta_beta1_cur[c(3,1)]),'pow_exp',alpha=c(1,1))
      }
      
      ##only look at missing one
      R_tilde_inv_r=backsolve(t(L_cur),forwardsolve(L_cur,t(r)))        
      
      ##here I only compute the missing thing
      Y_sample_mean_hat_full_GP=r%*%R_tilde_inv_Y_cur
      
      pred_cov_missing=(sigma_2_cur+sigma_2_0_cur)*R_star- sigma_2_cur*r%*%R_tilde_inv_r
      
      L_pred=t(chol(pred_cov_missing))
      
      
      Y_sample_missing=Y_sample_mean_hat_full_GP+L_pred%*%rnorm(n_missing)
      
      #system.time(
      # for(i_missing in 1:n_missing){
      #   sigma_2_pred[i_missing]=sigma_2_cur+sigma_2_0_cur- sigma_2_cur*r[i_missing,]%*%R_tilde_inv_r[,i_missing]
      # }
      #)
      # system.time(
      #   for(i_missing in 1:n_missing){
      #     R_tilde_inv_r=backsolve(t(L_cur),forwardsolve(L_cur,(r[i_missing,])))            
      #     
      #     sigma_2_pred[i_missing]=sigma_2_cur+sigma_2_0_cur- sigma_2_cur*r[i_missing,]%*%R_tilde_inv_r
      #   }
      # )
      
      
      ##this is only from the marginal
      # Y_sample_missing_marginal=Y_sample_mean_hat_full_GP[missing_index]+sqrt(sigma_2_pred)*rnorm(n_missing)
      
      if(i_M>M_0){
        #Uz_hat_sum=Uz_hat_sum+Uz_hat
        #Y_sample_sum=Y_sample_sum+Y_sample
        Y_sample_mean_sum_full_GP=Y_sample_mean_sum_full_GP+Y_sample_mean_hat_full_GP
        
        # sigma_2_pred_sum=sigma_2_pred_sum+sigma_2_pred
        # sigma_pred_sum=sigma_pred_sum+sqrt(sigma_2_pred)
        
        if( (i_M-M_0)<=M_025 ){
          lower_interval_95_full_GP[,i_M-M_0]=Y_sample_missing
          upper_interval_95_full_GP[,i_M-M_0]=Y_sample_missing
        }else{
          max_lower=rowMaxs(lower_interval_95_full_GP)
          min_upper=rowMins(upper_interval_95_full_GP)
          
          index_lower=cbind(1:length(missing_index),max_lower)
          index_upper=cbind(1:length(missing_index),min_upper)
          
          max_lower_values=lower_interval_95_full_GP[index_lower]
          min_upper_values=upper_interval_95_full_GP[index_upper]
          
          
          #Y_here=Y_sample[missing_index]
          
          Y_here=Y_sample_missing
          
          index_fill_lower=which(Y_here<max_lower_values)
          lower_interval_95_full_GP[index_lower[index_fill_lower,]]=Y_here[index_fill_lower]
          
          index_fill_upper=which(Y_here> min_upper_values)
          
          upper_interval_95_full_GP[index_upper[index_fill_upper,]]=Y_here[index_fill_upper]
          
          #index  Y_sample[missing_index]
        }
        
      }
      
      
      
      print(c(exp(log_beta_eta_beta1_cur),sigma_2_0_cur))
      
      
    }
  )
  
  
  Y_sample_mean_avg_full_GP_missing=Y_sample_mean_sum_full_GP/(M-M_0)
  Y_sample_mean_avg_full_GP=matrix(0,n1,n2)
  Y_sample_mean_avg_full_GP[missing_index]=Y_sample_mean_avg_full_GP_missing
  Y_sample_mean_avg_full_GP[-missing_index]=Y[-missing_index]
  
  Y_sample_mean_avg_full_GP=matrix(Y_sample_mean_avg_full_GP,n1,n2)
  #image2D(Y_sample_mean_avg_full_GP)
  
  # sqrt(mean( (Y_sample_mean_avg_full_GP[missing_index]- Y_full[missing_index])^2))
  # sqrt(mean( (Y_sample_mean_avg_full_GP[missing_index]- Y_truth[missing_index])^2))
  
  ##these two are  close
  # sqrt((mean( abs(Y_sample_mean_avg_full_GP[missing_index]-Y_avg_sample_mean[missing_index]))))
  
  # sigma_2_pred_avg=sigma_2_pred_sum/(M-M_0)
  # sigma_pred_avg=sigma_pred_sum/(M-M_0)
  # 
  # LB=Y_sample_mean_avg_full_GP[missing_index]+sqrt(sigma_2_pred_avg)*qnorm(0.025)
  # UB=Y_sample_mean_avg_full_GP[missing_index]+sqrt(sigma_2_pred_avg)*qnorm(0.975)
  
  lower_interval_95_values_full_GP=rowMaxs(lower_interval_95_full_GP,value=T)
  upper_interval_95_values_full_GP=rowMins(upper_interval_95_full_GP,value=T)
  # 
  # sum( Y_full[missing_index]<upper_interval_95_values_full_GP& Y_full[missing_index]>lower_interval_95_values_full_GP)/length(upper_interval_95_values_full_GP)
  # 
  # sum( Y_full[missing_index]<UB& Y_full[missing_index]>LB)/length(upper_interval_95_values_full_GP)
  
  #upper_interval_95_values_full_GP[missing_index]-Y_sample_mean_avg_full_GP[missing_index]
  
  #upper_interval_95_full_GP-lower_interval_95_full_GP
  
  #sum( Y_full[missing_index]<UB& Y_full[missing_index]>LB)/length(LB)
  
  ###difference of the interval. These intervals are for missing so 
  # mean( abs(upper_interval_95_values_full_GP-upper_interval_95_values))
  # mean( abs(lower_interval_95_values_full_GP-lower_interval_95_values))
  
  #sd(upper_interval_95_values)
  
  #sd(lower_interval_95_values)
  
  
  # sum(LB-lower_interval_95_values>0)/length(lower_interval_95_values)
  # 
  # sum(UB- upper_interval_95_values<0)/length(lower_interval_95_values)
  # 
  # mean(UB-LB)
  
  # mean((param_record_full_GP[M_0:M,2])/param_record_full_GP[M_0:M,4])
  # 
  # mean(param_record_full_GP[M_0:M,2])
  
  # max_value=max(Y_full,Y,Y_avg_sample_mean,na.rm=T)
  # min_value=min(Y_full,Y,Y_avg_sample_mean,na.rm=T)
  # Y_avg_sample_mean_mat=matrix(Y_avg_sample_mean,n1,n2)
  # missing_index_record = rep(0,n1*n2)
  # missing_index_record[missing_index] = 1
  # results_full_GP = cbind(param_record_full_GP, as.vector(Y_sample_mean_avg_full_GP),as.vector(Y_full), as.vector(Y_truth),missing_index_record,lower_interval_95_values_full_GP,upper_interval_95_values_full_GP)
  # colnames(results_full_GP) = c('beta',  'eta', 'beta1', 'sigma_2_0', 'Y_sample_mean_avg_full_GP', 'Y_full','Y_truth','missing_index_record','lower_interval_95_values_full_GP', 'upper_interval_95_values_full_GP')
  # 
  filename = paste0("./result/Parallel/" ,"size_25_",kernel_type,"_", missing_pattern ,"_",missing_per,"/", "Results_full_GP_seed_",index,".RData")
  # write.csv(results_full_GP, file = filename,row.names=F)
  
  save(param_record_full_GP, Y_sample_mean_avg_full_GP, Y_full, Y_truth, missing_index, lower_interval_95_values_full_GP, upper_interval_95_values_full_GP, file = filename)
}
stopCluster(cl)  
