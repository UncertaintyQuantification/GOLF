#library(ggplot2)
library(plot3D)
#load("Data_Heaton/AllSatelliteTemps.RData")
library(RobustGaSP)
#library(FastGaSP)
library(Rcpp)
library(RcppEigen)
 sourceCpp(file='../../src/functions_Feb_26_2020.cpp')

 seed_here=0
set.seed(seed_here)

disk_missing = function(center, radius){
  m = matrix(0,n1,n2)
  g = expand.grid(1:nrow(m), 1:nrow(m))
  g$d2 = sqrt((g$Var1-center[1])^2 + (g$Var2-center[2])^2)
  g$inside = g$d2<=radius
  return(which(g$inside==T))
}

kernel = 'matern_5_2' #c('matern_5_2', 'exp')
missing_per = 0.2 #c(0.5,0.2)
missing_pattern = "disk" #c("random","disk")



# set.seed(1)


items = data.frame(kernel, missing_per, missing_pattern)

###and show filtering can be done using two-dim filtering 
num_simulation=100

n1=100
n2=100
N=n1*n2
d_real=n1*0.3  ## d_real = 30
if (items$missing_pattern=="random"){
  missing_index=sample(1:(n1*n2),size=items$missing_per*n1*n2)
} else if (items$missing_pattern=="disk"){
  if (items$missing_per==0.5){
    radius = 40
    missing_index = disk_missing(c(n1/2,n2/2),radius)
  }
  if (items$missing_per==0.2){
    radius = 25
    missing_index = disk_missing(c(n1/2,n2/2),radius)
  }
}
# length(missing_index)/N

separable_truth_data_all = matrix(0,100,N) # each row is a simulated data set
nonseparable_truth_data_all = matrix(0,100,N) # each row is a simulated data set

separable_full_data_all = matrix(0,100,N) # each row is a simulated data set
nonseparable_full_data_all = matrix(0,100,N) # each row is a simulated data set

# separable_sample = F

##let me do equally space first
input1=seq(0,1,1/(n1-1)) ##row
input2=seq(0,1,1/(n2-1)) ##column

input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)

input=cbind(rep(input1,n2),input_2_vec)

##sample the truth
#beta=c(2,5)
#kernel_type='exp'

kernel_type=as.character(items$kernel)



for (data_index in 1:num_simulation){
  set.seed(data_index)
  if(kernel_type=='matern_5_2'){
    beta=c(4,2)
    R1=matern_5_2_funct(abs(outer(input1,input1,'-')),beta[1])
    R2=matern_5_2_funct(abs(outer(input2,input2,'-')),beta[2])
  }else if(kernel_type=='exp'){
    beta=c(1/4,1)
    # beta=c(1/10,1)
    R1=pow_exp_funct(abs(outer(input1,input1,'-')),beta[1],alpha_i=1)
    R2=pow_exp_funct(abs(outer(input2,input2,'-')),beta[2],alpha_i=1)
    
  }
  
  # R1=matern_5_2_funct(abs(outer(input1,input1,'-')),beta[1])
  # R2=matern_5_2_funct(abs(outer(input2,input2,'-')),beta[2])
  
  ##cholesky
  eigen_R1=eigen(R1)
  L1=eigen_R1$vector[,1:d_real]%*%diag(sqrt(eigen_R1$values[1:d_real]))
  L2=t(chol(R2))
  
  Z=matrix(rnorm(d_real*n2),d_real,n2)
  
  Z2=t(L2%*%t(Z)) 
  eta=0.1
  
  Y_truth=L1%*%Z2
  separable_truth_data_all[data_index,] = as.vector(Y_truth)
  Y_full=L1%*%Z2+eta*matrix(rnorm(n1*n2),n1,n2)
  separable_full_data_all[data_index,] = as.vector(Y_full)
}


for (data_index in 1:num_simulation){
  set.seed(data_index)
  

  if(kernel_type=='matern_5_2'){
    beta=3
    #beta=1
    
    R1=matern_5_2_funct(abs(outer(input1,input1,'-')),beta[1])
  }else if(kernel_type=='exp'){
    beta=c(1/4)
    R1=pow_exp_funct(abs(outer(input1,input1,'-')),beta[1],alpha_i=1)
  }
  R1_eigen=eigen(R1)
  L1=R1_eigen$vector[,1:d_real]%*%diag( sqrt(R1_eigen$values[1:d_real]))
  
  # seed_here=1
  # set.seed(seed_here)
  ##case 1. nonseparable kernel with a mean function  with add mean
  
  Z=matrix(rnorm(d_real*n2),d_real,n2)
  Z2=matrix(0,d_real,n2)
  for(i in 1:d_real){
    if(kernel_type=='matern_5_2'){
      beta=(i)*1
      R2=matern_5_2_funct(abs(outer(input2,input2,'-')),beta)
    }else if(kernel_type=='exp'){
      beta=(i/4)
      R2=pow_exp_funct(abs(outer(input2,input2,'-')),beta,alpha_i=1)
    }
    L2=t(chol(R2))
    Z2[i,]=t(L2%*%(Z[i,]))
  }
  eta=0.1
  
  ##add mean
  #B2_real=matrix(seq(0,3,3/(n1-1)),1,n1)
  H2=matrix(1,n2,1)
  ##Y_truth=L1%*%Z2+t(H2%*%B2_real)   ##have mean
  ##for simplicity, let's just do zero mean
  Y_truth = L1%*%Z2
  nonseparable_truth_data_all[data_index, ]=as.vector(Y_truth)                  ##have no mean
  
  nonseparable_full_data_all[data_index,]=as.vector(Y_truth+eta*matrix(rnorm(n1*n2),n1,n2))
  
}

#save(missing_index, separable_truth_data_all, separable_full_data_all, file="All Separable Data.RData")

#save(missing_index, nonseparable_truth_data_all, nonseparable_full_data_all, file="All Nonseparable Data.RData")


# 
image2D(matrix(nonseparable_full_data_all[5,],n1,n2))
image2D(matrix(nonseparable_truth_data_all[5,],n1,n2))
# 
image2D(matrix(separable_full_data_all[2,],n1,n2))
image2D(matrix(separable_truth_data_all[2,],n1,n2))


disk_missing = function(center, radius){
  #m = matrix(0,n1,n2)
  #g = expand.grid(1:nrow(m), 1:nrow(m))
  #m = matrix(0,n1,n2)
  g = expand.grid(1:n1, 1:n2)
  
  g$d2 = sqrt((g$Var1-center[1])^2 + (g$Var2-center[2])^2)
  g$inside = g$d2<=radius
  return(which(g$inside==T))
}

missing_per=0.2
if (missing_per==0.2){
  radius = 25
  missing_index = disk_missing(c(n1/2,n2/2),radius)
}
length(missing_index)/(n1*n2)

###compare separable with nonseparable when the truth data is separable 
d=d_real

return_list_nonseparable_record=as.list(1:num_simulation)
return_list_separable_record=as.list(1:num_simulation)

record_nonseparable_error=matrix(0,num_simulation,4)
record_separable_error=matrix(0,num_simulation,4)

#record_nonseparable_coverage=rep(0,M)
#record_separable_coverage=rep(0,M)
#record_nonseparable_length=rep(0,M)
#record_separable_length=rep(0,M)

record_simulation_time=system.time(
  for(i_M in 1:num_simulation){
    ##make data
    print(i_M)
    Y_full=as.matrix(separable_full_data_all[i_M,],n1,n2)
    Y_truth=as.matrix(separable_truth_data_all[i_M,],n1,n2)
    Y=matrix(NA,n1,n2)
    Y[-missing_index]=Y_full[-missing_index]
    #image2D((Y),input1,input2)
    missing_index_all=missing_index
    
    input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)
    
    input=cbind(rep(input1,n2),input_2_vec)
    
    ##number of d
    #d=n1
    ##try to change it to explore with less d and more d
    #d=d_real
    #d=10
    
    
    ##default jointly robust prior
    a=1/2-2
    b=1
    
    CL = rep(0,2)    ###CL is also used in the prior so I make it a model parameter
    
    CL[1] = (max(input[,1])-min(input[,1]))/n1
    CL[2] = (max(input[,2])-min(input[,2]))/n2
    
    ##prior parameter CL[1:2], a, b
    prior_par=c(CL,a,b)
    
    sd_log_beta=sd_log_eta=40/n2
    sd_log_beta1=40/n1
    
    step_size_kernel_par=c(sd_log_beta,sd_log_eta,sd_log_beta1)
    
    log_beta1_cur=3
    log_beta_eta_cur=c(rep(0,d),rep(-3,d))  ##this is for the second one
    initial_kernel_par=c(log_beta_eta_cur,log_beta1_cur)
    
    delta_x_2=input2[2:n2]-input2[1:(n2-1)]
    
    R0_1=abs(outer(input1,input1,'-'))
    
    ##already defined before
    #kernel_type='exp'
    #kernel_type='matern_5_2'
    ##okay reimplement this
    missing_index_all_mat=matrix(0,length(missing_index_all),2)
    missing_index_all_mat[,1]=as.integer(missing_index_all%%n1-1)
    
    missing_index_all_mat[,2]=as.integer(floor(missing_index_all/n1))
    
    Y_sample=Y
    
    ##use column mean to initiate the initial value
    index_missing_col=list()
    #
    for(i in 1:n1){
      index_missing_col[[i]]=which(is.na(Y[i,])==T)
      Y_sample[i, index_missing_col[[i]]]=mean(Y[i,],na.rm=T)
    }
    n_missing=length(missing_index)
    
    ##
    for(i in 1:n1){
      j=i;
      while(sum(is.na(Y_sample[i,]))>0){
        j=j+1
        Y_sample[i, ]=mean(Y[j,],na.rm=T)
        
      }
    }
    
    #image2D(Y_sample)
    
    Y_sample_1=Y_sample
    
    
    
    M=5000
    M_0=M/5
    interval_prop=0.025
    need_interval=T
    
    H2=matrix(1,n2,1)
    
    ##typically we will suggest user to add trend if the found it from the training data
    #plot(rowMeans(Y,na.rm=T))  #this indicate we should add a trend
    
    #have mean or not
    have_mean=F
    
    #system.time(
    #for(ii in 1:1){
    return_list_nonseparable=GOLF_2D_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                                          M, M_0, d, missing_index, missing_index_all_mat,
                                          prior_par,initial_kernel_par,step_size_kernel_par,seed_here,
                                          need_interval,interval_prop,have_mean, H2)
    #}
    #)
    return_list_nonseparable_record[[i_M]]=return_list_nonseparable
    
    record_nonseparable_error[i_M,1]=sqrt(mean(((return_list_nonseparable[[1]][missing_index])-Y_full[missing_index])^2))
    record_nonseparable_error[i_M,2]=sqrt(mean(((return_list_nonseparable[[1]][missing_index])-Y_truth[missing_index])^2))
    ##coverage
    record_nonseparable_error[i_M,3]=length(which(return_list_nonseparable[[4]]<=Y_full[missing_index] & return_list_nonseparable[[5]]>=Y_full[missing_index]))/length(missing_index)
    ##length of interval
    record_nonseparable_error[i_M,4]=mean(return_list_nonseparable[[5]]-return_list_nonseparable[[4]])
    
    
    ###change some of the parameter
    ###need to change the step size when correlation are shared. We need a smaller step size
    sd_log_beta_separable=.05  ##this may depend on the number of n
    sd_log_beta1_separable=.05
    sd_log_eta_separable=0.2
    step_size_kernel_par_separable=(c(sd_log_beta_separable,sd_log_eta_separable,sd_log_beta1_separable))
    ###the number of initial parameters should also be changed
    log_beta1_separable_cur=3
    log_beta_eta_separable_cur=c(0,-3)  ##this is for the second one
    initial_kernel_par_separable=(c(log_beta_eta_separable_cur,log_beta1_separable_cur))
    
    
    
    #system.time(
    #for(ii in 1:1){
    return_list_separable=GOLF_2D_separable_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                                                 M, M_0, d, missing_index, missing_index_all_mat,
                                                 prior_par,initial_kernel_par_separable,step_size_kernel_par_separable,seed_here,
                                                 need_interval,interval_prop,have_mean, H2)
    #}
    #)
    
    return_list_separable_record[[i_M]]=return_list_separable
    
    record_separable_error[i_M,1]=sqrt(mean(((return_list_separable[[1]][missing_index])-Y_full[missing_index])^2))
    record_separable_error[i_M,2]=sqrt(mean(((return_list_separable[[1]][missing_index])-Y_truth[missing_index])^2))
    ##coverage
    record_separable_error[i_M,3]=length(which(return_list_separable[[4]]<=Y_full[missing_index] & return_list_separable[[5]]>=Y_full[missing_index]))/length(missing_index)
    ##length of interval
    record_separable_error[i_M,4]=mean(return_list_separable[[5]]-return_list_separable[[4]])
    
    print(c(record_nonseparable_error[i_M,],record_separable_error[i_M,]))
    
  }
)

#colMeans(record_nonseparable_error[1:72,])
#colMeans(record_separable_error[1:72,])


save(record_nonseparable_error, record_separable_error, return_list_separable_record,return_list_nonseparable_record, file="separable_prediction_d_real_30_d_30.Rdata")

colMeans(record_nonseparable_error)
colMeans(record_separable_error)

#####
d=10

return_list_nonseparable_record=as.list(1:num_simulation)
return_list_separable_record=as.list(1:num_simulation)

record_nonseparable_error=matrix(0,num_simulation,4)
record_separable_error=matrix(0,num_simulation,4)

#record_nonseparable_coverage=rep(0,M)
#record_separable_coverage=rep(0,M)
#record_nonseparable_length=rep(0,M)
#record_separable_length=rep(0,M)

record_simulation_time=system.time(
for(i_M in 1:num_simulation){
  ##make data
  print(i_M)
  Y_full=as.matrix(separable_full_data_all[i_M,],n1,n2)
  Y_truth=as.matrix(separable_truth_data_all[i_M,],n1,n2)
  Y=matrix(NA,n1,n2)
  Y[-missing_index]=Y_full[-missing_index]
  #image2D((Y),input1,input2)
  missing_index_all=missing_index
  
  input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)
  
  input=cbind(rep(input1,n2),input_2_vec)
  
  ##number of d
  #d=n1
  ##try to change it to explore with less d and more d
  #d=d_real
  #d=10
  
  
  ##default jointly robust prior
  a=1/2-2
  b=1
  
  CL = rep(0,2)    ###CL is also used in the prior so I make it a model parameter
  
  CL[1] = (max(input[,1])-min(input[,1]))/n1
  CL[2] = (max(input[,2])-min(input[,2]))/n2
  
  ##prior parameter CL[1:2], a, b
  prior_par=c(CL,a,b)
  
  sd_log_beta=sd_log_eta=40/n2
  sd_log_beta1=40/n1
  
  step_size_kernel_par=c(sd_log_beta,sd_log_eta,sd_log_beta1)
  
  log_beta1_cur=3
  log_beta_eta_cur=c(rep(0,d),rep(-3,d))  ##this is for the second one
  initial_kernel_par=c(log_beta_eta_cur,log_beta1_cur)
  
  delta_x_2=input2[2:n2]-input2[1:(n2-1)]
  
  R0_1=abs(outer(input1,input1,'-'))
  
  ##already defined before
  #kernel_type='exp'
  #kernel_type='matern_5_2'
  ##okay reimplement this
  missing_index_all_mat=matrix(0,length(missing_index_all),2)
  missing_index_all_mat[,1]=as.integer(missing_index_all%%n1-1)
  
  missing_index_all_mat[,2]=as.integer(floor(missing_index_all/n1))
  
  Y_sample=Y
  
  ##use column mean to initiate the initial value
  index_missing_col=list()
  #
  for(i in 1:n1){
    index_missing_col[[i]]=which(is.na(Y[i,])==T)
    Y_sample[i, index_missing_col[[i]]]=mean(Y[i,],na.rm=T)
  }
  n_missing=length(missing_index)
  
  ##
  for(i in 1:n1){
    j=i;
    while(sum(is.na(Y_sample[i,]))>0){
      j=j+1
      Y_sample[i, ]=mean(Y[j,],na.rm=T)
      
    }
  }
  
  #image2D(Y_sample)
  
  Y_sample_1=Y_sample
  
  
  
  M=5000
  M_0=M/5
  interval_prop=0.025
  need_interval=T
  
  H2=matrix(1,n2,1)
  
  ##typically we will suggest user to add trend if the found it from the training data
  #plot(rowMeans(Y,na.rm=T))  #this indicate we should add a trend
  
  #have mean or not
  have_mean=F
  
  #system.time(
    #for(ii in 1:1){
      return_list_nonseparable=GOLF_2D_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                               M, M_0, d, missing_index, missing_index_all_mat,
                               prior_par,initial_kernel_par,step_size_kernel_par,seed_here,
                               need_interval,interval_prop,have_mean, H2)
    #}
  #)
  return_list_nonseparable_record[[i_M]]=return_list_nonseparable
  
  record_nonseparable_error[i_M,1]=sqrt(mean(((return_list_nonseparable[[1]][missing_index])-Y_full[missing_index])^2))
  record_nonseparable_error[i_M,2]=sqrt(mean(((return_list_nonseparable[[1]][missing_index])-Y_truth[missing_index])^2))
  ##coverage
  record_nonseparable_error[i_M,3]=length(which(return_list_nonseparable[[4]]<=Y_full[missing_index] & return_list_nonseparable[[5]]>=Y_full[missing_index]))/length(missing_index)
  ##length of interval
  record_nonseparable_error[i_M,4]=mean(return_list_nonseparable[[5]]-return_list_nonseparable[[4]])
  
  
  ###change some of the parameter
  ###need to change the step size when correlation are shared. We need a smaller step size
  sd_log_beta_separable=.05  ##this may depend on the number of n
  sd_log_beta1_separable=.05
  sd_log_eta_separable=0.2
  step_size_kernel_par_separable=(c(sd_log_beta_separable,sd_log_eta_separable,sd_log_beta1_separable))
  ###the number of initial parameters should also be changed
  log_beta1_separable_cur=3
  log_beta_eta_separable_cur=c(0,-3)  ##this is for the second one
  initial_kernel_par_separable=(c(log_beta_eta_separable_cur,log_beta1_separable_cur))
  
  
  
  #system.time(
    #for(ii in 1:1){
      return_list_separable=GOLF_2D_separable_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                                                   M, M_0, d, missing_index, missing_index_all_mat,
                                                   prior_par,initial_kernel_par_separable,step_size_kernel_par_separable,seed_here,
                                                   need_interval,interval_prop,have_mean, H2)
    #}
  #)
  
  return_list_separable_record[[i_M]]=return_list_separable
  
  record_separable_error[i_M,1]=sqrt(mean(((return_list_separable[[1]][missing_index])-Y_full[missing_index])^2))
  record_separable_error[i_M,2]=sqrt(mean(((return_list_separable[[1]][missing_index])-Y_truth[missing_index])^2))
  ##coverage
  record_separable_error[i_M,3]=length(which(return_list_separable[[4]]<=Y_full[missing_index] & return_list_separable[[5]]>=Y_full[missing_index]))/length(missing_index)
  ##length of interval
  record_separable_error[i_M,4]=mean(return_list_separable[[5]]-return_list_separable[[4]])
  
  print(c(record_nonseparable_error[i_M,],record_separable_error[i_M,]))
  
}
)

#sqrt(mean(record_nonseparable_error[,1]^2))
#sqrt(mean(record_separable_error[,1]^2))


save(record_nonseparable_error, record_separable_error, return_list_separable_record,return_list_nonseparable_record, file="separable_prediction_d_real_30_d_10.Rdata")


###compare separable with nonseparable when the truth data is separable 

d=20

return_list_nonseparable_record=as.list(1:num_simulation)
return_list_separable_record=as.list(1:num_simulation)

record_nonseparable_error=matrix(0,num_simulation,4)
record_separable_error=matrix(0,num_simulation,4)

#record_nonseparable_coverage=rep(0,M)
#record_separable_coverage=rep(0,M)
#record_nonseparable_length=rep(0,M)
#record_separable_length=rep(0,M)

record_simulation_time=system.time(
  for(i_M in 1:num_simulation){
    ##make data
    print(i_M)
    Y_full=as.matrix(separable_full_data_all[i_M,],n1,n2)
    Y_truth=as.matrix(separable_truth_data_all[i_M,],n1,n2)
    Y=matrix(NA,n1,n2)
    Y[-missing_index]=Y_full[-missing_index]
    #image2D((Y),input1,input2)
    missing_index_all=missing_index
    
    input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)
    
    input=cbind(rep(input1,n2),input_2_vec)
    
    ##number of d
    #d=n1
    ##try to change it to explore with less d and more d
    #d=d_real
    #d=10
    
    
    ##default jointly robust prior
    a=1/2-2
    b=1
    
    CL = rep(0,2)    ###CL is also used in the prior so I make it a model parameter
    
    CL[1] = (max(input[,1])-min(input[,1]))/n1
    CL[2] = (max(input[,2])-min(input[,2]))/n2
    
    ##prior parameter CL[1:2], a, b
    prior_par=c(CL,a,b)
    
    sd_log_beta=sd_log_eta=40/n2
    sd_log_beta1=40/n1
    
    step_size_kernel_par=c(sd_log_beta,sd_log_eta,sd_log_beta1)
    
    log_beta1_cur=3
    log_beta_eta_cur=c(rep(0,d),rep(-3,d))  ##this is for the second one
    initial_kernel_par=c(log_beta_eta_cur,log_beta1_cur)
    
    delta_x_2=input2[2:n2]-input2[1:(n2-1)]
    
    R0_1=abs(outer(input1,input1,'-'))
    
    ##already defined before
    #kernel_type='exp'
    #kernel_type='matern_5_2'
    ##okay reimplement this
    missing_index_all_mat=matrix(0,length(missing_index_all),2)
    missing_index_all_mat[,1]=as.integer(missing_index_all%%n1-1)
    
    missing_index_all_mat[,2]=as.integer(floor(missing_index_all/n1))
    
    Y_sample=Y
    
    ##use column mean to initiate the initial value
    index_missing_col=list()
    #
    for(i in 1:n1){
      index_missing_col[[i]]=which(is.na(Y[i,])==T)
      Y_sample[i, index_missing_col[[i]]]=mean(Y[i,],na.rm=T)
    }
    n_missing=length(missing_index)
    
    ##
    for(i in 1:n1){
      j=i;
      while(sum(is.na(Y_sample[i,]))>0){
        j=j+1
        Y_sample[i, ]=mean(Y[j,],na.rm=T)
        
      }
    }
    
    #image2D(Y_sample)
    
    Y_sample_1=Y_sample
    
    
    
    M=5000
    M_0=M/5
    interval_prop=0.025
    need_interval=T
    
    H2=matrix(1,n2,1)
    
    ##typically we will suggest user to add trend if the found it from the training data
    #plot(rowMeans(Y,na.rm=T))  #this indicate we should add a trend
    
    #have mean or not
    have_mean=F
    
    #system.time(
    #for(ii in 1:1){
    return_list_nonseparable=GOLF_2D_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                                          M, M_0, d, missing_index, missing_index_all_mat,
                                          prior_par,initial_kernel_par,step_size_kernel_par,seed_here,
                                          need_interval,interval_prop,have_mean, H2)
    #}
    #)
    return_list_nonseparable_record[[i_M]]=return_list_nonseparable
    
    record_nonseparable_error[i_M,1]=sqrt(mean(((return_list_nonseparable[[1]][missing_index])-Y_full[missing_index])^2))
    record_nonseparable_error[i_M,2]=sqrt(mean(((return_list_nonseparable[[1]][missing_index])-Y_truth[missing_index])^2))
    ##coverage
    record_nonseparable_error[i_M,3]=length(which(return_list_nonseparable[[4]]<=Y_full[missing_index] & return_list_nonseparable[[5]]>=Y_full[missing_index]))/length(missing_index)
    ##length of interval
    record_nonseparable_error[i_M,4]=mean(return_list_nonseparable[[5]]-return_list_nonseparable[[4]])
    
    
    ###change some of the parameter
    ###need to change the step size when correlation are shared. We need a smaller step size
    sd_log_beta_separable=.05  ##this may depend on the number of n
    sd_log_beta1_separable=.05
    sd_log_eta_separable=0.2
    step_size_kernel_par_separable=(c(sd_log_beta_separable,sd_log_eta_separable,sd_log_beta1_separable))
    ###the number of initial parameters should also be changed
    log_beta1_separable_cur=3
    log_beta_eta_separable_cur=c(0,-3)  ##this is for the second one
    initial_kernel_par_separable=(c(log_beta_eta_separable_cur,log_beta1_separable_cur))
    
    
    
    #system.time(
    #for(ii in 1:1){
    return_list_separable=GOLF_2D_separable_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                                                 M, M_0, d, missing_index, missing_index_all_mat,
                                                 prior_par,initial_kernel_par_separable,step_size_kernel_par_separable,seed_here,
                                                 need_interval,interval_prop,have_mean, H2)
    #}
    #)
    
    return_list_separable_record[[i_M]]=return_list_separable
    
    record_separable_error[i_M,1]=sqrt(mean(((return_list_separable[[1]][missing_index])-Y_full[missing_index])^2))
    record_separable_error[i_M,2]=sqrt(mean(((return_list_separable[[1]][missing_index])-Y_truth[missing_index])^2))
    ##coverage
    record_separable_error[i_M,3]=length(which(return_list_separable[[4]]<=Y_full[missing_index] & return_list_separable[[5]]>=Y_full[missing_index]))/length(missing_index)
    ##length of interval
    record_separable_error[i_M,4]=mean(return_list_separable[[5]]-return_list_separable[[4]])
    
    print(c(record_nonseparable_error[i_M,],record_separable_error[i_M,]))
    
  }
)

#colMeans(record_nonseparable_error[1:72,])
#colMeans(record_separable_error[1:72,])


save(record_nonseparable_error, record_separable_error, return_list_separable_record,return_list_nonseparable_record, file="separable_prediction_d_real_30_d_20.Rdata")


