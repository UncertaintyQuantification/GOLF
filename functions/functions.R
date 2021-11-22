# Intial_Y_sample_KF<-function(Y,grid_with_full_obs_Y,log_beta_eta,delta_x_2){
#   beta=exp(log_beta_eta[1])
#   eta=exp(log_beta_eta[2])
#   
#   lambda=1.0*beta
#   
#   GG_cur=Construct_G_exp( delta_x_2,  lambda)
#   W0_cur=Construct_W0_exp(1, lambda)
#   W_cur=Construct_W_exp(1, delta_x_2, lambda, W0_cur) 
#   
# 
#   Y_sample_return=Y
#     
#   for(i_row in 1:n1){
#     index_missing_row=which(is.na(Y[i_row,])==T)
#     index_obs=rep(1,n2)
#     index_obs[index_missing_row]=0
#     #lambda2=sqrt(5.0)*exp(log_beta_eta_cur[i_row])
#     ##matern 2.5
#     #lambda2=sqrt(5.0)*exp(log_beta_propose)
#     ##exp
# 
# 
#     C_R_K_Q_cur=Get_C_R_K_Q_pred(index_obs, GG_cur,W_cur, W0_cur,   eta );
#     
#     if(length(index_missing_row)>0){
#       Y_here=as.vector(Y[i_row,-index_missing_row])
#     }else{
#       Y_here=as.vector(Y[i_row,])
#     }
#     
#     need_mean=TRUE;
#     need_sample=FALSE;
#     need_var=FALSE;
#     ###this needs to be changed
#     KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_exp(C_R_K_Q_cur, W0_cur, GG_cur, W_cur,eta,
#                                                                           index_obs,Y_here,need_sample,need_mean,need_var)
#     Y_sample_return[i_row,]= KF_smoother_sample_list[[1]]
#     
#     #plot(Y[i_row,])
#     #lines(Y_sample_return[i_row,],type='l')
#   }
#   return(Y_sample_return)
#   
# }
# 
# Intial_Y_sample_LM<-function(Y,grid_with_full_obs_Y,input2){
# 
# 
#   
#   Y_sample_return=Y
#   
#   for(i_row in 1:n1){
#     index_missing_row=which(is.na(Y[i_row,])==T)
#     
#     m=lm(Y[i_row,-index_missing_row]~input2[-index_missing_row])
#     Y_sample_return[i_row,index_missing_row]=m$coefficients[1]+input2[index_missing_row]*m$coefficients[2]
#     
#   }
#   return(Y_sample_return)
#   
# }

##this is to update the posterior distribution of B given GP prior
##this is to implement at the beginning of the cor
##let me make this one dim

Construct_Cov_KF<-function(log_beta_eta,delta_x,index_obs,kernel_type='exp'){
  #num_kernel=length(log_beta_eta)/2
  return_list=list()
  #for(i_ker in 1:num_kernel){

  if(kernel_type=='exp'){
    lambda=1.0*exp(log_beta_eta[1])
    eta=exp(log_beta_eta[2])
    
    #return_list=list()
    return_list$GG=Construct_G_exp( delta_x,  lambda)
    return_list$W0=Construct_W0_exp(1, lambda)
    return_list$W=Construct_W_exp(1, delta_x, lambda, return_list$W0) 
  
    
    return_list$C_R_K_Q=Get_C_R_K_Q_pred(index_obs, return_list$GG,return_list$W,
                                             return_list$W0, eta )
  }else if(kernel_type=='matern_5_2'){
    lambda=sqrt(5)*exp(log_beta_eta[1])
    eta=exp(log_beta_eta[2])
    
    return_list$GG=Construct_G_matern_5_2( delta_x,  lambda)
    return_list$W0=Construct_W0_matern_5_2(1, lambda)
    return_list$W=Construct_W_matern_5_2(1, delta_x, lambda, return_list$W0) 
    
    
    return_list$C_R_K_Q=Get_C_R_K_Q_pred(index_obs, return_list$GG,return_list$W,
                                         return_list$W0, eta )
    
  }
  
  #}
  
  return(return_list)
  
}

Get_eigen_1<-function(R0_1,log_beta1,kernel_type='exp',d){
  if(kernel_type=='exp'){
    R1=pow_exp_funct(R0_1,exp(log_beta1),1.0)
    eigen_R1=eigs(R1,k=d,symmetric=T) ## only compute the first d
    
    
  }else if(kernel_type=='matern_5_2'){
    R1=matern_5_2_funct(R0_1,exp(log_beta1))
    eigen_R1=eigs(R1,k=d,symmetric=T) ## only compute the first d
  }
  return(eigen_R1)
  
}

##this is to update at the end
# Update_beta_post_S_2_all<-function(Y_sample,Y_sample_tilde,Cov_2,tilde_eta_vec,log_beta_eta,sigma_2_0_cur){
#   return_list=list()
#   #return_list$log_det=rep(0,d)
#   return_list$S_2=rep(0,d)
#   return_list$log_post_cur=rep(0,d)
#   
#   Trace_Y_sample_tilde_2=0
#   for(i_row in 1:d){
#     Trace_Y_sample_tilde_2=Trace_Y_sample_tilde_2+sum(Y_sample_tilde[i_row,]*Y_sample_tilde[i_row,])
#     
#     L_inv_y_tilde=Get_L_inv_y(Cov_2[[i_row]]$GG,tilde_eta_vec[i_row], as.matrix(Cov_2[[i_row]]$C_R_K_Q[[4]],n2,1), Cov_2[[i_row]]$C_R_K_Q[[3]], t(Y_sample_tilde[i_row,] ))
#     
#     return_list$S_2[i_row]=sum(L_inv_y_tilde*(L_inv_y_tilde) )*tilde_eta_vec[i_row] ###eigen_R1_values[i_row]*eta[i_row]
#     #return_list$log_det[i_row]=sum(log(Cov_2[[i_row]]$C_R_K_Q[[4]]))
#     log_det=sum(log(Cov_2[[i_row]]$C_R_K_Q[[4]]))
#     
#     
#     log_JR_prior_cur=log_approx_ref_prior(log_beta_eta[c(i_row,i_row+d)],0,T,CL[2],a,b)
#     ##terms related to 
#     return_list$log_post_cur[i_row]=n2/2*(log(tilde_eta_vec_cur[i_row]))-1/2*log_det-return_list$S_2[i_row]/(2*sigma_2_0_cur)+log_JR_prior_cur+sum(log_beta_eta[c(i_row,i_row+d)])
#     
#   }
#   
#   return(return_list)
# }

Get_beta_post_S_2_log_det<-function(Y_sample_tilde_cur_i_row_t,log_beta_eta,tilde_eta,Cov_here,add_log_prior=T){
  
  
  L_inv_y_tilde=Get_L_inv_y(Cov_here$GG,tilde_eta, as.matrix(Cov_here$C_R_K_Q[[4]],n2,1), Cov_here$C_R_K_Q[[3]], Y_sample_tilde_cur_i_row_t)
  
  S_2=sum(L_inv_y_tilde*(L_inv_y_tilde) )*tilde_eta ###eigen_R1_values[i_row]*eta[i_row]
  
  log_det=sum(log(Cov_here$C_R_K_Q[[4]]))
  #log_JR_prior=0
  if(add_log_prior){
     log_JR_prior=log_approx_ref_prior(log_beta_eta,0,T,CL[2],a,b)
     log_post=n2/2*(log(tilde_eta))-1/2*log_det-S_2/(2*sigma_2_0_cur)+log_JR_prior+sum(log_beta_eta)
  }else{
    log_post=n2/2*(log(tilde_eta))-1/2*log_det-S_2/(2*sigma_2_0_cur)
    
  }
  ##terms related to 

  return(c(log_post,S_2,log_det))
}


Get_beta_post_S_2_log_det_ori<-function(Y_sample_tilde_cur_i_row_t,log_beta_eta,Cov_here){
  
  beta=exp(log_beta_eta)[1]
  eta=exp(log_beta_eta)[2]
  
  L_inv_y_tilde=Get_L_inv_y(Cov_here$GG,eta, as.matrix(Cov_here$C_R_K_Q[[4]],n2,1), Cov_here$C_R_K_Q[[3]], Y_sample_tilde_cur_i_row_t)
  
  S_2=sum(L_inv_y_tilde*(L_inv_y_tilde) )*eta ###eigen_R1_values[i_row]*eta[i_row]
  
  log_det=sum(log(Cov_here$C_R_K_Q[[4]]))

  log_JR_prior=log_approx_ref_prior(log_beta_eta,0,T,CL[2],a,b)
  
  ##terms related to 
  log_post=n2/2*(log(eta))-1/2*log_det-S_2/(2*sigma_2_0_cur)+log_JR_prior+sum(log_beta_eta)
  return(c(log_post,S_2,log_det))
  
}
  


Update_B_post_GP<-function(Cov_B_here,S_r_here, X,log_beta_eta,Y_sample_tilde_B,index_obs){
  #q_B=length(log_beta_eta_B_cur)/2
  

  #for(i_B in 1:q_B){
    
    #S_r=as.numeric(t(X)%*%X)
    beta=exp(log_beta_eta[1])
    eta= exp(log_beta_eta[2])
    
    ###only this part needs to be reconstructed if the return list is known
    

    S_xy=t(X)%*%Y_sample_tilde_B
    S2_xy=(S_xy)%*%t(S_xy)
    L_inv_y_tilde_B=Get_L_inv_y(Cov_B_here$GG, eta , as.matrix(Cov_B_here$C_R_K_Q[[4]],n2,1), Cov_B_here$C_R_K_Q[[3]], t(S_xy))
    
    S2_B_cur=(S2_xy-t(L_inv_y_tilde_B)%*%L_inv_y_tilde_B*eta)/S_r_here
    log_det_B=sum(log(Cov_B_here$C_R_K_Q[[4]]))
    
    log_JR_prior_B=log_approx_ref_prior(log_beta_eta,0,T,CL[2],a,b)
    
    log_post_B=-1/2*log_det_B+n2/2*(eta)-S2_B_cur+log_JR_prior_B+sum(log_beta_eta)
    
  #}
  
  return(log_post_B)
  
  
}

# Cov_B_cur=Construct_Cov_KF(log_beta_eta_B_cur,delta_x_2,index_obs)
# 
# log_post_B_cur=rep(0,2)  
# for(i_B in 1:2){
# 
#   S_r=as.numeric(t(H1[,i_B])%*%H1[,i_B])
#   
#   
#   lambda_B=1.0*exp(log_beta_eta_B_cur[i_B])
#   
#   GG_2_B_cur=Construct_G_exp( delta_x_2,  lambda_B)
#   W0_2_B_cur=Construct_W0_exp(1, lambda_B)
#   W_2_B_cur=Construct_W_exp(1, delta_x_2, lambda_B, W0_2_B_cur) 
#   
#   
#   C_R_K_Q_B_cur=Get_C_R_K_Q_pred(index_obs, GG_2_B_cur,W_2_B_cur,
#                                  W0_2_B_cur, 
#                                  exp(log_beta_eta_B_cur[i_B+2]) );
#   
#   
#   ###only this part needs to be reconstruct
#   
#   Y_sample_tilde_B=Y_sample-(U1_cur)%*%(z_sample_matrix)
#   S_xy=t(H1[,i_B])%*%Y_sample_tilde_B
#   S2_xy=(S_xy)%*%t(S_xy)
#   L_inv_y_tilde_B=Get_L_inv_y(GG_2_B_cur, exp(log_beta_eta_B_cur[i_B+2]) , as.matrix(C_R_K_Q_B_cur[[4]],n2,1), C_R_K_Q_B_cur[[3]], t(S_xy))
#   
#   S2_B_cur=(S2_xy-t(L_inv_y_tilde_B)%*%L_inv_y_tilde_B*exp(log_beta_eta_B_cur[i_B+2]))/S_r
#   log_det_B_cur=sum(log(C_R_K_Q_B_cur[[4]]))
#   
#   log_JR_prior_B_cur=log_approx_ref_prior(log_beta_eta_B_cur[c(i_B,i_B+2)],0,T,CL[2],a,b)
#   
#   log_post_B_cur[i_B]=-1/2*log_det_B_cur+n2/2*(exp(log_beta_eta_B_cur[i_B+2]))-S2_B_cur+log_JR_prior_B_cur+sum(log_beta_eta_B_cur[c(i_B,i_B+2)])
# 
# }
# 

branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  ##########################################################################
  #
  # BRANIN FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2)
  # a = constant (optional), with default value 1
  # b = constant (optional), with default value 5.1/(4*pi^2)
  # c = constant (optional), with default value 5/pi
  # r = constant (optional), with default value 6
  # s = constant (optional), with default value 10
  # t = constant (optional), with default value 1/(8*pi)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}

braninsc <- function(xx)
{
  ##########################################################################
  #
  # BRANIN FUNCTION, RESCALED
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  x1bar <- 15*x1 - 5
  x2bar <- 15 * x2
  
  term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
  term2 <- (10 - 10/(8*pi)) * cos(x1bar)
  
  y <- (term1^2 + term2 - 44.81) / 51.95
  return(y)
}
