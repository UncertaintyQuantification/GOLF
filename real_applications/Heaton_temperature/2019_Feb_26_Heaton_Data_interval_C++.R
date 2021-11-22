library(ggplot2)
library(plot3D)
load("Data_Heaton/AllSatelliteTemps.RData")
library(RobustGaSP)
library(FastGaSP)
library(Rcpp)
library(RcppEigen)
sourceCpp(file='../../src/functions_Feb_26_2020.cpp') 
library(maps)#Install maps package if not done before
library(rARPACK)
source(file='../../functions/functions.R')
library(Rfast)


seed_here=0
set.seed(seed_here)

dat = all.sat.temps

Masktemp_mat=matrix(dat[,3],500,300)
Truetemp_mat=matrix(dat[,4],500,300)

#image2D(Masktemp_mat)
#image2D(Truetemp_mat)


sd(dat$TrueTemp, na.rm=T) # standard deviation of True temperatures

## get the lattitude and longitude
lat_unique = unique(dat$Lat) # from 34.29519 to 37.06811
#lat_unique
lon_unique = unique(dat$Lon ) # from -95.91153 to -91.28381
#lon_unique
xy.grid = expand.grid(lat_unique,lon_unique)

temp_lat=colMeans(Truetemp_mat,na.rm = T)
#plot(lat_unique,temp_lat)

temp_lon=rowMeans(Truetemp_mat,na.rm = T)
#plot(lon_unique,temp_lon)



index_lon=1:500
index_lat=300:1


Masktemp_mat_sub=Masktemp_mat[index_lon,index_lat]
Truetemp_mat_sub=Truetemp_mat[index_lon,index_lat]

lon_unique_sub=lon_unique[index_lon]
lat_unique_sub=lat_unique[index_lat]

#image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='',xlab='',main='Observed temperature',
#        mgp=c(3,0.5,0))
#image2D(Truetemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='',xlab='',mgp=c(3,0.5,0),
#        main='Full temperature')


n1=length(lat_unique_sub)
n2=length(lon_unique_sub)

input1=lat_unique_sub ##row
input2=lon_unique_sub##column

input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)

input=cbind(rep(input1,n2),input_2_vec)


Y_full=t(Truetemp_mat_sub)
grid_with_full_obs_Y_full=which(is.na(Y_full)==F)

Y=t(Masktemp_mat_sub)
grid_with_full_obs_Y=which(is.na(Y)==F)

missing_index_all=(1:(n1*n2))[-grid_with_full_obs_Y]  ##this needs to be updated

missing_index=setdiff(grid_with_full_obs_Y_full, grid_with_full_obs_Y)


n_missing=length(missing_index)
n_missing_all=length(missing_index_all)

mean_not_missing=mean(Y[grid_with_full_obs_Y])



Y_sample=Y

##use column mean to initiate the initial value
index_missing_col=list()
#
for(i in 1:n1){
  index_missing_col[[i]]=which(is.na(Y[i,])==T)
  Y_sample[i, index_missing_col[[i]]]=mean(Y[i,],na.rm=T)
  #Y_sample[i, index_missing_col[[i]]]=mean(Y[i,],na.rm=T)+  0.1*sd(Y_sample[grid_with_full_obs_Y])*rnorm(length(index_missing_col[[i]]))
}

#Y_sample[-grid_with_full_obs_Y]= mean(Y_sample[grid_with_full_obs_Y])+0.1*sd(Y_sample[grid_with_full_obs_Y])*rnorm(n1*n2-length(grid_with_full_obs_Y))
#Y_sample[-grid_with_full_obs_Y]= sample(Y_sample[grid_with_full_obs_Y],n_missing_all)

Y_sample_1=Y_sample

sum(is.na(Y_sample)) ##should see zero
#image2D(Y_sample)


##number of sample
#M=40000
#M_0=8000
M=6000
M_0=1200

interval_prop=0.025
need_interval=T

##number of d
d=150


##default jointly robust prior
a=1/2-2
b=1

CL = rep(0,2)    ###CL is also used in the prior so I make it a model parameter

CL[1] = (max(input[,1])-min(input[,1]))/n1
CL[2] = (max(input[,2])-min(input[,2]))/n2

##prior parameter CL[1:2], a, b
prior_par=c(CL,a,b)

##the more sample you expect a short step size? (because the likelihood is more concentrate)
#sd_log_beta=sd_log_eta=20/n2
#sd_log_beta1=10/n1

sd_log_beta=sd_log_eta=40/n2
sd_log_beta1=40/n1

step_size_kernel_par=c(sd_log_beta,sd_log_eta,sd_log_beta1)






##initial value
 log_beta1_cur=3
 log_beta_eta_cur=c(rep(0,d),rep(-3,d))  ##this is for the second one
#log_beta1_cur=-1+2*runif(1)
#log_beta_eta_cur=c(-1+2*runif(d),rep(-3,d))  ##this is for the second one



initial_kernel_par=c(log_beta_eta_cur,log_beta1_cur)

delta_x_2=input2[2:n2]-input2[1:(n2-1)]

R0_1=abs(outer(input1,input1,'-'))

kernel_type='exp'
#kernel_type='matern_5_2'
##okay reimplement this
missing_index_all_mat=matrix(0,length(missing_index_all),2)
missing_index_all_mat[,1]=as.integer(missing_index_all%%n1-1)

missing_index_all_mat[,2]=as.integer(floor(missing_index_all/n1))

# return_list=GOLF_2D_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
#            M,M_0, d, missing_index, missing_index_all_mat,
#           prior_par,initial_kernel_par,step_size_kernel_par,seed_here)

H2=matrix(1,n2,1)
have_mean=T
record_time=system.time(
  for(ii in 1:1){
    return_list=GOLF_2D_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                             M, M_0, d, missing_index, missing_index_all_mat,
                             prior_par,initial_kernel_par,step_size_kernel_par,seed_here,
                             need_interval,interval_prop,have_mean,H2 )
  }
)

#save(return_list, record_time,file="March_17_mean_d_150_M_40000_M0_8000_M_0_revised_C++_seed_0_40_40_uniformly_sampled_initial_pars.RData")


Y_avg_mean=return_list[[1]]
#image2D(Y_avg_mean)
#image2D(Y_sample)
rmse=sqrt(mean( (Y_avg_mean[missing_index]-Y_full[missing_index])^2))

diff_set=setdiff( missing_index_all,missing_index)

index_pred=rep(0,length(missing_index))

i_missing_all=1
for(i_missing in 1:length(missing_index)){
  while(missing_index_all[i_missing_all]!=missing_index[i_missing]){
    i_missing_all=i_missing_all+1
  }
  index_pred[i_missing]=i_missing_all
  
}
max(abs(missing_index_all[index_pred]-missing_index))



lower_95=return_list[[4]][index_pred]
upper_95=return_list[[5]][index_pred]
#upper_95-rowMins(return_list[[4]][index_pred,],value=T)
prop_95=sum(lower_95<=Y_full[missing_index] &upper_95>=Y_full[missing_index])/length(missing_index)
#length(which(lower_95<=Y_full[missing_index] &upper_95>=Y_full[missing_index]))/length(missing_index)
length_95=mean(upper_95-lower_95)

rmse
prop_95
length_95
record_time[3]/60

plot(return_list[[3]]/n1,type='l')

position=2
return_list[[3]][position,]
return_list[[4]][position,]
#sort(return_list[[5]][position,])

plot(log(return_list[[2]][M_0:M,1]))
plot(log(return_list[[2]][M_0:M,2]))
plot(log(return_list[[2]][M_0:M,3]))
plot(log(return_list[[2]][M_0:M,4]))

plot(log(return_list[[2]][M_0:M,151]))
plot(log(return_list[[2]][M_0:M,152]))
plot(log(return_list[[2]][M_0:M,300]))

plot(log(return_list[[2]][M_0:M,301]))

plot((return_list[[2]][M_0:M,302]))





########
###change some of the initial parameters, test for sensitivity
M=15000
M_0=3000
###need to change the step size when correlation are shared. We need a smaller step size
sd_log_beta=.02  ##this may depend on the number of n
sd_log_beta1=.02
sd_log_eta=0.1
step_size_kernel_par=(c(sd_log_beta,sd_log_eta,sd_log_beta1))
###the number of initial parameters should also be changed
log_beta1_cur=3
log_beta_eta_cur=c(0,-3)  ##this is for the second one
initial_kernel_par=(c(log_beta_eta_cur,log_beta1_cur))



system.time(
  for(ii in 1:1){
    return_list_separable=GOLF_2D_separable_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type,
                                                 M, M_0, d, missing_index, missing_index_all_mat,
                                                 prior_par,initial_kernel_par,step_size_kernel_par,seed_here,
                                                 need_interval,interval_prop,have_mean, H2)
  }
)

##
Y_avg_mean_separable=return_list_separable[[1]]
image2D(Y_avg_mean_separable)
image2D(Y)

sqrt(mean(((Y_avg_mean_separable[missing_index])-Y_full[missing_index])^2))

lower_95_separable=return_list_separable[[4]][index_pred]
upper_95_separable=return_list_separable[[5]][index_pred]

sum(lower_95_separable<=Y_full[missing_index] &upper_95_separable>=Y_full[missing_index])/length(missing_index)
#length(which(lower_95<=Y_full[missing_index] &upper_95>=Y_full[missing_index]))/length(missing_index)
mean(upper_95_separable-lower_95_separable)



##plot it to see whether the step size is okay
 par(mfrow=c(1,3))
 plot(log(return_list_separable[[2]][M_0:M,1]),type='l')
 plot(log(return_list_separable[[2]][M_0:M,2]),type='l')
 plot(log(return_list_separable[[2]][M_0:M,3]),type='l')
# dev.off()

 
 
 ##56 24
 max(Y_avg_mean)
 min(Y_avg_mean)
 max(Truetemp_mat,na.rm=T)
 min(Truetemp_mat,na.rm=T)
 
 png("pred_temp_golf.png",height=4,width=5.5,units="in",res=200)
 image2D(t(Y_avg_mean),lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
         main='GOLF',zlim=c(24,56))
 dev.off()
 
 
 png("obs_temp.png",height=4,width=5.5,units="in",res=200)
 image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
         main='Observed temperature',zlim=c(24,56))
 dev.off()
 
 png("full_temp.png",height=4,width=5.5,units="in",res=200)
 image2D(Truetemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
         main='Full temperature',zlim=c(24,56))
 dev.off()
 
 #par(mfrow=c(1,3))
 pdf('full_Data.pdf',height=4.5,width=4.5)
 image2D(Truetemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
         main='Full data',zlim=c(24,56),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,cex.sub=1.5,colkey = list(plot = FALSE, side = 3))
 dev.off()
 pdf('Observed_data.pdf',height=4.5,width=4.5)
 
 image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
         main='Observed data',zlim=c(24,56),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,cex.sub=1.5,colkey = list(plot = FALSE, side = 3))
 dev.off()
 pdf('Predictive_mean.pdf',height=4.5,width=4.7)
 image2D(t(Y_avg_mean),lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
         main='Predictive mean',zlim=c(24,56),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,cex.sub=1.5,colkey = list(plot = TRUE, side = 4, cex.axis=1.5))
 dev.off()
 
 #colkey(side = 4, add = TRUE, clim = c(24,56))
 
