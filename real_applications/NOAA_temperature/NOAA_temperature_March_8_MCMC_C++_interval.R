library(Rcpp)
library(RcppEigen)
library(plot3D)
library(maps)#Install maps package if not done before
sourceCpp(file='../../src/functions_March_8_2020.cpp') 

gridded_data=readRDS("NOAA_gridded_temperature/NOAA_gridded_data_1999_2018_GOLF.rds")

gridded_data_temp_array=array(gridded_data$temperature,c(72,36,240 ))


##na count 
record_na=rep(0,36)
for(i in 1:36){  
  record_na[i]=sum(is.na(gridded_data_temp_array[,i,]))
}

gridded_data_temp_array_sub=gridded_data_temp_array[37:72,6:33,]
Y_full=matrix(gridded_data_temp_array_sub,36*28,240)

# #as.vector(gridded_data_temp_array[,,1])-gridded_data$temperature[,1]
#   
# Y_full=gridded_data$temperature[(3*72+1):(72*36-72),]
# 
# index_here=NULL
# for(i_index in 1:32){
#   index_here=c(index_here,(1:36)+72*(i_index-1))
# }
# Y_full=Y_full[index_here, ]
time=gridded_data$months
time[150]
#Lat= seq(-87.5, 87.5, length=36)
#Lat= seq(-72.5, 82.5, length=32)
Lat= seq(-62.5, 72.5, length=28)

#Lon=seq(2.5, 357.5, length=72)
#Lon=seq(2.5, 177.5, length=36)
Lon=seq(182.5, 357.5, length=36)

i_plot=180
mapmat=matrix(Y_full[,i_plot],36,28)
#mapmat=matrix(gridded_data$temperature[,i_plot],72,36)

max(mapmat,na.rm=T)
min(mapmat,na.rm=T)
int=seq(-7,5,length.out=81)

rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="full temperature anomalies",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})

##sample 20 % as missing value to test
seed_here=0
set.seed(seed_here)
additional_proportion_missing=0.5

n11=length(Lon)
n12=length(Lat)
n1=n11*n12
n2=length(gridded_data$months)

index_with_obs_full= which(!is.na(Y_full))
n_obs_full=length(index_with_obs_full)
missing_index=sample(index_with_obs_full,n_obs_full*additional_proportion_missing)
Y=Y_full
Y[missing_index]=NA
missing_index_all=which(is.na(Y))
n_missing=length(missing_index)
n_missing_all=length(missing_index_all)

mapmat=matrix(Y[,i_plot],36,28)
max(mapmat,na.rm=T)
min(mapmat,na.rm=T)
int=seq(-8,10,length.out=81)

rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="obs temperature anomalies",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})


#input11=Lat 
#input12=Lon 

input11=Lon 
input12=Lat


input2=c(1:n2)

delta_x_2=input2[2:n2]-input2[1:(n2-1)]


Y_sample=Y

##use row mean to initiate the initial value
index_missing_col=list()
record_all_NA_row=NULL
#
for(i in 1:n1){
  index_missing_col[[i]]=which(is.na(Y[i,])==T)
  if(is.na(mean(Y[i,],na.rm=T)) ){
    record_all_NA_row=c(record_all_NA_row,i)
  }else{
    Y_sample[ i, index_missing_col[[i]]]=mean(Y[i,],na.rm=T)
  }
}

index_missing_row=list()

for(j in 1:n2){
  index_missing_row[[j]]=which(is.na(Y_sample[,j])==T)
  Y_sample[index_missing_row[[j]],j]=mean(Y[,j],na.rm=T)
  
}


Y_sample_1=Y_sample

sum(is.na(Y_sample)) ##should see zero


##number of sample
M=3000 ###converge very fast for this example
M_0=600
interval_prop=0.025
need_interval=T


##number of d
n11
n12
d1=n11*3/4
d2=n12*3/4


#d1_plus_d2=d1+d2
d=d1*d2
d_vec=as.vector(c(d1,d2))
d_vec=as.integer(d_vec)

##default jointly robust prior
a=1/2-2
b=1


CL1 = rep(0,2)    ###CL is also used in the prior so I make it a model parameter


CL1[1] = (max(input11)-min(input11))/n11
CL1[2] = (max(input12)-min(input12))/n12

CL2=(max(input2)-min(input2))/n2
CL=c(CL1,CL2)

sd_log_beta=sd_log_eta=40/n2
sd_log_beta11=40/n11^2
sd_log_beta12=40/n12^2 ##we are using matern here so 
sd_log_beta1=c(sd_log_beta11,sd_log_beta12)


prior_par=c(CL,a,b)

step_size_kernel_par=c(sd_log_beta,sd_log_eta,sd_log_beta1)

log_beta11_cur=log(1/CL1[1])
log_beta12_cur=log(1/CL1[2])
log_beta1_cur=c(log_beta11_cur,log_beta12_cur)

log_beta_eta_cur=c(rep(log(1/mean(delta_x_2)),d),rep(-3,d))  ##this is for the second one

 # log_beta11_cur=log_beta12_cur=-3
 # log_beta1_cur=c(log_beta11_cur,log_beta12_cur)
 # 
 # log_beta_eta_cur=c(rep(0,d),rep(-3,d))  ##this is for the second one
# 
initial_kernel_par=c(log_beta_eta_cur,log_beta1_cur)

delta_x_2=input2[2:n2]-input2[1:(n2-1)]

###let me use exponential here
R0_11=abs(outer(input11,input11,'-'))
R0_12=abs(outer(input12,input12,'-'))
R0_1=as.list(1:2)
R0_1[[1]]=R0_11
R0_1[[2]]=R0_12

kernel_type_spatial='matern_5_2'
kernel_type_temporal='exp'
##okay reimplement this
missing_index_all_mat=matrix(0,length(missing_index_all),2)
missing_index_all_mat[,1]=as.integer(missing_index_all%%n1-1)

missing_index_all_mat[,2]=as.integer(floor(missing_index_all/n1))


H2=cbind(rep(1,n2),input2)

#H2=as.matrix(input2)
q2=dim(H2)[2]
have_mean=T

record_time_GOLF=system.time(
  for(ii in 1:1){
    return_list=GOLF_3D_MCMC(Y_sample, R0_1,delta_x_2,  kernel_type_spatial,kernel_type_temporal,
                             M, M_0, d_vec, missing_index, missing_index_all_mat,
                             prior_par,initial_kernel_par,step_size_kernel_par,seed_here,
                             need_interval,interval_prop,have_mean,H2)
    

  }
)

#sqrt(mean( (return_list[[1]][missing_index]-Y_full[missing_index])^2))
#sqrt(mean( (return_list[[1]][missing_index]-Y_full[missing_index])^2))


##using built in R to sample
#user   system  elapsed 
#2578.975    8.578 2587.252 

##using not built in R
#user   system  elapsed 
#2627.155    8.719 2635.596 

#> sqrt(mean( (Y_avg_mean[missing_index]-Y_full[missing_index])^2))
#[1] 0.3256277
Y_avg_mean=return_list[[1]]
#image2D(Y_avg_mean)
#image2D(Y_sample)
sqrt(mean( (Y_avg_mean[missing_index]-Y_full[missing_index])^2))

diff_set=setdiff( missing_index_all,missing_index)

missing_index_sort=sort(missing_index)
index_pred=rep(0,length(missing_index_sort))

i_missing_all=1
for(i_missing in 1:length(missing_index_sort)){
  while(missing_index_all[i_missing_all]!=missing_index_sort[i_missing]){
    i_missing_all=i_missing_all+1
  }
  index_pred[i_missing]=i_missing_all
  
}
max(abs(missing_index_all[index_pred]-missing_index))

lower_95=return_list[[4]][index_pred]
upper_95=return_list[[5]][index_pred]
#upper_95-rowMins(return_list[[4]][index_pred,],value=T)
sum(lower_95<=Y_full[missing_index_sort] &upper_95>=Y_full[missing_index_sort])/length(missing_index_sort)
#length(which(lower_95<=Y_full[missing_index] &upper_95>=Y_full[missing_index]))/length(missing_index)
mean(upper_95-lower_95)
#dim(return_list[[3]])
#dim(return_list[[4]])

dim(return_list[[2]])
#

plot(log(return_list[[2]][M_0:M,1]))

plot(log(return_list[[2]][M_0:M,2*d+1]),type='l')
plot(log(return_list[[2]][M_0:M,2*d+2]),type='l')




##1. GP to interpolate the grid for each missing year
library(RobustGaSP)
##here the grid is 72*36
input_12_vec=matrix(t(matrix(rep(Lat,n11),n12,n11)),n12*n11)

input=cbind(rep(Lon,n12),input_12_vec)


predict_rgasp_list=as.list(1:n2)
record_time_spatial_rgasp=system.time(
  for(i_month in 1:n2){
    print(i_month)
    
    index_obs=which(is.na(Y[,i_month])==F)
    Y_obs=Y[index_obs,i_month]
    input_obs=input[index_obs,]
    
    #model.rgasp=rgasp(design=input_obs,response=as.vector(Y_obs),nugget.est=T)
    model.rgasp=rgasp(design=input_obs,response=as.vector(Y_obs),nugget.est=T,num_initial_values=1)
    
    pred.rgasp.all=predict(model.rgasp,input)
    predict_rgasp_list[[i_month]]=pred.rgasp.all
    
    #image2D(matrix(pred.rgasp.all$mean,36,28))
    #image2D(matrix(Y_full[,i_month],36,28))
    
  }
)

n_missing_record=rep(0,n2)

sum_of_squares_record_rgasp=rep(0,n2)
coverage_num_record_rgasp=rep(0,n2)
total_length_interval_record_rgasp=rep(0,n2)

for(i_month in 1:n2){
  index_missing_ori=which( (missing_index>= 1+(i_month-1)*n1)&(missing_index<= (i_month)*n1))
  missing_index_here=missing_index[index_missing_ori]-(i_month-1)*n1
  n_missing_record[i_month]=length(missing_index_here)
  sum_of_squares_record_rgasp[i_month]=sum( (  predict_rgasp_list[[i_month]]$mean[missing_index_here]-Y_full[missing_index[index_missing_ori]])^2)
  coverage_num_record_rgasp[i_month]=length(which(Y_full[missing_index[index_missing_ori]]<=   predict_rgasp_list[[i_month]]$upper95[missing_index_here] &Y_full[missing_index[index_missing_ori]]>=   predict_rgasp_list[[i_month]]$lower95[missing_index_here]))
  total_length_interval_record_rgasp[i_month]=sum(abs(  predict_rgasp_list[[i_month]]$upper95-  predict_rgasp_list[[i_month]]$lower95))
  
}



#sqrt(sum_of_squares_record[1:n2]/n_missing_record[1:n2])

sqrt( sum(sum_of_squares_record_rgasp[1:n2])/sum(n_missing_record[1:n2]) )
sum(coverage_num_record_rgasp)/sum(n_missing_record)
sum(total_length_interval_record_rgasp)/sum(n_missing_record)
record_time_spatial_rgasp
#total_length_interval_record_rgasp/n_missing_record


##1. GP to interpolate the grid for each missing year
library(RobustGaSP)
##here the grid is 72*36
input_12_vec=matrix(t(matrix(rep(Lat,n11),n12,n11)),n12*n11)

input=cbind(rep(Lon,n12),input_12_vec)


predict_rgasp_list=as.list(1:n2)
record_time_spatial_rgasp2=system.time(
  for(i_month in 1:n2){
    print(i_month)
    
    index_obs=which(is.na(Y[,i_month])==F)
    Y_obs=Y[index_obs,i_month]
    input_obs=input[index_obs,]
    
    #model.rgasp=rgasp(design=input_obs,response=as.vector(Y_obs),nugget.est=T)
    model.rgasp=rgasp(design=input_obs,response=as.vector(Y_obs),nugget.est=T,num_initial_values=2)
    
    pred.rgasp.all=predict(model.rgasp,input)
    predict_rgasp_list[[i_month]]=pred.rgasp.all
    
    #image2D(matrix(pred.rgasp.all$mean,36,28))
    #image2D(matrix(Y_full[,i_month],36,28))
    
  }
)

n_missing_record=rep(0,n2)

sum_of_squares_record_rgasp2=rep(0,n2)
coverage_num_record_rgasp2=rep(0,n2)
total_length_interval_record_rgasp2=rep(0,n2)

for(i_month in 1:n2){
  index_missing_ori=which( (missing_index>= 1+(i_month-1)*n1)&(missing_index<= (i_month)*n1))
  missing_index_here=missing_index[index_missing_ori]-(i_month-1)*n1
  n_missing_record[i_month]=length(missing_index_here)
  sum_of_squares_record_rgasp2[i_month]=sum( (  predict_rgasp_list[[i_month]]$mean[missing_index_here]-Y_full[missing_index[index_missing_ori]])^2)
  coverage_num_record_rgasp2[i_month]=length(which(Y_full[missing_index[index_missing_ori]]<=   predict_rgasp_list[[i_month]]$upper95[missing_index_here] &Y_full[missing_index[index_missing_ori]]>=   predict_rgasp_list[[i_month]]$lower95[missing_index_here]))
  total_length_interval_record_rgasp2[i_month]=sum(abs(  predict_rgasp_list[[i_month]]$upper95-  predict_rgasp_list[[i_month]]$lower95))
  
}



#sqrt(sum_of_squares_record[1:n2]/n_missing_record[1:n2])

sqrt( sum(sum_of_squares_record_rgasp2[1:n2])/sum(n_missing_record[1:n2]) )
sum(coverage_num_record_rgasp2)/sum(n_missing_record)
sum(total_length_interval_record_rgasp2)/sum(n_missing_record)
record_time_spatial_rgasp2



par(mfrow=c(1,2))    
i_plot=229
#mapmat=matrix(Y_sample[,i_plot],36,32)
mapmat=matrix(Y_full[,i_plot],36,28)

max(mapmat,na.rm=T)
min(mapmat,na.rm=T)
int=seq(-6,6,length.out=81)
pdf("full_temp.pdf",height=6,width=6)

rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Temperature anomalies in  Jan 2018",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})
dev.off()
mapmat_golf=matrix(Y_avg_mean[,i_plot],36,28)
max(mapmat_golf,na.rm=T)
min(mapmat_golf,na.rm=T)
int=seq(-6,6,length.out=81)

pdf("golf_temp.pdf",height=6,width=6)

rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat_golf, color.palette=rgb.palette, levels=int,
               plot.title=title(main="GOLF",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})
dev.off()



#png("obs_temp.png",height=4,width=5.5,units="in",res=200)
pdf("spatial_temp.pdf",height=6,width=6)

mapmat_spatial2=matrix(predict_rgasp_list[[i_plot]]$mean,36,28)
max(mapmat_spatial2,na.rm=T)
min(mapmat_spatial2,na.rm=T)
int=seq(-6,6,length.out=81)

rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat_spatial2, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Spatial model",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})
dev.off()

# ###spatial by dicekriging
# library(DiceKriging)
# 
# predict_dk_list=as.list(1:n2)
# record_time_spatial_dk=system.time(
#   for(i_month in 1:n2){
#     print(i_month)
#     
#     index_obs=which(is.na(Y[,i_month])==F)
#     Y_obs=Y[index_obs,i_month]
#     input_obs=input[index_obs,]
#     
#     model.dk=km(design=input_obs,response=as.vector(Y_obs),nugget.estim=T)
#     pred.dk.all=predict(model.dk,input,type='UK')
#     predict_dk_list[[i_month]]=pred.dk.all
#     
#     #image2D(matrix(pred.rgasp.all$mean,36,28))
#     #image2D(matrix(Y_full[,i_month],36,28))
#     
#   }
# )
# 
# 
# 
# #n_missing_record=rep(0,n2)
# 
# sum_of_squares_record_dk=rep(0,n2)
# coverage_num_record_dk=rep(0,n2)
# total_length_interval_record_dk=rep(0,n2)
# 
# for(i_month in 1:n2){
#   index_missing_ori=which( (missing_index>= 1+(i_month-1)*n1)&(missing_index<= (i_month)*n1))
#   missing_index_here=missing_index[index_missing_ori]-(i_month-1)*n1
#   #n_missing_record[i_month]=length(missing_index_here)
#   sum_of_squares_record_dk[i_month]=sum( (  predict_dk_list[[i_month]]$mean[missing_index_here]-Y_full[missing_index[index_missing_ori]])^2)
#   coverage_num_record_dk[i_month]=length(which(Y_full[missing_index[index_missing_ori]]<=   predict_dk_list[[i_month]]$upper95[missing_index_here] &Y_full[missing_index[index_missing_ori]]>=   predict_dk_list[[i_month]]$lower95[missing_index_here]))
#   total_length_interval_record_dk[i_month]=sum(abs(  predict_dk_list[[i_month]]$upper95-  predict_dk_list[[i_month]]$lower95))
#   
# }
# 
# sqrt( sum(sum_of_squares_record_dk[1:n2])/sum(n_missing_record[1:n2]) )
# sum(coverage_num_record_dk)/sum(n_missing_record)
# sum(total_length_interval_record_dk)/sum(n_missing_record)
# record_time_spatial_dk

###lagp
library(laGP)
input_12_vec=matrix(t(matrix(rep(Lat,n11),n12,n11)),n12*n11)

input1=cbind(rep(Lon,n12),input_12_vec)
input_all=matrix(0,n1*n2,3)

for(i in 1:n2){
  input_all[(1+(i-1)*n1):(i*n1), ]=cbind(input1,rep(i,n1))
}

index_obs=setdiff(index_with_obs_full,missing_index)

input_training=input_all[index_obs,]
input_testing=input_all[missing_index,]
output_training=as.vector(Y_full[index_obs])
output_testing=as.vector(Y_full[missing_index])

record_time_lagp=system.time(
for(i in 1:1){
    pred_lagp <- aGPsep(input_training,output_training, input_testing, method="nn",verb=0)
}
)
#sd(output_testing)

UB_lagp=pred_lagp$mean+sqrt(pred_lagp$var)*qnorm(0.975)
LB_lagp=pred_lagp$mean+sqrt(pred_lagp$var)*qnorm(0.025)
sqrt(mean((pred_lagp$mean - output_testing)^2))

sum(LB_lagp<=Y_full[missing_index] &UB_lagp>=Y_full[missing_index])/length(missing_index)

mean(UB_lagp-LB_lagp)
record_time_lagp


###FRK



time_date = as.Date(as.yearmon(time))
time_long = rep(time_date, each=dim(Y)[1])
lon_long = rep(Lon,length(Lat)*length(time))
lat_long = rep(rep(Lat,each=length(Lon)), length(time))
refined_data = data.frame(temperature = as.vector(Y_sample), lon = lon_long, lat = lat_long, time = time_long)

#######################
#   Training
#######################
library(FRK)
library(sp)
library(spacetime)
library(zoo)
library(INLA)

# create the long-time-format data frame
time_date = as.Date(as.yearmon(time))
time_long = rep(time_date, each=dim(Y)[1])
lon_long = rep(Lon,length(Lat)*length(time))
lat_long = rep(rep(Lat,each=length(Lon)), length(time))
refined_data = data.frame(temperature = as.vector(Y_sample), lon = lon_long, lat = lat_long, time = time_long)



STObj <- stConstruct(x = refined_data, # dataset
                     space = c("lon","lat"), # spatial fields
                     time="time", # time field
                     interval=TRUE) # time reflects an interval

STObj$std <- 1 # the standard deviation of the measurement error

pred_BAUs <- auto_BAUs(manifold=STplane(), # spatio-temporal process on the plane
                       data=STObj, # data
                       cellsize = c(5,1,1), # BAU cell size
                       type="grid", # grid or hex?
                       # nonconvex_hull = F, # whether to use INLA to generate a non-convex hull
                       convex=-0.1, # parameter for hull construction
                       tunit="months") # time unit

pred_BAUs$fs = 1 # fine-scale variation

##120 spatial basis?
G_spatial <- auto_basis(manifold = plane(), # spatial functions on the plane
                        data=as(STObj,"Spatial"), # remove the temporal dimension
                        nres = 2, # controls the number of basis functions. The most important parameter.
                        # max_basis = 20, # max number of basis functions
                        type = "bisquare", # bisquare basis functions
                        regular = 1, # regular basis functions
                        scale_aperture = 1)  # range from 1 to 1.5

G_temporal <- local_basis(manifold = real_line(), # functions on the real line
                          type = "Gaussian", # Gaussian functions
                          loc = matrix(seq(1,240,by=40)), # locations of functions
                          scale = rep(1,length(seq(1,240,by=40)))) # scales of functions


G <- TensorP(G_spatial,G_temporal) # take the tensor product
G

f <- temperature ~ lon + lat # fixed effects part

record_time_FRK=system.time({
  S <- FRK(f = f,               # formula
           data = list(STObj), # (list of) data
           basis = G,           # basis functions
           BAUs = pred_BAUs,         # BAUs
           n_EM = 500,            # max. no. of EM iterations
           tol = 0.01,
           est_error = FALSE)          # tol. on change in log-likelihood
})
pred_FRK <- predict(S)


#######################
# prediction
#######################

#Y_pred_FRK = matrix(0,1008,240)

#CL_L_FRK = matrix(0,1008,240)

#CL_U_FRK = matrix(0,1008,240)

BAUs_pred_df <- data.frame(pred_FRK) # convert to data frame

BAUs_pred_df$sd_obs <- sqrt(BAUs_pred_df$sd^2 + S@Ve[1,1])
BAUs_pred_df$conflo <- BAUs_pred_df$mu - 1.96*BAUs_pred_df$sd_obs
BAUs_pred_df$confhi <- BAUs_pred_df$mu + 1.96*BAUs_pred_df$sd_obs

real_location_2 = (BAUs_pred_df$lon %in% Lon) & (BAUs_pred_df$lat %in% Lat)
sum(real_location_2)
BAUs_pred_df_part = BAUs_pred_df[real_location_2,]
Y_pred_FRK = matrix(BAUs_pred_df_part$mu,1008,241)
CL_L_FRK = matrix(BAUs_pred_df_part$conflo,1008,241)
CL_U_FRK = matrix(BAUs_pred_df_part$confhi,1008,241)
Y_pred_FRK = Y_pred_FRK[,-1]
CL_L_FRK = CL_L_FRK[,-1]
CL_U_FRK = CL_U_FRK[,-1]

# RMSE
sqrt(mean((Y_pred_FRK[missing_index]-Y_full[missing_index])^2))

#sqrt(mean((Y_sample[missing_index]-Y_full[missing_index])^2))

# CVG
length(which(CL_L_FRK[missing_index]<=Y_full[missing_index] & CL_U_FRK[missing_index]>=Y_full[missing_index]))/length(missing_index)

# INT
mean(CL_U_FRK[missing_index]-CL_L_FRK[missing_index])
