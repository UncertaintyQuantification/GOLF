#colMeans(record_nonseparable_error)
#colMeans(record_separable_error)
num_simulation=100

record_nonseparable_RMSE=matrix(0,num_simulation,6)
record_separable_RMSE=matrix(0,num_simulation,6)
record_nonseparable_truth_RMSE=matrix(0,num_simulation,6)
record_separable_truth_RMSE=matrix(0,num_simulation,6)

record_nonseparable_interval=matrix(0,num_simulation,6)
record_separable_interval=matrix(0,num_simulation,6)

record_nonseparable_length=matrix(0,num_simulation,6)
record_separable_length=matrix(0,num_simulation,6)

j=1
i=5*j
load(paste0("nonseparable_prediction_d_real_30_d_",i,".Rdata"))
record_nonseparable_RMSE[,j]=record_nonseparable_error[,1]
record_nonseparable_truth_RMSE[,j]=record_nonseparable_error[,2]
record_nonseparable_interval[,j]=record_nonseparable_error[,3]
record_nonseparable_length[,j]=record_nonseparable_error[,4]

record_separable_RMSE[,j]=record_separable_error[,1]
record_separable_truth_RMSE[,j]=record_separable_error[,2]
record_separable_interval[,j]=record_separable_error[,3]
record_separable_length[,j]=record_separable_error[,4]

for(j in 1:5){
   i=10*j
   load(paste0("nonseparable_prediction_d_real_30_d_",i,".Rdata"))
   record_nonseparable_RMSE[,j+1]=record_nonseparable_error[,1]
   record_nonseparable_truth_RMSE[,j+1]=record_nonseparable_error[,2]
   record_nonseparable_interval[,j+1]=record_nonseparable_error[,3]
   record_nonseparable_length[,j+1]=record_nonseparable_error[,4]
   
   record_separable_RMSE[,j+1]=record_separable_error[,1]
   record_separable_truth_RMSE[,j+1]=record_separable_error[,2]
   record_separable_interval[,j+1]=record_separable_error[,3]
   record_separable_length[,j+1]=record_separable_error[,4]
}

d=c(5,10,20,30,40,50)
pdf("RMSE_nonseparable_data.pdf",height=4,width=4)
plot(d,sqrt(colMeans(record_nonseparable_RMSE^2)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,mgp=c(2.5,1,0),ylim=c(0.05,0.70),ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_RMSE^2)),type='b',col='red',pch=17)

lines(d,sqrt(colMeans(record_nonseparable_truth_RMSE^2)),type='b',col='blue',pch=20,lty=4,ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_truth_RMSE^2)),type='b',col='red',pch=17,lty=4)

legend("topright", legend=c("nonseparable pred data", "separable pred data", "nonseparable pred mean", "separable pred mean"),
      pch=c(20,17,20,17),lty=c(1,1,4,4),col=c('blue','red','blue','red'))
dev.off()

pdf("coverage_nonseparable_data.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_interval)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.5,0.95),ylab=expression(P[CI](95~'%')))
lines(d,(colMeans(record_separable_interval)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("bottomright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()


pdf("length_nonseparable_data.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_length)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.45,0.65),ylab=expression(L[CI](95~'%')))
lines(d,(colMeans(record_separable_length)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("topright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()

################################################################################

record_nonseparable_RMSE2=matrix(0,num_simulation,6)
record_separable_RMSE2=matrix(0,num_simulation,6)
record_nonseparable_truth_RMSE2=matrix(0,num_simulation,6)
record_separable_truth_RMSE2=matrix(0,num_simulation,6)

record_nonseparable_interval2=matrix(0,num_simulation,6)
record_separable_interval2=matrix(0,num_simulation,6)

record_nonseparable_length2=matrix(0,num_simulation,6)
record_separable_length2=matrix(0,num_simulation,6)

j=1
i=5*j
load(paste0("separable_prediction_d_real_30_d_",i,".Rdata"))
record_nonseparable_RMSE2[,j]=record_nonseparable_error[,1]
record_nonseparable_truth_RMSE2[,j]=record_nonseparable_error[,2]
record_nonseparable_interval2[,j]=record_nonseparable_error[,3]
record_nonseparable_length2[,j]=record_nonseparable_error[,4]

record_separable_RMSE2[,j]=record_separable_error[,1]
record_separable_truth_RMSE2[,j]=record_separable_error[,2]
record_separable_interval2[,j]=record_separable_error[,3]
record_separable_length2[,j]=record_separable_error[,4]

for(j in 1:5){
  i=10*j
  load(paste0("separable_prediction_d_real_30_d_",i,".Rdata"))
  record_nonseparable_RMSE2[,j+1]=record_nonseparable_error[,1]
  record_nonseparable_truth_RMSE2[,j+1]=record_nonseparable_error[,2]
  record_nonseparable_interval2[,j+1]=record_nonseparable_error[,3]
  record_nonseparable_length2[,j+1]=record_nonseparable_error[,4]
  
  record_separable_RMSE2[,j+1]=record_separable_error[,1]
  record_separable_truth_RMSE2[,j+1]=record_separable_error[,2]
  record_separable_interval2[,j+1]=record_separable_error[,3]
  record_separable_length2[,j+1]=record_separable_error[,4]
}

d=c(5,10,20,30,40,50)
pdf("RMSE_separable_data.pdf",height=4,width=4)
plot(d,sqrt(colMeans(record_nonseparable_RMSE2^2)),mgp=c(2.5,1,0),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.05,0.7),ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_RMSE2^2)),type='b',col='red',pch=17)

lines(d,sqrt(colMeans(record_nonseparable_truth_RMSE2^2)),type='b',col='blue',pch=20,lty=4,ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_truth_RMSE2^2)),type='b',col='red',pch=17,lty=4)

legend("topright", legend=c("nonseparable pred data", "separable pred data", "nonseparable pred mean", "separable pred mean"),
       pch=c(20,17,20,17),lty=c(1,1,4,4),col=c('blue','red','blue','red'))
dev.off()

pdf("coverage_separable_data.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_interval2)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.5,0.95),ylab=expression(P[CI](95~'%')))
lines(d,(colMeans(record_separable_interval2)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("bottomright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()


pdf("length_separable_data.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_length2)),type='b',col='blue',mgp=c(2.5,1,0),
     pch=20,ylim=c(0.45,0.65),ylab=expression(L[CI](95~'%')))
lines(d,(colMeans(record_separable_length2)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("topright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()


###plot a figure 
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
#library(ggplot2)
library(plot3D)
#load("Data_Heaton/AllSatelliteTemps.RData")
library(RobustGaSP)

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
   
   ##cholsky
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

i_plot=1
Y_full=matrix(nonseparable_full_data_all[i_plot,],n1,n2)
Y_truth=matrix(nonseparable_truth_data_all[i_plot,],n1,n2)
Y=matrix(NA,n1,n2)
Y[-missing_index]=Y_full[-missing_index]


image2D(Y_truth,xlab=expression(x[1]),ylab=expression(x[2]))
image2D(Y_full)


load(paste0("nonseparable_prediction_d_real_30_d_",30,".Rdata"))

Y_pred=matrix(return_list_nonseparable_record[[i_plot]][[1]],n1,n2)

image2D(Y_pred)

#sqrt(mean( (Y_pred[missing_index]-Y_full[missing_index])^2))

#sqrt(mean( (Y_pred[missing_index]-Y_truth[missing_index])^2))

#png("pred_temp_golf.png",height=4,width=5.5,units="in",res=200)
max(Y,Y_full,Y_pred,na.rm=T)
min(Y,Y_full,Y_pred,na.rm=T)

png("obs_nonseparable_d_30.png",height=4,width=4,units="in",res=200)
image2D(t(Y),xlab=expression(x[2]),ylab=expression(x[1]),mgp=c(2,0.5,0),
        main='Simulated observations',zlim=c(-1.3,1.7))
dev.off()
png("truth_nonseparable_d_30.png",height=4,width=4,units="in",res=200)
image2D(t(Y_truth),xlab=expression(x[2]),ylab=expression(x[1]),mgp=c(2,0.5,0),
        main='Simulated mean',zlim=c(-1.3,1.7))
dev.off()

png("pred_golf_nonseparable_d_30.png",height=4,width=4,units="in",res=200)

image2D(t(Y_pred),xlab=expression(x[2]),ylab=expression(x[1]),mgp=c(2,0.5,0),
        main='GOLF',zlim=c(-1.3,1.7))
dev.off()




