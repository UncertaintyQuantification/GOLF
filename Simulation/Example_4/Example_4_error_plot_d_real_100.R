num_simulation=100

record_nonseparable_RMSE=matrix(0,num_simulation,7)
record_separable_RMSE=matrix(0,num_simulation,7)
record_nonseparable_truth_RMSE=matrix(0,num_simulation,7)
record_separable_truth_RMSE=matrix(0,num_simulation,7)

record_nonseparable_interval=matrix(0,num_simulation,7)
record_separable_interval=matrix(0,num_simulation,7)

record_nonseparable_length=matrix(0,num_simulation,7)
record_separable_length=matrix(0,num_simulation,7)

j=1
i=5
load(paste0("nonseparable_prediction_d_real_100_d_",i,".Rdata"))
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
  load(paste0("nonseparable_prediction_d_real_100_d_",i,".Rdata"))
  record_nonseparable_RMSE[,j+1]=record_nonseparable_error[,1]
  record_nonseparable_truth_RMSE[,j+1]=record_nonseparable_error[,2]
  record_nonseparable_interval[,j+1]=record_nonseparable_error[,3]
  record_nonseparable_length[,j+1]=record_nonseparable_error[,4]
  
  record_separable_RMSE[,j+1]=record_separable_error[,1]
  record_separable_truth_RMSE[,j+1]=record_separable_error[,2]
  record_separable_interval[,j+1]=record_separable_error[,3]
  record_separable_length[,j+1]=record_separable_error[,4]
}


j=7
i=100
load(paste0("nonseparable_prediction_d_real_100_d_",i,".Rdata"))
record_nonseparable_RMSE[,j]=record_nonseparable_error[,1]
record_nonseparable_truth_RMSE[,j]=record_nonseparable_error[,2]
record_nonseparable_interval[,j]=record_nonseparable_error[,3]
record_nonseparable_length[,j]=record_nonseparable_error[,4]

record_separable_RMSE[,j]=record_separable_error[,1]
record_separable_truth_RMSE[,j]=record_separable_error[,2]
record_separable_interval[,j]=record_separable_error[,3]
record_separable_length[,j]=record_separable_error[,4]


d=c(5,10,20,30,40,50,100)
pdf("RMSE_nonseparable_data_d_real_100.pdf",height=4,width=4)
plot(d,sqrt(colMeans(record_nonseparable_RMSE^2)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,mgp=c(2.5,1,0),ylim=c(0.05,0.70),ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_RMSE^2)),type='b',col='red',pch=17)

lines(d,sqrt(colMeans(record_nonseparable_truth_RMSE^2)),type='b',col='blue',pch=20,lty=4,ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_truth_RMSE^2)),type='b',col='red',pch=17,lty=4)

legend("topright", legend=c("nonseparable pred data", "separable pred data", "nonseparable pred mean", "separable pred mean"),
       pch=c(20,17,20,17),lty=c(1,1,4,4),col=c('blue','red','blue','red'))
dev.off()



pdf("coverage_nonseparable_data_d_real_100.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_interval)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.5,0.95),ylab=expression(P[CI](95~'%')))
lines(d,(colMeans(record_separable_interval)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("bottomright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()


pdf("length_nonseparable_data_d_real_100.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_length)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.45,0.65),ylab=expression(L[CI](95~'%')))
lines(d,(colMeans(record_separable_length)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("topright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()


################################################################################

record_nonseparable_RMSE2=matrix(0,num_simulation,7)
record_separable_RMSE2=matrix(0,num_simulation,7)
record_nonseparable_truth_RMSE2=matrix(0,num_simulation,7)
record_separable_truth_RMSE2=matrix(0,num_simulation,7)

record_nonseparable_interval2=matrix(0,num_simulation,7)
record_separable_interval2=matrix(0,num_simulation,7)

record_nonseparable_length2=matrix(0,num_simulation,7)
record_separable_length2=matrix(0,num_simulation,7)

j=1
i=5*j
load(paste0("separable_prediction_d_real_100_d_",i,".Rdata"))
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
  load(paste0("separable_prediction_d_real_100_d_",i,".Rdata"))
  record_nonseparable_RMSE2[,j+1]=record_nonseparable_error[,1]
  record_nonseparable_truth_RMSE2[,j+1]=record_nonseparable_error[,2]
  record_nonseparable_interval2[,j+1]=record_nonseparable_error[,3]
  record_nonseparable_length2[,j+1]=record_nonseparable_error[,4]
  
  record_separable_RMSE2[,j+1]=record_separable_error[,1]
  record_separable_truth_RMSE2[,j+1]=record_separable_error[,2]
  record_separable_interval2[,j+1]=record_separable_error[,3]
  record_separable_length2[,j+1]=record_separable_error[,4]
}

j=7
i=100
load(paste0("separable_prediction_d_real_100_d_",i,".Rdata"))
record_nonseparable_RMSE2[,j]=record_nonseparable_error[,1]
record_nonseparable_truth_RMSE2[,j]=record_nonseparable_error[,2]
record_nonseparable_interval2[,j]=record_nonseparable_error[,3]
record_nonseparable_length2[,j]=record_nonseparable_error[,4]

record_separable_RMSE2[,j]=record_separable_error[,1]
record_separable_truth_RMSE2[,j]=record_separable_error[,2]
record_separable_interval2[,j]=record_separable_error[,3]
record_separable_length2[,j]=record_separable_error[,4]

d=c(5,10,20,30,40,50,100)
pdf("RMSE_separable_data_d_real_100.pdf",height=4,width=4)
plot(d,sqrt(colMeans(record_nonseparable_RMSE2^2)),mgp=c(2.5,1,0),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.05,0.7),ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_RMSE2^2)),type='b',col='red',pch=17)

lines(d,sqrt(colMeans(record_nonseparable_truth_RMSE2^2)),type='b',col='blue',pch=20,lty=4,ylab='RMSE')
lines(d,sqrt(colMeans(record_separable_truth_RMSE2^2)),type='b',col='red',pch=17,lty=4)

legend("topright", legend=c("nonseparable pred data", "separable pred data", "nonseparable pred mean", "separable pred mean"),
       pch=c(20,17,20,17),lty=c(1,1,4,4),col=c('blue','red','blue','red'))
dev.off()

pdf("coverage_separable_data_100.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_interval2)),type='b',col='blue',mgp=c(2.5,1,0),pch=20,ylim=c(0.5,0.95),ylab=expression(P[CI](95~'%')))
lines(d,(colMeans(record_separable_interval2)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("bottomright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()


pdf("length_separable_data_100.pdf",height=4,width=4)
plot(d,(colMeans(record_nonseparable_length2)),type='b',col='blue',mgp=c(2.5,1,0),
     pch=20,ylim=c(0.45,0.65),ylab=expression(L[CI](95~'%')))
lines(d,(colMeans(record_separable_length2)),type='b',col='red',pch=17)
abline(a=0.95,b=0,lty=2)
legend("topright", legend=c("nonseparable", "separable"),
       pch=c(20,17),lty=c(1,1),col=c('blue','red'))
dev.off()

