library(ggplot2)
library(plot3D)
load("Data_Heaton/AllSatelliteTemps.RData")



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

# temp_lat=colMeans(Truetemp_mat,na.rm = T)
# plot(lat_unique,temp_lat)
# 
# temp_lon=rowMeans(Truetemp_mat,na.rm = T)
# plot(lon_unique,temp_lon)



##let me makes the test est
#this seems to be right for plot
index_lon=1:500
index_lat=300:1
#index_lon=1:200
#index_lat=1:150

Masktemp_mat_sub=Masktemp_mat[index_lon,index_lat]
Truetemp_mat_sub=Truetemp_mat[index_lon,index_lat]

lon_unique_sub=lon_unique[index_lon]
lat_unique_sub=lat_unique[index_lat]

# image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub)
# image2D(Truetemp_mat_sub,lon_unique_sub,lat_unique_sub)

#matplot(Truetemp_mat_sub)

#plot(Truetemp_mat_sub[,3])
#plot(Truetemp_mat_sub[,4])


n1=length(lat_unique_sub)
n2=length(lon_unique_sub)

input1=lat_unique_sub ##row
input2=lon_unique_sub##column

input_2_vec=matrix(t(matrix(rep(input2,n1),n2,n1)),n2*n1)

input=cbind(rep(input1,n2),input_2_vec)

#Y_transpose_full_ori=matrix(NA,n2,n1)

#grid_with_full_obs=which(is.na(Truetemp_mat_sub)==F)

#Y_transpose_full_ori=Truetemp_mat_sub

Y_full=t(Truetemp_mat_sub)
grid_with_full_obs_Y_full=which(is.na(Y_full)==F)

Y=t(Masktemp_mat_sub)
grid_with_full_obs_Y=which(is.na(Y)==F)



record_coefficient1=rep(0,500)
record_bound1=matrix(0,2,500)
for(i in 1:500){
  # m=lm(Y[,i]~seq(1,300)) 
  m=lm(Y[,i]~input1-1)
  record_coefficient1[i]=m$coefficients
  record_bound1[,i]=confint(m, level=0.95)
}

pdf("plot_trend_b_1.pdf",height=4,width=5)

plot(input2,record_coefficient1,type='p',xlab='longitude',ylab=expression(hat(b)[1]),
     xlim=c(min(input2)+0.165,max(input2)-0.165),
     mgp=c(2,0.5,0),pch=20,cex=1)
dev.off()
#polygon( c(input2,rev(input2)),c(record_bound1[1,],
#                                               rev(record_bound1[2,])),col = "grey80", border = F)
#lines(input2,record_coefficient1,type='l')




record_coefficient2=rep(0,300)
record_bound2=matrix(0,2,300)
for(i in 1:300){
  
  #m_full=lm(Y_full[,i]~seq(1,300)) 
  m=lm(Y[i,]~input2-1)
  record_coefficient2[i]=m$coefficients
  
  record_bound2[,i]=confint(m, level=0.95)
}


pdf("plot_trend_b_2.pdf",height=4,width=5)
plot(input1,record_coefficient2,type='p',xlab='latitude',ylab=expression(hat(b)[2]),
     xlim=c(min(input1)+0.095,max(input1)-0.095),
     mgp=c(2,0.5,0),pch=20,cex=1)
dev.off()



# ##
# pdf("plot_trend_b_2.pdf",height=3,width=6)
# plot((input1),(record_coefficient2),type='l',axes=F,xlab='',ylab='',
#      xlim=c(min(input1)+0.095,max(input1)-0.095),
#      mgp=c(2.5,0.5,0))
# box()
# #ylab=expression(hat(b)[2]),
# axis(side = 3, las = 1,padj=.5)
# mtext('latitude',side=3,padj=-3)
# axis(side = 2, las = 1,hadj=.9)
# #text(side=2,padj=-2.3, expression(hat(b)[2]))
# 
# mtext(expression(hat(b)[2]),las = 1,side=2,padj=0,adj=4.8)
# 
# dev.off()

#image2D((input2%*%t(record_coefficient2)),lon_unique_sub,lat_unique_sub,ylab='',xlab='',mgp=c(3,0.5,0))




png("obs_temp.png",height=4,width=5.5,units="in",res=200)
image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
        main='Observed temperature',zlim=c(24,56))
dev.off()

png("full_temp.png",height=4,width=5.5,units="in",res=200)
image2D(Truetemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
        main='Full temperature',zlim=c(24,56))
dev.off()


# pdf("obs_temp.pdf",height=4,width=5.5)
# #image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',main='Observed temperature',
# #        mgp=c(2,0.5,0))
# image2D(Masktemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',
#         mgp=c(2,0.5,0),zlim=c(24,56))
# 
# dev.off()


# pdf("full_temp.pdf",height=4,width=5.5)
# image2D(Truetemp_mat_sub,lon_unique_sub,lat_unique_sub,ylab='latitude',xlab='longitude',mgp=c(2,0.5,0),
#         main='Full temperature',zlim=c(24,56))
# dev.off()





