library(ggplot2)
library(gridExtra )
library(ggpubr)

load("Data_Heaton/AllSatelliteTemps.RData")


##################
# Figure 1
##################

dat = all.sat.temps

range(dat$Lon)

rotate <- function(x) t(apply(x, 2, rev))

Masktemp_mat=rotate(t(matrix(dat[,3],500,300)))
Truetemp_mat=rotate(t(matrix(dat[,4],500,300)))

## get the lattitude and longitude
lat_unique = unique(dat$Lat) # from 34.29519 to 37.06811
lat_unique
lon_unique = unique(dat$Lon ) # from -95.91153 to -91.28381
lon_unique
xy.grid = expand.grid(lat_unique,lon_unique)

##let me makes a smaller test est
index_lon=1:500
index_lat=1:300

Masktemp_mat_sub=Masktemp_mat[index_lon,index_lat]
Truetemp_mat_sub=Truetemp_mat[index_lon,index_lat]

lon_unique_sub=lon_unique[index_lon]
lat_unique_sub=lat_unique[index_lat]

n1=length(lat_unique_sub)
n2=length(lon_unique_sub)

input1=lat_unique_sub ##row
input2=lon_unique_sub##column


Y=t(Masktemp_mat_sub)

######
par(mfrow=c(1,1))

record_coefficient=rep(0,500)
for(i in 1:500){
  
  m=lm(Y[,i]~input1-1)
  record_coefficient[i]=m$coefficients[1]
}

lon_reg_df = data.frame(x = input2, y = record_coefficient)

plot_1 = ggplot(lon_reg_df) +
  geom_point(aes(x = x, y = y), size = 0.5)+
  xlab("longitude")+
  ylab(expression(hat(b)[1]))


record_coefficient2=rep(0,300)
for(i in 1:300){  
  m=lm(Y[i,]~input2-1)
  record_coefficient2[i]=m$coefficients[1]
}

lat_reg_df = data.frame(x = input1, y = record_coefficient2)

plot_2 = ggplot(lat_reg_df) +
  geom_point(aes(x = x, y = y), size = 0.5)+
  xlab("latitude")+
  ylab(expression(hat(b)[2]))

ggarrange(plot_1, plot_2, nrow=1, ncol=2)
