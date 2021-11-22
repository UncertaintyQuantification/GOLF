da1=scan("NOAA_gridded_temperature/NOAAGlobalTemp.gridded.v4.0.1.201901.asc")
length(da1)
da1[1:20]

dim(da1)
tm1=seq(1,4329386-1, by=2594)
tm2=seq(2,4329386, by=2594)
length(tm1)
length(tm2)
mm1=da1[tm1] #Extract months
yy1=da1[tm2] #Extract years
rw1<-paste(yy1, sep="-", mm1) #Combine YYYY with MM
head(mm1)
head(yy1)
rw1
tm3=cbind(tm1,tm2)
tm4=as.vector(t(tm3))
head(tm4)
#[1] 1 2 2595 2596 5189 5190
da2<-da1[-tm4] #Remote the months and years data from the scanned data
length(da2)/(36*72)
#
da3<-matrix(da2,ncol=1669) #Generate the space-time data
#2592 (=36*72) rows and 1669 months 


colnames(da3)<-rw1
lat1=seq(-87.5, 87.5, length=36)
lon1=seq(2.5, 357.5, length=72)
LAT=rep(lat1, each=72)
LON=rep(lon1,36)
gpcpst=cbind(LAT, LON, da3)
head(gpcpst)
dim(gpcpst)

##first column of gpcpst is latitude and second column is the longitude 
##each column of the rest of the column is a monthly temperature anomaly at a certain grid
##each row is a grid for all the later month

rw1[1057]
rw1[1634]

library(maps)#Install maps package if not done before
Lat= seq(-87.5, 87.5, length=36)
Lon=seq(2.5, 357.5, length=72)
mapmat=matrix(gpcpst[,1636],nrow=72)
mapmat[which(mapmat==-999.9)]=NA
#column 1634 corresponding to Dec 2015
#Covert the vector into a lon-lat matrix for R map plotting
#mapmat=pmax(pmin(mapmat,6),-6)
#Matrix flipping is not needed since the data go from 2.5 to 375.5
plot.new()
#par(mar=c(4,5,3,0))
int=seq(-10,10,length.out=81)
rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
##this is a typo in the original code, and it should be this one
#mapmat= mapmat[,seq(length(mapmat[1,]),1)]
#mapmat= t(cbind(t(mapmat[37:72,]),t(mapmat[1:36,])))

filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="NOAAGlobalTemp Anomalies Feb 2016 [deg C]",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})

##let me create the yearly data from 1969 to 2018
rw1[1069]
rw1[1668]


matrix_1969_2018=gpcpst[,1071:1670]  ###note here there are two difference as the first two column are coordinate

n_missing_month=rep(0,12*50)
matrix_1969_2018[which(matrix_1969_2018==-999.9)]=NA
for(i_month in 1:600){
  n_missing_month[i_month]=sum(is.na(matrix_1969_2018[,i_month]))
}

n_missing_month/(32*76)

NOAA_gridded_data=list()

NOAA_gridded_data$temperature=matrix_1969_2018
NOAA_gridded_data$lat=gpcpst[,1]
NOAA_gridded_data$lon=gpcpst[,2]

NOAA_gridded_data$months=rw1[1069:1668]


saveRDS(NOAA_gridded_data,file="NOAA_gridded_data_1969_2018.rds")



##let me create the yearly data from 1989 to 2018
rw1[1309]
rw1[1668]
matrix_1989_2018=gpcpst[,1311:1670]  ###note here there are two difference as the first two column are coordinate


n_missing_month=rep(0,12*30)
matrix_1989_2018[which(matrix_1989_2018==-999.9)]=NA
for(i_month in 1:360){
  n_missing_month[i_month]=sum(is.na(matrix_1989_2018[,i_month]))
}

n_missing_month/(32*76)

NOAA_gridded_data=list()

NOAA_gridded_data$temperature=matrix_1989_2018
NOAA_gridded_data$lat=gpcpst[,1]
NOAA_gridded_data$lon=gpcpst[,2]

NOAA_gridded_data$months=rw1[1309:1668]


saveRDS(NOAA_gridded_data,file="NOAA_gridded_data_1989_2018.rds")

