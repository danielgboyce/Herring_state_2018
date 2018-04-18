data<-data.frame(lon=rnorm(50,mean=-66.24,sd=.1),
                 lat=rnorm(50,mean=43.2,sd=.1))
setwd("c:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/webdrogue")
write.table(data,'data.txt',row.names=FALSE,col.names=FALSE,sep='\t')



data<-data.frame(lon=rnorm(50,mean=-65.2,sd=.05),
                 lat=rnorm(50,mean=45.2,sd=.05))
setwd("c:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/webdrogue")
write.table(data,'data2.txt',row.names=FALSE,col.names=FALSE,sep='\t')



library(maps)
library(sp)
library(rgdal)
library(maptools)
library(raster)
cst<-readShapePoly('N:/data/shapefiles/naturalearthdata_ne_50m_ocean_poly/ne_50m_ocean.shp',proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cst<-crop(cst,extent(-70,-60,40,50))
plot(cst,col='lightblue')

lon1<--70
lon2<--63
lat1<-43
lat2<-46
map('world',xlim=c(-70,-63),ylim=c(43,46),fill=TRUE,col='gray')
map.axes()
dat<-expand.grid(lon=seq(-68,-67,length.out=100),
                 lat=seq(43,44,length.out=100),
                 year=seq(1950,1951,1))

dat<-expand.grid(lon=seq(lon1,lon2,length.out=100),
                 lat=seq(lat1,lat2,length.out=100),
                 year=seq(1950,1951,1))

dat$day<-200
dat$hour<-12
dat$minute<-10
dat$second<-10

crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,cst)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)

dat<-dat[,!((names(dat) %in% c('dum')))]
points(dat$lon,dat$lat,pch=16,col='red')

setwd("c:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/webdrogue")
write.table(dat,'webtidedata.txt',row.names=FALSE,col.names=FALSE,sep='\t')
