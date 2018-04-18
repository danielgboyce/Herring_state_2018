library(rgdal)
library(sp)
library(maptools)
library(raster)
library(spatstat)
library(gdistance)
library(maps)
library(fossil)
library(lubridate)
library(plyr)
library(RODBC)
datadir<-'C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/herring/tagging'
datadir<-'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/herring/tagging'
setwd(datadir)

#########################################################################
#EXTRACT TABLES FROM ACCESS DATABASE; NEEDS TO RUN ON 32 BIT
channel<-odbcConnectAccess('C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/herring/tagging/herring_tagging_sabs_2003.mdb')
ret<-sqlQuery(channel, paste('select * from HerringTagReturns'))
rel<-sqlQuery(channel, paste('select * from TagReleases'))
names(rel)<-tolower(names(rel))
names(ret)<-tolower(names(ret))

dat<-merge(rel,ret,by=c('tag_number'),all=FALSE)
save(dat,file='herring_tagging_bio.RData')
#########################################################################

setwd(datadir)
load('herring_tagging_bio.RData')
#FIX NAMES TO SHORTEN
names(dat)<-gsub('_','.',names(dat))
names(dat)<-gsub('release','rel',names(dat))
names(dat)<-gsub('return','ret',names(dat))
names(dat)<-gsub('longitude','lon',names(dat))
names(dat)<-gsub('latitude','lat',names(dat))
names(dat)<-gsub('location','loc',names(dat))
names(dat)<-gsub('tag.number','tagnum',names(dat))

#FORMAT LON/LAT TO DECIMAL DEGREES
dat$ret.lon<-dat$ret.lon/100
dat$ret.lon<-ifelse(dat$ret.lon==0,NA,dat$ret.lon)
dg<-as.numeric(substr(dat$ret.lon,1,2))
mn<-as.numeric(substr(dat$ret.lon,3,5))/0.6
dat$ret.lon<-(dg+mn)*-1

dat$ret.lat<-dat$ret.lat/100
dat$ret.lat<-ifelse(dat$ret.lat==0,NA,dat$ret.lat)
dg<-as.numeric(substr(dat$ret.lat,1,2))
mn<-as.numeric(substr(dat$ret.lat,3,5))/0.6
dat$ret.lat<-(dg+mn)

dat$rel.lon<-dat$rel.lon/100
dat$rel.lon<-ifelse(dat$rel.lon==0,NA,dat$rel.lon)
dg<-as.numeric(substr(dat$rel.lon,1,2))
mn<-as.numeric(substr(dat$rel.lon,3,5))/0.6
dat$rel.lon<-(dg+mn)*-1

dat$rel.lat<-dat$rel.lat/100
dat$rel.lat<-ifelse(dat$rel.lat==0,NA,dat$rel.lat)
dg<-as.numeric(substr(dat$rel.lat,1,2))
mn<-as.numeric(substr(dat$rel.lat,3,5))/0.6
dat$rel.lat<-(dg+mn)

#FORMAT DATES
dat$rel.date<-strptime(substr(dat$rel.date,1,10),'%Y-%m-%d')
dat$rel.year<-as.numeric(substr(dat$rel.date,1,4))
dat$rel.month<-as.numeric(substr(dat$rel.date,6,7))
dat$rel.day<-yday(dat$rel.date)
dat$rel.mday<-as.numeric(substr(dat$rel.date,9,10))

dat$ret.date<-strptime(substr(dat$ret.date,1,10),'%Y-%m-%d')
dat$ret.year<-as.numeric(substr(dat$ret.date,1,4))
dat$ret.month<-as.numeric(substr(dat$ret.date,6,7))
dat$ret.day<-yday(dat$ret.date)
dat$ret.mday<-as.numeric(substr(dat$ret.date,9,10))


#CALCUALTE DAYS AT LARGE (DAL)
dat$dal<-as.numeric(difftime(dat$ret.date,dat$rel.date,units='days'))
dat<-subset(dat,dal>=0)

#SET DATES TO NULL ELSE PLYR WON'T WORK
dat$rel.date<-NULL
dat$ret.date<-NULL

#ADD DISTANCE TRAVELLED AND DAYS TRAVELLED
f<-function(d){
d$dist<-deg.dist(d$rel.lon,d$rel.lat,d$ret.lon,d$ret.lat)
return(d)
}
l<-dlply(dat,.(tagnum),.fun=f,.progress='text')
dat<-data.frame(do.call('rbind',l))

#LOOPS THROUGH AND ADDS RETURN LON/LAT IN INSTANCES WHERE MISSING AND RETURN LOCATIONS ARE KNOWN
f<-function(d){

if(is.na(d$ret.lon)==TRUE){

if(d$ret.loc=='Wolves'){   d$ret.lon<--66.73;   d$ret.lat<-44.93}
if(d$ret.loc=='German Bank'){d$ret.lon<--66.4;   d$ret.lat<-43.33}
if(d$ret.loc=='Scots Bay'){ d$ret.lon<--66.39;   d$ret.lat<-45.3}
if(d$ret.loc=='Bliss Islands'){ d$ret.lon<--66.85;   d$ret.lat<-45.02}
if(d$ret.loc=='Chebucto Head'){ d$ret.lon<--63.5;   d$ret.lat<-44.5}
if(d$ret.loc=='Deer Islands'){ d$ret.lon<--68.66;   d$ret.lat<-44.24}
if(d$ret.loc=='Campobello'){ d$ret.lon<--66.92;   d$ret.lat<-44.89}
if(d$ret.loc=='Long Island Shore'){ d$ret.lon<--73.14;   d$ret.lat<-40.79}
if(d$ret.loc=='Beaver Harbour'){ d$ret.lon<--66.74;   d$ret.lat<-45.03}
if(d$ret.loc=='Margaretville'){ d$ret.lon<--65.06;   d$ret.lat<-45.05}
if(d$ret.loc=='Lurcher'){ d$ret.lon<--66.3;   d$ret.lat<-45.73}
if(d$ret.loc %in% c('Grand Manan','North Head, Grand Manan')){d$ret.lon<--66.8;   d$ret.lat<-44.71}

} else { NULL }
return(d)
}
l<-dlply(dat,.(tagnum),.fun=f,.progress='text')
dat<-data.frame(do.call('rbind',l))
dat<-subset(dat,dal>=0 & is.na(ret.lon)==FALSE & is.na(ret.lat)==FALSE& is.na(rel.lon)==FALSE & is.na(rel.lat)==FALSE)



plot(log10(dat$dal),(dat$dist),pch=15)

par(mfrow=c(1,2))
map('world',col='gray',fill=TRUE,xlim=c(-73,-61),ylim=c(40,46))
points(dat$rel.lon,dat$rel.lat,pch=16,col='red')
map('world',col='gray',fill=TRUE,xlim=c(-73,-61),ylim=c(40,46))
points(dat$ret.lon,dat$ret.lat,pch=16,col='green')

######################################################

######################################################
#IMPORT TAGGING DATA FROM M FOWLER (1970'S)
tag<-load('herring_tagging_fowler.RData')
tag<-herring
names(tag)<-tolower(names(tag))


#FIX NAMES TO SHORTEN
names(tag)<-gsub('_','.',names(tag))
names(tag)<-gsub('longitude','rel.lon',names(tag))
names(tag)<-gsub('latitude','rel.lat',names(tag))
names(tag)<-gsub('clong','ret.lon',names(tag))
names(tag)<-gsub('clat','ret.lat',names(tag))

#FORMAT LON/LAT TO DECIMAL DEGREES
tag$ret.lon<-tag$ret.lon/100
tag$ret.lon<-ifelse(tag$ret.lon==0,NA,tag$ret.lon)
dg<-as.numeric(substr(tag$ret.lon,1,2))
mn<-as.numeric(substr(tag$ret.lon,3,5))/0.6
tag$ret.lon<-(dg+mn)*-1

tag$ret.lat<-tag$ret.lat/100
tag$ret.lat<-ifelse(tag$ret.lat==0,NA,tag$ret.lat)
dg<-as.numeric(substr(tag$ret.lat,1,2))
mn<-as.numeric(substr(tag$ret.lat,3,5))/0.6
tag$ret.lat<-(dg+mn)

tag$rel.lon<-tag$rel.lon/100
tag$rel.lon<-ifelse(tag$rel.lon==0,NA,tag$rel.lon)
dg<-as.numeric(substr(tag$rel.lon,1,2))
mn<-as.numeric(substr(tag$rel.lon,3,5))/0.6
tag$rel.lon<-(dg+mn)*-1

tag$rel.lat<-tag$rel.lat/100
tag$rel.lat<-ifelse(tag$rel.lat==0,NA,tag$rel.lat)
dg<-as.numeric(substr(tag$rel.lat,1,2))
mn<-as.numeric(substr(tag$rel.lat,3,5))/0.6
tag$rel.lat<-(dg+mn)

#FORMAT DATES
tag$rel.date<-strptime(gsub(' ','',paste(tag$year,'-',tag$month,'-',tag$day)),'%Y-%m-%d')
tag$rel.year<-as.numeric(substr(tag$rel.date,1,4))
tag$rel.month<-as.numeric(substr(tag$rel.date,6,7))
tag$rel.day<-yday(tag$rel.date)
tag$rel.mday<-as.numeric(substr(tag$rel.date,9,10))

tag$ret.date<-strptime(gsub(' ','',paste(tag$cyear,'-',tag$cmonth,'-',tag$cday)),'%Y-%m-%d')
tag$ret.year<-as.numeric(substr(tag$ret.date,1,4))
tag$ret.month<-as.numeric(substr(tag$ret.date,6,7))
tag$ret.day<-yday(tag$ret.date)
tag$ret.mday<-as.numeric(substr(tag$ret.date,9,10))

tag$dal<-as.numeric(difftime(tag$ret.date,tag$rel.date,units='days'))
tag<-subset(tag,dal>=0 & is.na(ret.lon)==FALSE & is.na(ret.lat)==FALSE& is.na(rel.lon)==FALSE & is.na(rel.lat)==FALSE)

tag$rel.date<-NULL
tag$ret.date<-NULL

#ADD DISTANCE TRAVELLED AND DAYS TRAVELLED
f<-function(d){
    d$dist<-deg.dist(d$rel.lon,d$rel.lat,d$ret.lon,d$ret.lat)
return(d)
}
tag<-ddply(tag,.(tagnum),.fun=f,.progress='text')

plot(log10(tag$dal),(tag$dist),pch=15)
plot(log10(tag$dal),log10(tag$dist),pch=15)


par(mfrow=c(1,2))
map('world',col='gray',fill=TRUE,xlim=c(-73,-57),ylim=c(40,50))
points(tag$rel.lon,tag$rel.lat,pch=16,col='red')
map.axes()
map('world',col='gray',fill=TRUE,xlim=c(-73,-57),ylim=c(40,50))
points(tag$ret.lon,tag$ret.lat,pch=16,col='green')
map.axes()

a<-subset(tag,rel.lat>=43 & rel.lat<=44.2 & rel.lon>=-67 & rel.lon<=-65 & rel.month>=8 & rel.month<=11  & ret.month>=8 & ret.month<=11)
par(mfrow=c(1,2))
map('world',col='gray',fill=TRUE,xlim=c(-73,-57),ylim=c(40,50))
points(a$rel.lon,a$rel.lat,pch=16,col='red')
map.axes()
map('world',col='gray',fill=TRUE,xlim=c(-73,-57),ylim=c(40,50))
points(a$ret.lon,a$ret.lat,pch=16,col='green')
map.axes()


dat$db<-'sabs'
tag$db<-'bio'
dat$rel.date<-NULL
dat$ret.date<-NULL
data<-rbind.fill(dat,tag)

data<-data[,!(names(data) %in% c('set.number','rel.comments','date.entered','ret.vessel.name','ret.gear.type','ret.fish.plant','ret.name','ret.address','ret.comments','species.res','species.com','program','gear','gear.com','effort.cnt','effort.min','cruise','time','stn','nafo','depth','tempdep','tagtype','reward','special','mtags','length','lentype','weight','sex','kept','cspecies.res','cspecies.com','cyear','cmonth','cday','cnafo','areadesc','locdet','cgear','cgear.com','cvessel','nation','flag','cdepth','ctemp','ctempdep','clength','clentype','cweight','cwgttype','csex','age','amat','found','measure','comments','temp','month','year','month','day','dal2'))]


#save(data,file='herring_tagging_data_spera.RData')
load('herring_tagging_data_spera.RData')

#LOAD EXCLUSIVE ECONOMIC ZONE FOR CANADA/US
eez.us<-readOGR('N:/data/stock_assessments/data/fishing_shapefiles/eez_oceans_intersect','EEZ_IHO_union_v2')
eez.us<-spTransform(eez.us,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
eez<-subset(eez.us,EEZ%in% c('United States Exclusive Economic Zone','Canadian Exclusive Economic Zone'))
eezus<-subset(eez.us,EEZ%in% c('United States Exclusive Economic Zone'))

plot(eez,col=ifelse(eez$EEZ=='Canadian Exclusive Economic Zone','dodgerblue3','firebrick3'),xlim=c(-75,-60),ylim=c(40,50),axes=TRUE)
points(data$ret.lon,data$ret.lat,pch=16,col=alpha('green',.4),cex=.5)
points(data$rel.lon,data$rel.lat,pch=16,col=alpha('green',.4),cex=.5)


crds<-SpatialPoints(data.frame(lon=data$ret.lon,lat=data$ret.lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,eezus)
data$dum<-gint[,1]
data$dum<-ifelse(is.na(data$dum)==FALSE,'us','can')

#GET TAGS DURING SPAWNING PERIOD IN SW NS
plot(0,0,xlim=c(-71,-58),ylim=c(42,49))
plot(cst,col='lightblue',add=TRUE)
points(data$rel.lon,data$rel.lat,pch=16,col='red')
a<-subset(data,rel.mday %in% c(8,9,10))
points(a$rel.lon,a$rel.lat,pch=16,col='green')

lon1<--67
lon2<--65.5
lat1<-43.2
lat2<-44.25
aa<-subset(a,rel.lon>lon1 & rel.lon<lon2 & rel.lat>lat1 & rel.lat<lat2)
points(aa$ret.lon,aa$ret.lat,pch=16,col='blue')

library(scales)
par(mfrow=c(3,3))
mnth<-c(11,12,1,2,3,4,5,6,7)
for(i in 1:length(mnth)){
d<-subset(aa,ret.month==mnth[i])
plot(0,0,xlim=c(-71,-58),ylim=c(42,49))
plot(cst,col='lightblue',add=TRUE)
points(d$ret.lon,d$ret.lat,pch=16,col=alpha('red3',.4),cex=2)
legend('top',paste('month=',mnth[i]))
}

plot(0,0,xlim=c(-71,-58),ylim=c(42,49))
plot(cst,col='lightblue',add=TRUE)
points(data$rel.lon,data$rel.lat,pch=16,col='red')
#a<-subset(data,rel.mday %in% c(8,9,10))
#points(a$rel.lon,a$rel.lat,pch=16,col='green')

#GET TAGS RELEASED IN PASSAMAQUODDY REGION
lon1<--67.5
lon2<--66.5
lat1<-44.4
lat2<-45.25
pas<-subset(data,rel.lon>lon1 & rel.lon<lon2 & rel.lat>lat1 & rel.lat<lat2)

pasr<-subset(pas,ret.lon>lon1 & ret.lon<lon2 & ret.lat>lat1 & ret.lat<lat2)
pasnr<-subset(pas,(ret.lon<lon1 | ret.lon>lon2) & (ret.lat<lat1 | ret.lat>lat2))
points(pasnr$ret.lon,pasnr$ret.lat,pch=16,col='green')

par(mfrow=c(3,3))
mnth<-c(11,12,1,2,3,4,5,6,7)
for(i in 1:length(mnth)){
d<-subset(pas,ret.month==mnth[i])
plot(0,0,xlim=c(-71,-58),ylim=c(42,49))
plot(cst,col='lightblue',add=TRUE)
points(d$ret.lon,d$ret.lat,pch=16,col=alpha('red3',.4),cex=2)
legend('top',paste('month=',mnth[i]))
}


#OVERLAY LAND MASSES
cst<-readShapePoly('N:/data/shapefiles/naturalearthdata_ne_50m_ocean_poly/ne_50m_ocean.shp',proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cst<-crop(cst,extent(-75,-55,35,55))
proj4string(cst)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Convert the lake's boundary to a raster, with values of 1 for
## cells within the lake and values of 0 for cells on land
#bound <- owin(poly=data.frame(x=ex.poly$x, y=ex.poly$y))
#Simple example of a polygon and points.
Poly<-cst
R <- raster(extent(Poly), nrow=300,  ncol=300) ## 2nd to RasterLaye
RR <- rasterize(Poly, R)                       ## ...
RR[is.na(RR)]<-0  ## Set cells on land to "0"

## gdistance requires that you 1st prepare a sparse "transition matrix"
## whose values give the "conductance" of movement between pairs of
## adjacent and next-to-adjacent cells (when using directions=16)
tr1 <- transition(RR, transitionFunction=mean, directions=16)
tr1 <- geoCorrection(tr1,type="c")


#a<-subset(tag,dist>100)
#a<-a[1:10,]
a<-data
crds1<-SpatialPoints(data.frame(lon=a$rel.lon,lat=a$rel.lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
crds2<-SpatialPoints(data.frame(lon=a$ret.lon,lat=a$ret.lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

b<-costDistance(tr1,crds1,crds2)[1,]
a$tdist<-b/1000

plot(a$dist,a$tdist)
abline(a=0,b=1)
m<-subset(a,tdist<dist)
points(m$dist,m$tdist,pch=16,col='red')


data<-a

a<-subset(a,tdist==Inf)

#GET VECTOR PATHWAYS BETWEEN ALL RELEASE AND RETURN TAGS
sl<-list()
for(i in 1:dim(data)[1]){
print(i)
sl[[i]]<-shortestPath(tr1, crds1[i],crds2[i], output="SpatialLines")
}
save(sl,file='herring_tagging_vectors.RData')

#PLOT ALL VECTOR PATHWAYS
pdf('herring_tracks.pdf',height=9,width=10)
#plot(RR,col=c('white','lightblue'),xlim=c(-74,-57),ylim=c(40,50),las=1)
map('world',fill=TRUE,col='white',xlim=c(-74,-57),ylim=c(40,50),las=1)
lapply(sl, function(X) plot(X, col=alpha("red",.1), add=TRUE, lwd=2))
plot(crds1, pch=16, col=alpha("gold3",.5), cex=.1, add=TRUE)
plot(crds2, pch=16, col=alpha("dodgerblue",.5), cex=.1, add=TRUE)
dev.off()

library(scales)

## View the selected paths
plot(RR)
plot(crds1, pch=16, col="gold", cex=1.5, add=TRUE)
plot(crds2, pch=16, col="red", cex=1.5, add=TRUE)
SL12 <- shortestPath(tr1, crds1,crds2, output="SpatialLines")
SL13 <- shortestPath(tr1, crds1,crds2, output="SpatialLines")
SL23 <- shortestPath(tr1, crds1,crds2, output="SpatialLines")
lapply(list(SL12, SL13, SL23), function(X) plot(X, col="red", add=TRUE, lwd=2))



m$delta<-m$tdist-m$dist
m2<-subset(m,delta< -377)
m2<-m2[1,]

crds1<-SpatialPoints(data.frame(lon=m2$rel.lon,lat=m2$rel.lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
crds2<-SpatialPoints(data.frame(lon=m2$ret.lon,lat=m2$ret.lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

mm<-costDistance(tr1,crds1,crds2)[1,]/1000

plot(RR)
plot(crds1, pch=16, col="gold", cex=1.5, add=TRUE)
plot(crds2, pch=16, col="red", cex=1.5, add=TRUE)
SL12 <- shortestPath(tr1, crds1,crds2, output="SpatialLines")
SL13 <- shortestPath(tr1, crds1,crds2, output="SpatialLines")
SL23 <- shortestPath(tr1, crds1,crds2, output="SpatialLines")
lapply(list(SL12, SL13, SL23), function(X) plot(X, col="red", add=TRUE, lwd=2))
deg.dist(m2$rel.lon,m2$rel.lat,m2$ret.lon,m2$ret.lat)
