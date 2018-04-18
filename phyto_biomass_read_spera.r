library(raster)
library(maps)
library(plyr)
library(lubridate)
library(sp)
library(rgdal)
library(maptools)

#SPECIFY DIRECTORIES
datadir<-'N:/data/chl_phenology/data'
#figsdir<-'C:/Users/sailfish/Documents/aalldocuments/literature/postdoc_2013/chl_phenology/figures'
setwd(datadir)

#READ IN NAFO SHAPEFILE TO DEFINE SPATIAL BOUNDARIES
nafo<-readOGR('N:/data/dynamic_trophic_control/Assessments/data/fishing_shapefiles/nafo','Divisions')#works for Alaska - most others don't
nafo<-spTransform(nafo,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fourx<-subset(nafo,ZONE=='4X')

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/ACTIVE/spera/Figures')
pdf('div4x_map.pdf',height=4,width=4)
map('world',col='gray',fill=TRUE,lwd=.1,border=NA,xlim=c(-68,-63),ylim=c(42,46))
plot(fourx,col='royalblue',border=NA,add=TRUE)
axis(1,seq(-69,-63,1),las=1,cex.axis=.75,lwd=.1)
axis(2,seq(42,46,1),las=1,cex.axis=.75,lwd=.1)
box(lwd=.1)
dev.off()


cst<-readShapePoly('N:/data/shapefiles/naturalearthdata_ne_50m_ocean_poly/ne_50m_ocean.shp',proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cst<-crop(cst,extent(-70,-60,40,50))

#BOUNDING BOX FOR SUBSETTING
lon1<--69
lon2<--64
lat1<-42
lat2<-46

dat<-read.csv('N:/data/phytoplankton/functional_groups/alvain_dominance__monthly_9km_1997_2010_spera.csv',header=FALSE,col.names=c('lon','lat','proc','diat','slc','month','year'))
dat<-subset(dat,lat>=lat1 & lat<=lat2 & lon>=lon1 & lon<=lon2)

map('world',xlim=c(-70,-60),ylim=c(40,50))
map.axes()
points(dat$lon,dat$lat,pch=16,col='red')
points(dat$lon,dat$lat,pch=16,col='green')

#dat<-subset(dat,!(lon> -65.25 & lat<45) & lon< -62 & lat<45.5)
crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,cst)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)
pfg<-dat
setwd('N:/data/Spera/final')
write.csv(pfg,'phyto_functional_monthly_1997_2010_spera.csv',row.names=FALSE)

dat<-read.csv('N:/data/phytoplankton/size_groups/brewin_phyto_sizegroups_monthly_9km_1997_2007_spera.csv',header=FALSE,col.names=c('lon','lat','pic','nano','micro','month','year'))
dat<-subset(dat,lat>=lat1 & lat<=lat2 & lon>=lon1 & lon<=lon2)
#dat<-subset(dat,!(lon> -65.25 & lat<45) & lon< -62 & lat<45.5)
crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,cst)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)
psg<-dat
setwd('N:/data/Spera/final')
write.csv(psg,'phyto_size_monthly_1997_2007_spera.csv',row.names=FALSE)

####################################

####################################
#READS IN, BINS, FORMATS, REMOTE SENSING DATA
readfun<-function(nm){
print(nm)
dat<-read.csv(gsub(' ','',paste(nm,'.csv')),header=F,col.names=c('lon','lat','chl','year','day'))
gc()
if(nm=='czcs_chla_1day_9km_spera'){
dat<-subset(dat,year>=1979 & year<=1983)#UNCERTAIN OUTSIDE THIS RANGE
} else NULL
dat<-subset(dat,lat>=lat1 & lat<=lat2 & lon>=lon1 & lon<=lon2)
#dat<-subset(dat,!(lon> -65.25 & lat<45) & lon< -62 & lat<45.5)
crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,cst)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)
}
setwd(datadir)
md<-readfun('modis_chla_1day_4km_spera')
sw<-readfun('swfs_chla_1day_9km_spera')
mr<-readfun('meris_chla_1day_4km_spera')
cz<-readfun('czcs_chla_1day_9km_spera')


library(maps)
map('world',ylim=c(35,46),xlim=c(-70,-60))
points(md$lon,md$lat,pch=16,cex=.5,col='red')
a<-subset(md,year==2005)


#######################################################################

#######################################################################

olay<-function(dat){
dat<-subset(dat,lat>=lat1 & lat<=lat2 & lon>=lon1 & lon<=lon2)
#dat<-subset(dat,!(lon> -65.25 & lat<45) & lon< -62 & lat<45.5)
crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,cst)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)
return(dat)
}

#READ IN IN SITU DATA
setwd(datadir)
bdat<-read.csv('Boyce_etal_2012_integrated_CHL_1890_2010.csv',header=T,skip=30)
bdat<-subset(bdat,select=c('Longitude','Latitude','Chlc','Year','Dayofyear','Field'))
names(bdat)<-c('lon','lat','chl','year','day','field')
bdat$db<-rep('BOY',dim(bdat)[1])
bdat<-olay(bdat)

icdat<-read.csv('ICES_chl.csv',header=TRUE)
icdat<-subset(icdat,select=c('lon','lat','chl','year','doy'))
names(icdat)<-c('lon','lat','chl','year','day')
icdat$db<-rep('ICES',dim(icdat)[1])
icdat<-olay(icdat)

bcdat<-read.csv('BIOCHEM_chl.csv',header=TRUE)
bcdat<-subset(bcdat,select=c('lon','lat','chl','year','doy'))
names(bcdat)<-c('lon','lat','chl','year','day')
bcdat$db<-rep('BIOCHEM',dim(bcdat)[1])
bcdat<-olay(bcdat)



#REMOVE DUPLICATE ENTRIES FROM ALL IN SITU DATABASES; BUT NEED TO DEFINE ORDER IN WHICH TO KEEP/REMOVE OBSERVATIONS; TRUST FU THE LEAST SO DULPICATED FU'S REMOVED FIRST
idat<-rbind.fill(subset(bdat,field=='ZD'),
                 subset(bdat,field=='CHL'),
                 bcdat,
                 icdat,
                 subset(bdat,field=='FU'))

idat$lon<-round(idat$lon,digits=2)
idat$lat<-round(idat$lat,digits=2)
idat$chl<-round(idat$chl,digits=3)
idat<-unique(idat,by=c('lat','lon','year','day','chl'))#~3196 rows duplicated
idat<-subset(idat,lon>=-180)


#ADDS BATHYMETRY AND DISTANCE TO COAST
#crds<-unique(subset(data,select=c('lon','lat')))
#r2gmt(crds,datafile='phenology.coords.gmt')
b1<-read.table('phenology.bathy.txt',header=FALSE,col.names=c('lon','lat','x','bathy'),sep='\t')
b1<-b1[,-3]
b2<-read.table('phenology.bathy2.txt',header=FALSE,col.names=c('lon','lat','bathy'),sep='\t')
bathy<-unique(rbind(b1,b2))

dist<-unique(rbind(read.table('phenology.distance.txt',header=FALSE,col.names=c('lon','lat','dist','lond','latd'),sep='\t'),read.table('phenology.distance.txt',header=FALSE,col.names=c('lon','lat','dist','lond','latd'),sep='\t')))
idat<-merge(idat,bathy,by=c('lon','lat'),all.x=TRUE,all.y=FALSE)
idat<-merge(idat,dist,by=c('lon','lat'),all.x=TRUE,all.y=FALSE)


map('world',xlim=c(-70,-62),ylim=c(39,46))
map.axes()
points(idat$lon,idat$lat,pch=16,cex=.5,col='red')




#CPR DATA
cprdat2<-read.csv('NOAACPRDATA.csv',header=TRUE,skip=6,col.names=c('cruise','sample','year','month','mday','hour','minute','lat','lon','chl','marmap.taxcode','abundance','sciname'))
cprdat2$lon<-cprdat2$lon*-1
cprdat2<-subset(cprdat2,is.na(chl)==FALSE)
cprdat2$date<-strptime(gsub(' ','',paste(cprdat2$mday,'/',cprdat2$month,'/',cprdat2$year)),"%d/%m/%Y")
cprdat2$day<-yday(cprdat2$date)
cprdat2<-subset(cprdat2,select=c('lon','lat','chl','year','day'))

setwd(datadir)
cprdat<-read.csv('CPR.2015.csv',header=TRUE)
names(cprdat)<-tolower(names(cprdat))
names(cprdat)[3:7]<-c('time','ampm','lat','lon','chl')
cprdat$date<-strptime(cprdat$sampledate,"%d/%m/%Y")
cprdat$year<-as.numeric(substr(cprdat$date,1,4))
cprdat$day<-yday(cprdat$date)
cprdat<-subset(cprdat,select=c('lon','lat','chl','year','day'))

cprdat<-rbind(cprdat,cprdat2)
cprdat<-olay(cprdat)

map('world',xlim=c(-70,-62),ylim=c(39,46))
map.axes()
points(cprdat$lon,cprdat$lat,pch=16,cex=.5,col='red')




###########################################################
###########################################################
###########################################################
sw$db<-rep('SW',dim(sw)[1])
mr$db<-rep('ME',dim(mr)[1])
md$db<-rep('MO',dim(md)[1])
cz$db<-rep('CZ',dim(cz)[1])
cprdat$db<-rep('CPR',dim(cprdat)[1])
idat$db<-rep('INSITU',dim(idat)[1])

#COMBINES ALL DATA TOGETHER
data<-rbind.fill(cz,md,mr,cprdat,idat,sw)


map('world',xlim=c(-70,-62),ylim=c(39,46))
map.axes()
points(data$lon,data$lat,pch=16,cex=.5,col='red')

#ADDS DATE AND MONTH
data$date<-strptime(gsub(' ','',paste(data$year,'-',data$day)),"%Y-%j")
data$month<-month(data$date)



#BINS INTO 3 YEAR INTERVALS
brks<-seq(1886,2015,3); lbls<-seq(min(brks)+1.5,max(brks)-1.5,3)
data$tbin3<-cut(data$year,breaks=brks,labels=lbls)

#BINS INTO 5 YEAR INTERVALS
brks<-seq(1885,2015,5);lbls<-seq(min(brks)+2.5,max(brks)-2.5,5)
data$tbin5<-cut(data$year,breaks=brks,labels=lbls)

#BINS INTO 10 YEAR INTERVALS
brks<-seq(1885,2015,10);lbls<-seq(min(brks)+5,max(brks)-5,10)
data$tbin10<-cut(data$year,breaks=brks,labels=lbls)



#CALCULATES DAYS BETWEEN EACH OBSERVATION AND BASELINE TIME (1890/01/01); UNITS ARE IN SECONDS, SO DIVIDED BY 60(SECONDS IN MIN), 60(MINS IN HR), 24(HRS IN DAY) TO GET DAYS
data$tday<-int_length(new_interval(ymd('1890-01-01'),ymd(data$date)))/60/60/24

#HAVE TO CONVERT DATE FORMATE FOR DDPLY TO WORK - NOT SURE WHY
data$date<-as.POSIXct(data$date)
data$month<-month(data$date)

#ADDS DAY OF MONTH
data$mday<-mday(data$date)


data$cl<-ifelse(data$db=='INSITU','red',NA)
data$cl<-ifelse(data$db=='MO','lightblue',data$cl)
data$cl<-ifelse(data$db=='ME','blue',data$cl)
data$cl<-ifelse(data$db=='CZ','gray',data$cl)
data$cl<-ifelse(data$db=='SW','purple',data$cl)
data$cl<-ifelse(data$db=='CPR','green',data$cl)
plot(data$tday,data$chl)


setwd('N:/data/Spera/final')
save(data,file='phytoplankton_biomass_all_spera.RData')
