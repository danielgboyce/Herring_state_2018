library(data.table)
library(raster)
library(sp)
library(maptools)
library(maps)
library(plyr)
library(lubridate)
library(rgdal)

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/fish_larvae/larv_herring_BOFsurv_DAN')
#setwd('C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/fish_larvae/larv_herring_BOFsurv_DAN')
bridgelogs<-read.table('bridgelogs.txt',header=TRUE,sep='\t')
latd<-as.numeric(substr(bridgelogs$STARTLAT,1,2))
latm<-round(as.numeric(substr(bridgelogs$STARTLAT,3,4))/60,digits=2)*100
lats<-as.numeric(substr(bridgelogs$STARTLAT,6,6))
lats<-ifelse(is.na(lats)==TRUE,0,lats)
bridgelogs$lat<-as.numeric(gsub(' ','',paste(latd,'.',latm,lats)))
lond<-as.numeric(substr(bridgelogs$STARTLON,1,2))
lonm<-round(as.numeric(substr(bridgelogs$STARTLON,3,4))/60,digits=2)*100
lons<-as.numeric(substr(bridgelogs$STARTLON,6,6))
lons<-ifelse(is.na(lons)==TRUE,0,lons)

bridgelogs$lon<-as.numeric(gsub(' ','',paste(lond,'.',lonm,lons)))*-1
date<-strptime(bridgelogs$STARTDATETIME,"%y-%m-%d")
bridgelogs$month<-month(date)
bridgelogs$day<-day(date)
bridgelogs$year<-year(date)
bridgelogs$yday<-yday(date)
bridgelogs$tday<-int_length(new_interval(ymd('1972-01-01'),ymd(date)))/60/60/24

samples<-read.table('samples.txt',header=TRUE,sep='\t')
catches<-read.table('catches.txt',header=TRUE,sep='\t')
lengths<-read.table('lengths.txt',header=TRUE,sep='\t')
stationloccodes<-read.table('stationloccodes.txt',header=TRUE,sep='\t')
species<-unique(read.csv('BoF_larval_herring_taxon_heirarchy.csv',header=TRUE))
names(species)[1]<-'SPEC'
species<-subset(species,sciname != 'CLUPEIDAE/OSMERIDAE F.')#2 ENTRIES FOR THIS ONE, SO DELETE ONE
fishproc<-read.table('fishproccodes.txt',header=FALSE,sep='\t',skip=1,col.names=c('FISHPROC','FISHPROCCODE'))#QUALITY INFO AOUT SAMPLE
surveytype<-read.table('surveytypecodes.txt',header=FALSE,sep='\t',skip=1,col.names=c('SURVEYTYPE','SURVEY'))#QUALITY INFO AOUT SAMPLE
invproccodes<-read.table('invproccodes.txt',header=FALSE,sep='\t',skip=1,col.names=c('INVPROC','INVPROCDESC'))#QUALITY INFO AOUT SAMPLE
#geartypecodes<-read.table('geartypecodes.txt',header=FALSE,sep='\t',skip=1,col.names=c('GEARTYPE','GEARDESC'))#QUALITY INFO AOUT SAMPLE
gearproccodes<-read.table('gearproccodes.txt',header=FALSE,sep='\t',skip=1,col.names=c('GEARPROC','GEARPROCDESC'))#QUALITY INFO AOUT SAMPLE

#dat<-merge(catches,samples,by=c("BRDGLG_ID","VESSEL","CRUISE","SETNO","SAMPLE"),all.x=TRUE,all.y=FALSE)
#dat<-subset(dat,is.na(VOLM3)==FALSE)

data<-merge(catches,samples,by=c("BRDGLG_ID","VESSEL","CRUISE","SETNO","SAMPLE"),all.x=TRUE,all.y=FALSE)
data<-merge(data,bridgelogs,by=c("BRDGLG_ID","VESSEL","CRUISE","SETNO"),all.x=TRUE,all.y=FALSE)
data<-merge(data,fishproc,by=c('FISHPROC'),all.x=TRUE,all.y=FALSE)
data<-merge(data,surveytype,by=c('SURVEYTYPE'),all.x=TRUE,all.y=FALSE)
data<-merge(data,species,by=c('SPEC'),all.x=TRUE,all.y=FALSE)
data<-merge(data,invproccodes,by=c('INVPROC'),all.x=TRUE,all.y=FALSE)
data<-merge(data,gearproccodes,by=c('GEARPROC'),all.x=TRUE,all.y=FALSE)
data<-subset(data,is.na(kingdom)==FALSE)
names(data)<-tolower(names(data))
data$rank<-tolower(data$rank)

data<-unique(data[,!(names(data) %in% c('stage1','stage2','stage3','stage4','deformed','flowserial','flowstart','flowend','flowdif','electronicflow','sinktime','netrateout','netratein','storage','jars','herring_index','observer_comments','rettime','hmlno sorter','subfishproc','subnoeggspec','subnofishspec','sorter_comments','cloudtype','cloudcover','shipspeed','shipheading','observer comments','btslide','hydrostn','gearproc','invproc','fishproc','hmlno','nofishspec','noinvspec','subinvproc','subnoinvspec','eggcura','fishcura','remcura','noeggspec','swellheight','sky','howlocobt','winddirection','swelldirection','windforce','barometer','invprocdesc','sorter','observer','comments'))])


#GETS COLLECTION OF ALL SAMPLES
sampdat<-unique(data[,!(names(data) %in% c('sciname','comname','aphiaid','rank','valid_name','kingdom','phylum','order','family','genus','tsn','species','spec','totno'))])



#d<-subset(data,spec==2001)
#FUNCTION ADDS ZEROES WHERE NECESSARY
zeroesfun<-function(d){
print(unique(d$spec))
d2<-merge(sampdat,d,by=c("surveytype", "brdglg_id", "vessel", "cruise",  "setno", "sample", "spec_stage", "netmesh", "towdur", "geardepth", "plankvol", "volm3","station", "startdatetime", "startlat","startlon","airtemp", "surftemp","bottemp", "startdepth", "enddepth","secchi", "geartype","lat", "lon", "month","day", "year", "yday", "tday", "fishproccode", "survey", "region", "gearprocdesc"),all.x=TRUE,all.y=FALSE)
d2$totno<-ifelse(is.na(d2$totno)==TRUE,0,d2$totno)
d2$sciname<-sort(unique(d$sciname))[1]
d2$comname<-sort(unique(d$comname))[1]
d2$spec<-unique(d$spec)
d2$aphiaid<-unique(d$aphiaid)
d2$rank<-unique(d$rank)
d2$valid_name<-unique(d$valid_name)
d2$kingdom<-unique(d$kingdom)
d2$phylum<-unique(d$phylum)
d2$order<-unique(d$order)
d2$family<-unique(d$family)
d2$genus<-unique(d$genus)
d2$tsn<-sort(unique(d$tsn))[1]
d2$species<-unique(d$species)
return(d2)
}
l<-dlply(data,.(spec),.fun=zeroesfun,.progress='text')
data<-rbindlist(l)



#################################################
#LENGTHS DATABASE
data2<-merge(lengths,catches,by=c("BRDGLG_ID","VESSEL","CRUISE","SETNO","SAMPLE",'SPEC_STAGE','SPEC'),all.x=TRUE,all.y=FALSE)
data2<-merge(data2,bridgelogs,by=c("BRDGLG_ID","VESSEL","CRUISE","SETNO"),all.x=TRUE,all.y=FALSE)
data2<-merge(data2,surveytype,by=c('SURVEYTYPE'),all.x=TRUE,all.y=FALSE)
data2<-merge(data2,species,by=c('SPEC'),all.x=TRUE,all.y=FALSE)
data2<-merge(data2,gearproccodes,by=c('GEARPROC'),all.x=TRUE,all.y=FALSE)
data2$rank<-tolower(data2$rank)
names(data2)<-tolower(names(data2))
data2<-subset(data2,is.na(rank)==FALSE)

(unique(data2$species))
unique(data2$comname)
summary(factor(data2$comname))

#READ IN NAFO SHAPEFILE TO DEFINE SPATIAL BOUNDARIES
nafo<-readOGR('N:/data/dynamic_trophic_control/Assessments/data/fishing_shapefiles/nafo','Divisions')#works for Alaska - most others don't
nafo<-spTransform(nafo,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fourx<-subset(nafo,ZONE=='4X')

olay<-function(dat){
crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,fourx)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)
return(dat)
}
#data<-olay(data)
#data2<-olay(data2)


#BOUNDING BOX FOR SUBSETTING
lon1<--69
lon2<--64
lat1<-42
lat2<-46

data<-subset(data,lat>=lat1 & lat<=lat2 & lon>=lon1 & lon<=lon2)
data2<-subset(data2,lat>=lat1 & lat<=lat2 & lon>=lon1 & lon<=lon2)


#OVERLAY LAND MASSES
cst<-readShapePoly('N:/data/shapefiles/naturalearthdata_ne_50m_ocean_poly/ne_50m_ocean.shp',proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cst<-crop(cst,extent(-70,-60,40,50))
proj4string(cst)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

crds<-SpatialPoints(data.frame(lon=data$lon,lat=data$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,cst)
data$dum<-gint[,1]
data<-subset(data,is.na(dum)==FALSE)
data<-data[,-49]

#ADD BATHYMETRY
library(gmt)
setwd('N:/data/Spera')
#r2gmt(subset(data,select=c('lon','lat')),datafile='crds.boflarv.gmt')
bathy<-read.table('boflarv.bathy.txt',header=FALSE,col.names=c('lon','lat','x','bathy'),sep='\t')
data$bathy<-bathy$bathy

setwd('N:/data/Spera/final')
#write.csv(data,'BoF_larval_herring_counts_spera.csv',row.names=FALSE)
save(data,file='BoF_larval_herring_counts_spera.RData')
save(data2,file='BoF_larval_herring_lengths_spera.RData')
write.csv(data2,'Bof_larval_herring_lengths_spera.csv',row.names=FALSE)


load('BoF_larval_herring_counts_spera.RData')






library(maps)
map('world',xlim=c(-69,-60),ylim=c(41,47))
map.axes()
points(data$lon,data$lat,pch=16,col='red',cex=.5)


setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/fish_larvae/larv_herring_BOFsurv_KEN')
data2<-read.csv('larv_herring_survey_data_spera.csv',header=TRUE)
#names(data)<-tolower(names(data))
names(data2)<-toupper(names(data2))

data2<-merge(data2,species,by='SPEC',all.x=TRUE,all.y=FALSE)
data2<-subset(data2,is.na(rank)==FALSE)
#data<-subset(data,is.na(common.name)==FALSE)#469 REMOVE MISSING NAMES (8401); UNABLE TO DETERMINE SPECIES


#ADD TAXONOMIC HEIRARCHIES
setwd('N:/data/dynamic_trophic_control/Surveys/ecnasap/species_codes_philgrayson/itis')
taxo<-read.csv('ecnasap_species_codes_taxo_heirarchies.csv')

taxo<-subset(taxo,region=='NS')
names(taxo)[2]<-c('spec')
taxo<-unique(taxo)
taxo$dup<-duplicated(taxo$spec)
taxo<-subset(taxo,dup==FALSE)
taxo<-subset(taxo,spec %in% unique(data$spec))
data2<-merge(data,taxo,by=c('spec'),all.x=TRUE,all.y=FALSE)


a<-subset(data2,is.na(kingdom)==TRUE)



aa<-subset(a,common.name=='WHITE BARRACUDINA')
t<-tapply(a$totnoperm3,a$common.name,sum)

data2<-merge(data,taxo,by=c('spec'),all.x=TRUE,all.y=FALSE)
a<-subset(data2,is.na(kingdom)==TRUE)




library(maps)
map('world',xlim=c(-69,-63),ylim=c(41,46))
map.axes()
points(a$lon,a$lat,pch=16,col='red',cex=.5)
points(data$lon,data$lat,pch=16,col='red',cex=.5)



library(rgdal)
library(raster)
nafo<-readOGR('N:/data/dynamic_trophic_control/Assessments/data/fishing_shapefiles/nafo','Divisions')#works for Alaska - most others don't
nafo<-spTransform(nafo,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(nafo,axes=TRUE,xlim=c(-70,-60),ylim=c(39,50))
fourx<-subset(nafo,ZONE=='4X')
fourx<-crop(fourx,extent(-66.9,-60,43,50))
plot(fourx,add=TRUE,col='red')

crds<-SpatialPoints(data.frame(lon=data$lon,lat=data$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


gint<-over(crds,fourx)
crds2<-data.frame(lon=data$lon,lat=data$lat,area=gint[,1])
crds2<-subset(crds2,is.na(area)==FALSE)
points(crds2$lon,crds2$lat,pch=16,col='red',cex=.5)



bathy200<-readOGR('N:/data/shapefiles/bathymetry_polygons/200m','ne_10m_bathymetry_K_200')#works for Alaska - most others don't
bathy200<-spTransform(bathy200,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
bathy1000<-readOGR('N:/data/shapefiles/bathymetry_polygons/1000m','ne_10m_bathymetry_J_1000')#works for Alaska - most others don't
bathy1000<-spTransform(bathy1000,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

abline(v=-67)
abline(v=-65)

library(scales)
map('world',xlim=c(-69,-50),ylim=c(35,50))
map.axes()
#plot(bathy200,add=TRUE,col='red')
plot(bathy1000,add=TRUE,col='purple')
plot(fourx,add=TRUE,col=alpha('blue',.25))

a<-subset(data,!(lon> -65 & lat<45) & lon< -62)
points(a$lon,a$lat,pch=16,col='red',cex=.5)
