library(lubridate)
library(data.table)
library(psych)
library(reshape2)
library(gplots)
library(forecast)
library(cluster)
library(vegan)
library(ggplot2)
library(hybridHclust)
library(raster)
library(fields)
library(gridExtra)
library(colorRamps)
library(mapdata)
library(scales)
library(MASS)
library(mgcv)
library(maps)
library(plyr)
library(plotrix)
library(lubridate)
library(fossil)
library(rgdal)
library(rgeos)
library(maptools)

#datadir<-'N:/cluster_2017/scratch/chl_phenology/data'
datadir<-'N://cluster_2017//scratch//spera//data'
#figsdir<-'N:/cluster_2017/scratch/spera/figs'
codedir<-'N://cluster_2017//scratch//spera//code'

#setwd(codedir)
#source('helper_functions.r')

#READ IN SHAPEFILE TO DEFINE AREA OF INTEREST
setwd('N:/cluster_2017/scratch/spera/data/shapefile_70perc')
plg70<-readShapePoly('plgdf70.shp')#COMMAND TO READ BACK IN
proj4string(plg70)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

setwd('N:/cluster_2017/scratch/spera/data/shapefile_75perc')
plg75<-readShapePoly('plgdf75.shp')#COMMAND TO READ BACK IN
proj4string(plg75)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

setwd('N:/cluster_2017/scratch/spera/data/shapefile_75perc_v2')
plg75v2<-readShapePoly('plgdf75v2.shp')#COMMAND TO READ BACK IN
proj4string(plg75v2)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


#READ IN NAFO SHAPEFILE TO DEFINE SPATIAL BOUNDARIES
nafo<-readOGR('N:/data/dynamic_trophic_control/Assessments/data/fishing_shapefiles/nafo','Divisions')#works for Alaska - most others don't
nafo<-spTransform(nafo,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fourx<-subset(nafo,ZONE=='4X')
fourxw<-subset(nafo,ZONE%in% c('4X','4W'))


##################################################
#OVERLAY MAIN SPAWNING AREA ON DATA SETS
ofun<-function(d,plg){
print(dim(d))
pts<-SpatialPoints(data.frame(x=d$lon,y=d$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#OVERLAY LAND POLYGON
gint<-over(pts,plg)
#SUBSET POINTS THAT FALL WITHIN LAND POLYGONS
d$dum<-gint[,1]
d<-subset(d,is.na(dum)==FALSE)
d$dum<-NULL
print(dim(d))
print(names(d))
return(d)
}


##################################################
#ADDS 3,5,10 YEAR BINS
datefun<-function(d){
#BINS INTO 3 YEAR INTERVALS
brks<-seq(1900,2017,3); lbls<-seq(min(brks)+1.5,max(brks)-1.5,3)
d$tbin3<-cut(d$year,breaks=brks,labels=lbls)

#BINS INTO 5 YEAR INTERVALS
brks<-seq(1897,2017,5);lbls<-seq(min(brks)+2.5,max(brks)-2.5,5)
d$tbin5<-cut(d$year,breaks=brks,labels=lbls)

#BINS INTO 10 YEAR INTERVALS
brks<-seq(1887,2017,10);lbls<-seq(min(brks)+5,max(brks)-5,10)
d$tbin10<-cut(d$year,breaks=brks,labels=lbls)
return(d)
}



#####################################################################

######## IMPORTS, FORMATS, EXTRACTS OBSERVATIONS OVER AREA OF INTEREST,
######## AND WRITE TO FILE

#####################################################################
library(data.table)
library(h5)
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
file = h5file(name ='pathfinder_sst_daily_1981_2012.mat', mode = "r")
dataset = list.datasets(file)
dat<-readDataSet(file['a'])
gc()
dat<-t(dat)
dat<-data.frame(dat)
gc()
names(dat)<-c('lon','lat','sst','wind','year','month','day')
dat$lon<-round(dat$lon,digits=1)
dat$lat<-round(dat$lat,digits=1)
#[1] 114819480         7

f<-function(d){
    return(data.frame(lon=unique(d$lon),
                      lat=unique(d$lat),
                      day=unique(d$day),
                      month=unique(d$month),
                      year=unique(d$year),
                      sst=mean(d$sst),
                      wind=mean(d$wind)))
}
l<-dlply(dat,.(lon,lat,year,month,day),.fun=f,.progress='text')
gc()
dat<-rbindlist(l)
save(dat,file='pathfinder_sst_wind_v2.RData')
setwd('N:/cluster_2017/scratch/spera/data/finaldat')
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
load('pathfinder_sst_wind_v2.RData')
crds<-unique(subset(dat,select=c('lon','lat')))
map('world',xlim=c(-80,-55),ylim=c(40,50))
points(crds$lon,crds$lat,pch=16)

date<-strptime(gsub(' ','',paste(dat$month,'/',dat$day,'/',dat$year)),"%m/%d/%Y")
dat$day<-yday(date)
gc()
nm2<-'pathfinder_sst_daily_1981_2012_spawnar.RData'
dat<-ofun(dat,plg70)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
save(dat,file=paste(nm2))

dat<-ofun(dat,plg75)
setwd('N:/cluster_2017/scratch/spera/data/finaldat')
save(dat,file=paste(nm2))



#PHYTOPLANKTON BIOMSS
nm<-'phytoplankton_biomass_all_spera.RData'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
load(paste(nm))
dat<-data
dat<-dat[,!(names(dat) %in% c('date','dum'))]
print(dim(dat))
names(dat)<-tolower(names(dat))
dat<-datefun(dat)
#OVERLAY AREA OF INTEREST
nm2<-'phytoplankton_biomass_all_spera_spawnar.RData'
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg70)
save(dat,file=paste(nm2))
#dat<-ofun(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#save(dat,file=paste(nm2))


#PHYTOPLANKTON BIOMSS
nm<-'phytoplankton_biomass_all_spera_v2.RData'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
load(paste(nm))
dat<-data
dat<-dat[,!(names(dat) %in% c('date','dum'))]
print(dim(dat))
names(dat)<-tolower(names(dat))
dat<-datefun(dat)
#OVERLAY AREA OF INTEREST
nm2<-'phytoplankton_biomass_all_spera_allar.RData'
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,fourx)
save(dat,file=paste(nm2))


#BOF LARVAL SURVEY
nm<-'BoF_larval_herring_counts_spera.RData'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
#dat<-read.csv(paste(nm),header=TRUE)
load(paste(nm))
dat<-data.frame(data)
names(dat)<-tolower(names(dat))
dat$day<-dat$yday
dat<-dat[,!(names(dat) %in% c('drdglg','vessel','startdatetime','startlat','startlon','secchi','spec','tsn','yday'))]
print(dim(dat))
dat<-datefun(dat)
#OVERLAY AREA OF INTEREST
#nm2<-gsub('.csv','_spawnar.csv',nm)
dat<-ofun(dat,plg70)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
#write.csv(dat,paste(nm2),row.names=FALSE)
nm2<-gsub('.RData','_spawnar.RData',nm)
save(dat,file=paste(nm2))
#dat(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#write.csv(dat,paste(nm2),row.names=FALSE)


#BOF LARVAL SURVEY
nm<-'BoF_larval_herring_counts_spera.RData'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
#dat<-read.csv(paste(nm),header=TRUE)
load(paste(nm))
dat<-data.frame(data)
names(dat)<-tolower(names(dat))
dat$day<-dat$yday
dat<-dat[,!(names(dat) %in% c('drdglg','vessel','startdatetime','startlat','startlon','secchi','spec','tsn','yday'))]
print(dim(dat))
dat<-datefun(dat)
nm2<-gsub('.RData','_allar.RData',nm)
dat<-ofun(dat,fourxw)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
#write.csv(dat,paste(nm2),row.names=FALSE)
nm2<-gsub('.RData','_allar.RData',nm)
save(dat,file=paste(nm2))



#NITRATE
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
nm<-'Atlas1999nuts_nit.csv'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
dat$year<-ifelse(dat$year<1900,dat$year+1900,dat$year)
date<-strptime(gsub(' ','',paste(dat$month,'/',dat$day,'/',dat$year)),"%m/%d/%Y")
dat$day<-yday(date)
dat<-subset(dat,depth<=20)
avfun<-function(d){return(data.frame(nit=mean(d$nitrate,na.rm=TRUE)))}
dat<-ddply(dat,.(lon,lat,day,month,year),.fun=avfun,.progress='text')
#ADD BIOCHEM DATA
nm2<-'biochem_nitrate_spera.csv'
dat2<-read.csv(paste(nm2),header=TRUE)
names(dat2)<-tolower(names(dat2))
names(dat2)<-ifelse(names(dat2) %in% c('doy'),'day',names(dat2))
dat2<-subset(dat2,select=c('lon','lat','day','month','year','value'))
names(dat2)[6]<-'nit'
dat<-rbind(dat,dat2)
avfun<-function(d){return(data.frame(nit=mean(d$nit,na.rm=TRUE)))}
dat<-ddply(dat,.(lon,lat,day,month,year),.fun=avfun,.progress='text')

dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
nm2<-gsub('.csv','_spawnar.csv',nm)
dat<-ofun(dat,plg70)
write.csv(dat,paste(nm2),row.names=FALSE)
#dat<-ofun(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#write.csv(dat,paste(nm2),row.names=FALSE)



#SILICATE
nm<-'Atlas1999nuts_sil.csv'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
dat$date<-strptime(dat$date,"%m/%d/%Y")
dat$day<-yday(dat$date)
dat$year<-year(dat$date)
dat<-dat[,!(names(dat) %in% c('date'))]
dat<-subset(dat,depth<=20)
avfun<-function(d){return(data.frame(sil=mean(d$silicate,na.rm=TRUE)))}
dat<-ddply(dat,.(lon,lat,day,month,year),.fun=avfun,.progress='text')

#ADD BIOCHEM DATA
nm2<-'biochem_silicate_spera.csv'
dat2<-read.csv(paste(nm2),header=TRUE)
names(dat2)<-tolower(names(dat2))
names(dat2)<-ifelse(names(dat2) %in% c('doy'),'day',names(dat2))
dat2<-subset(dat2,select=c('lon','lat','day','month','year','value'))
names(dat2)[6]<-'sil'
dat<-rbind(dat,dat2)
avfun<-function(d){return(data.frame(sil=mean(d$sil,na.rm=TRUE)))}
dat<-ddply(dat,.(lon,lat,day,month,year),.fun=avfun,.progress='text')
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv','_spawnar.csv',nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg70)
write.csv(dat,paste(nm2),row.names=FALSE)
#dat<-ofun(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#write.csv(dat,paste(nm2),row.names=FALSE)




#PHOSPHATE
nm<-'Atlas1999nuts_phos.csv'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
names(dat)<-ifelse(names(dat)=='longitude','lon',names(dat))
names(dat)<-ifelse(names(dat)=='latitude','lat',names(dat))
dat$date<-strptime(dat$date,"%m/%d/%Y")
dat$day<-yday(dat$date)
dat$year<-year(dat$date)
dat$month<-month(dat$date)
dat<-dat[,!(names(dat) %in% c('date'))]
dat<-subset(dat,depth<=20)
avfun<-function(d){return(data.frame(phos=mean(d$phosphate,na.rm=TRUE)))}
dat<-ddply(dat,.(lon,lat,day,month,year),.fun=avfun,.progress='text')

#ADD BIOCHEM DATA
nm2<-'biochem_phosphate_spera.csv'
dat2<-read.csv(paste(nm2),header=TRUE)
names(dat2)<-tolower(names(dat2))
names(dat2)<-ifelse(names(dat2) %in% c('doy'),'day',names(dat2))
dat2<-subset(dat2,select=c('lon','lat','day','month','year','value'))
names(dat2)[6]<-'phos'
dat<-rbind(dat,dat2)
avfun<-function(d){return(data.frame(phos=mean(d$phos,na.rm=TRUE)))}
dat<-ddply(dat,.(lon,lat,day,month,year),.fun=avfun,.progress='text')
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv','_spawnar.csv',nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg70)
write.csv(dat,paste(nm2),row.names=FALSE)
#dat<-ofun(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#write.csv(dat,paste(nm2),row.names=FALSE)





#PHYTO FUNCTIONAL GROUPS AND SIZE CLASSES
rdfun<-function(nm){
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
dat$day<-(dat$month*30)-15
dat<-dat[,!(names(dat) %in% c('dum'))]
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv','_spawnar.csv',nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg70)
write.csv(dat,paste(nm2),row.names=FALSE)
dat<-ofun(dat,plg75)
setwd('N:/cluster_2017/scratch/spera/data/finaldat')
write.csv(dat,paste(nm2),row.names=FALSE)
}
rdfun('phyto_size_monthly_1997_2007_spera.csv')
rdfun('phyto_functional_monthly_1997_2010_spera.csv')

#CPR COUNTS FOR ZOOP AND PHYTO
rdfun<-function(nm,plg,nm2){
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
names(dat)<-ifelse(names(dat)=='doy','day',names(dat))
dat<-dat[,!(names(dat) %in% c('aphiaid','kingdom','dom','quarter'))]
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv',nm2,nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg)
write.csv(dat,paste(nm2),row.names=FALSE)
#dat<-ofun(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#write.csv(dat,paste(nm2),row.names=FALSE)
}
rdfun('phyto_cpr_counts_spera.csv',plg70,'_spawnar.csv')
rdfun('zoop_cpr_counts_spera.csv',plg70,'_spawnar.csv')
rdfun('phyto_cpr_counts_spera.csv',fourxw,'_allar.csv')
rdfun('zoop_cpr_counts_spera.csv',fourxw,'_allar.csv')


#SUMMER RV WEIGHTS AND LENGTHS; CREATES SPAWNING AREA OVERLAY AND NON
nm<-'herring_weights_RV_survey_spera.csv'
nm<-'herring_lengths_RV_survey_spera.csv'

rdfun<-function(nm){
setwd('N:/data/spera/final')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
dat$date<-strptime(dat$sdate,"%y-%m-%d")
dat$day<-yday(dat$date)
dat$year<-year(dat$date)
dat$month<-month(dat$date)
dat<-dat[,!(names(dat) %in% c('gear1','gear2','gear3','gear4','gear5','vname','bdate','sdate','keys.y','date','remarks','wind','market','remarks.x','entr','authority','tsn','comments','region','kingdom_id','kingdom','phylum','class','infraorder','nann','edge','chkmrk','remarks.y','specimen_id','keys.x','slat','slong','elat','elong','subspecies','fmat','howd','speed','dist','type','gear','aux','etime','vessel','season','edate','start_depth','end_depth','vesel'))]
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv','_spawnar.csv',nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,fourx)
write.csv(dat,paste(nm2),row.names=FALSE)
#dat<-ofun(dat,plg75)
#setwd('N:/cluster_2017/scratch/spera/data/finaldat')
#setwd("N:/cluster_2017/scratch/spera/data/finaldat_v2")
#write.csv(dat,paste(nm2),row.names=FALSE)
}
rdfun('herring_weights_RV_survey_spera.csv')
rdfun('herring_lengths_RV_survey_spera.csv')
rdfun('haddock_weights_RV_survey_spera.csv')
rdfun('haddock_lengths_RV_survey_spera.csv')


#SAME BUT NO OVERLAY- ENTIRE SURVEY DOMAIN
rdfun<-function(nm){
setwd('N:/data/spera/final')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
dat$date<-strptime(dat$sdate,"%y-%m-%d")
dat$day<-yday(dat$date)
dat$year<-year(dat$date)
dat$month<-month(dat$date)
dat<-dat[,!(names(dat) %in% c('gear1','gear2','gear3','gear4','gear5','vname','bdate','sdate','keys.y','date','remarks','wind','market','remarks.x','entr','authority','tsn','comments','region','kingdom_id','kingdom','phylum','class','infraorder','nann','edge','chkmrk','remarks.y','specimen_id','keys.x','slat','slong','elat','elong','subspecies','fmat','howd','speed','dist','type','gear','aux','etime','vessel','season','edate','start_depth','end_depth','vesel'))]
dat<-datefun(dat)
print(dim(dat))
nm2<-gsub('.csv','_allar.csv',nm)
setwd("N:/cluster_2017/scratch/spera/data/finaldat_v2")
#dat<-ofun(dat,fourxw)
write.csv(dat,paste(nm2),row.names=FALSE)
}
rdfun('herring_weights_RV_survey_spera.csv')
rdfun('herring_lengths_RV_survey_spera.csv')
rdfun('haddock_weights_RV_survey_spera.csv')
rdfun('haddock_lengths_RV_survey_spera.csv')


#STRATIFICATION, TEMPERATURE SALINIT PROFILES, SST
nm<-('physical_stratification_spera.csv')
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
date<-strptime(gsub(' ','',paste(dat$year,'/',dat$doy)),format='%Y/%j')
dat$month<-month(date)
names(dat)<-ifelse(names(dat)=='doy','day',names(dat))
dat<-dat[,!(names(dat) %in% c('id','doy','dom'))]
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv','_spawnar.csv',nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg70)
write.csv(dat,paste(nm2),row.names=FALSE)
dat<-ofun(dat,plg75)
setwd('N:/cluster_2017/scratch/spera/data/finaldat')
write.csv(dat,paste(nm2),row.names=FALSE)





#CPR RICHNESS FOR ZOOP AND PHYTO
f<-function(plg,lbl){
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
nm<-'plank_cpr_richness_spera.csv'
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
names(dat)<-ifelse(names(dat)=='doy','day',names(dat))
dat<-dat[,!(names(dat) %in% c('aphiaid','kingdom','dom','quarter'))]
dat<-datefun(dat)
print(dim(dat))
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv',lbl,nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg)
write.csv(dat,file=paste(nm2),row.names=FALSE)
}
f(plg70,'_spawnar.csv')
f(fourxw,'_allar.csv')



#CPR TURNOVER FOR ZOOP AND PHYTO
#IMPORTANT: TRIED ALSO CALCULATING USING GENUS OR FAMILY RATHER THAT SPECIES AND VERY SIMILAR
f<-function(plg,lbl){
nm<-'phyto_cpr_counts_spera.csv'
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
names(dat)<-ifelse(names(dat)=='doy','day',names(dat))
dat<-dat[,!(names(dat) %in% c('aphiaid','kingdom','dom','quarter'))]
#OVERLAY AREA OF INTEREST
nm2<-gsub('.csv',lbl,nm)
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
dat<-ofun(dat,plg)
write.csv(dat,paste(nm2),row.names=FALSE)
return(dat)
}
dat2<-f(plg70,'_spawnar.csv')
dat2<-f(fourxw,'_allar.csv')

##NEED TO RUN THIS SEPARATELY FOR SPAWNAR AND ALLAR ITERATIONS
#CALCULATES TURNOVER RATE BETWEEN SUCCESSIVE TIME POINTS
phyto2<-subset(dat2,month>=8 & counts>0)
dt<-sort(unique(phyto2$year))
l2<-list()
   for(j in 1:(length(dt)-1)){
       dt1<-subset(phyto2,year==dt[j])
       dt2<-subset(phyto2,year==dt[j+1])
       total<-length(unique(c(dt1$scientificname,dt2$scientificname)))
       #GETS SPECIES AT FIRST TIME THAT ARE NOT PRESENT IN SECOND
       emmigrants<-dim(subset(dt1,!(scientificname %in% dt2$scientificname)))[1]
       #GETS SPECIES AT FIRST TIME THAT ARE NOT PRESENT IN SECOND
       immigrants<-dim(subset(dt2,!(scientificname %in% dt1$scientificname)))[1]
       to<-(immigrants+emmigrants)/total
       yrs<-as.numeric(dt[j+1]-dt[j])
       a1<-subset(dt1,select=c('lon','lat','day','month'))
       a2<-subset(dt2,select=c('lon','lat','day','month'))
       prp<-(dim(a2)[1]/dim(a1)[1])*100

       l2[[j]]<-data.frame(year=dt[j+1],
                           myear=mean(c(dt[j],dt[j+1]),digits=1),
                           to=to,
                           deltayear=yrs,
                           turnover=to/yrs,
                           prp=prp,
                           n1=dim(a1)[1],
                           n2=dim(a2)[1],
                           nmonth=length(unique(a2$month)),
                           ntot=dim(a1)[1]+dim(a2)[1])
}
pto<-data.frame(do.call('rbind',l2))
pto<-subset(pto,n1>=5 & n2>=5 & deltayear<2)

setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
write.csv(pto,file='phyto_turnover_spera_allar.csv',row.names=FALSE)
write.csv(pto,file='phyto_turnover_spera_spawnar.csv',row.names=FALSE)









##########################################################
f<-function(nm){
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
#IMPORT AND FORMAT
dat<-read.csv(paste(nm),header=TRUE)
names(dat)<-tolower(names(dat))
names(dat)<-ifelse(names(dat)=='doy','day',names(dat))
dat<-dat[,!(names(dat) %in% c('aphiaid','kingdom','dom','quarter'))]
#OVERLAY AREA OF INTEREST
dat<-ofun(dat,fourxw)
dat<-ofun(dat,plg70)


#CALCULATES TURNOVER RATE BETWEEN SUCCESSIVE TIME POINTS
zoop2<-subset(dat,month>=8 & counts>0)
dt<-sort(unique(zoop2$year))
l2<-list()
   for(j in 1:(length(dt)-1)){
       dt1<-subset(zoop2,year==dt[j])
       dt2<-subset(zoop2,year==dt[j+1])
       total<-length(unique(c(dt1$scientificname,dt2$scientificname)))
       #GETS SPECIES AT FIRST TIME THAT ARE NOT PRESENT IN SECOND
       emmigrants<-dim(subset(dt1,!(scientificname %in% dt2$scientificname)))[1]
       #GETS SPECIES AT FIRST TIME THAT ARE NOT PRESENT IN SECOND
       immigrants<-dim(subset(dt2,!(scientificname %in% dt1$scientificname)))[1]
       to<-(immigrants+emmigrants)/total
       yrs<-as.numeric(dt[j+1]-dt[j])
       a1<-subset(dt1,select=c('lon','lat','day','month'))
       a2<-subset(dt2,select=c('lon','lat','day','month'))
       prp<-(dim(a2)[1]/dim(a1)[1])*100

       l2[[j]]<-data.frame(year=dt[j+1],
                           myear=mean(c(dt[j],dt[j+1]),digits=1),
                           to=to,
                           deltayear=yrs,
                           turnover=to/yrs,
                           prp=prp,
                           n1=dim(a1)[1],
                           n2=dim(a2)[1],
                           nmonth=length(unique(a2$month)),
                           ntot=dim(a1)[1]+dim(a2)[1])
}
zto<-data.frame(do.call('rbind',l2))
zto<-subset(zto,n1>=5 & n2>=5 & deltayear<2)
return(zto)
}

plot(zto$myear,zto$turnover,pch=15)
plot(pto$myear,pto$to,pch=15)


zto<-f('zoop_cpr_counts_spera_spawnar.csv')
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
write.csv(zto,file='zoop_turnover_spera_spawnar.csv',row.names=FALSE)
zto<-f('zoop_cpr_counts_spera_allar.csv')
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
write.csv(zto,file='zoop_turnover_spera_allar.csv',row.names=FALSE)

pto<-f('phyto_cpr_counts_spera_spawnar.csv')
setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
write.csv(pto,file='phyto_turnover_spera_spawnar.csv',row.names=FALSE)







###############################################################################


#ofun('herring_surv_abund_2015_spera.csv')
#ofun('herring_w_atage_2015_spera.csv')
#ofun('temperature_spera_BP.csv')

load("phytoplankton_biomass_all_spera_spawnar.RData")

#AVERAGE VALUES PER CELL
avefun<-function(d){
    return(data.frame(chl=mean(d$chl,na.rm=TRUE),
                      clond=unique(d$clond),
                      clatd=unique(d$clatd),
                      cell=unique(d$newarea5),
                      day=unique(d$day),
                      year=unique(d$year),
                      month=unique(d$month),
                      db=unique(d$db),
                      bathy=mean(d$bathy,na.rm=TRUE),
                      dist=mean(d$dist,na.rm=TRUE),
                      tday=unique(d$tday),
                      mday=unique(d$mday),
                      cl=unique(d$cl),
                      tbin3=unique(d$tbin3),
                      tbin5=unique(d$tbin5)))
}
l<-dlply(data,.(db,newarea5,year,day),.fun=avefun,.progress='text')
data<-rbindlist(l)

save(data,file='phytoplankton_biomass_all_spera_gridded.RData')
load('phytoplankton_biomass_all_spera_gridded.RData')
